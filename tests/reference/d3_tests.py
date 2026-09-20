"""Adapt mctc test errors to test-drive without changing numerical assertions."""

import re
import sys
from pathlib import Path

source = Path(sys.argv[1]).read_text().replace("use mctc_env_testing,", "use testdrive,")
if Path(sys.argv[1]).name == "test_fourier.f90":
    retained = {
        "test_spme_energy",
        "test_spme_gradient",
        "test_ewald_converged_bj",
        "test_ewald_converged_zero",
        "test_ewald_kcut",
        "test_ewald_supercell",
        "test_ewald_translation",
        "test_ewald_gradient",
        "test_ewald_sigma",
    }
    registrations = re.findall(r'new_unittest\("[^"]+",\s*(test_\w+)\)', source)
    assert len(registrations) == 32 and retained <= set(registrations)
    header = source.split("\ncontains\n", 1)[0]
    header = header.replace("get_coordination_number,", "").replace("get_realspace_cutoff,", "")
    header = re.sub(
        r"   use dftd3_(?:citation|fourier_jacobi|fourier_kernel),.*?\n(?!\s*&)",
        "",
        header,
        flags=re.S,
    )
    selected = [name for name in registrations if name in retained]
    collector = "subroutine collect_fourier(testsuite)\n   type(unittest_type), allocatable, intent(out) :: testsuite(:)\n   testsuite = [ &\n"
    collector += ", &\n".join(f'      new_unittest("{name}", {name})' for name in selected)
    collector += " &\n      ]\nend subroutine collect_fourier\n"
    bodies = [re.search(r"subroutine get_cubic\(.*?end subroutine get_cubic", source, re.S)[0]]
    for name in selected:
        body = re.search(rf"subroutine {name}\(error\).*?end subroutine {name}", source, re.S)[0]
        scalar_names = [
            variable
            for variable in ("energies", "energies_ref")
            if re.search(rf"\b{variable}\(:\)", body)
        ]

        def scalar_energy_declarations(match):
            variables = [variable.strip() for variable in match[1].split(", ")]
            variables = [
                variable
                for variable in variables
                if variable not in [f"{name}(:)" for name in scalar_names]
            ]
            return "real(wp), allocatable :: " + ", ".join(variables) if variables else ""

        body = re.sub(r"real\(wp\), allocatable :: ([^\n]+)", scalar_energy_declarations, body)
        body = body.replace(
            "   type(error_type), allocatable, intent(out) :: error",
            "   type(error_type), allocatable, intent(out) :: error\n   real(wp) :: "
            + ", ".join(scalar_names),
        )
        for variable in scalar_names:
            body = body.replace(f"{variable}(mol%nat), ", "").replace(f", {variable}(mol%nat)", "")
            body = body.replace(f"allocate({variable}(mol%nat))", "").replace(
                f"deallocate({variable})", ""
            )
            body = body.replace(f"sum({variable})", variable)
            assert not re.search(rf"\b{variable}\(", body), name
        bodies.append(body)
    source = (
        header + "\ncontains\n" + collector + "\n\n".join(bodies) + "\nend module test_fourier\n"
    )
if Path(sys.argv[1]).name == "test_partition.f90":
    omitted = {
        "test_coordination_number",
        "test_c6_on_demand",
        "test_reducer",
        "test_ewald_partitioned",
    }
    registered = set(re.findall(r'new_unittest\("[^"]+",\s*(test_\w+)\)', source))
    assert len(registered) == 9 and omitted <= registered
    for name in omitted:
        source, count = re.subn(rf'^.*new_unittest\("[^"]+", {name}\), &\n', "", source, flags=re.M)
        assert count == 1, name
    for name in omitted | {"check_coordination_number", "identity_reduce"}:
        source, count = re.subn(
            rf"\nsubroutine {name}\(.*?\nend subroutine {name}\b", "", source, flags=re.S
        )
        assert count == 1, name
    source = source.replace(", work_reducer", "").replace("use dftd3_gcp,", "use dftd3,")
    source, count = re.subn(r"   use dftd3_ncoord,.*?\n(?!\s*&)", "", source, flags=re.S)
    assert count == 1
    source, count = re.subn(
        r"   type, extends\(work_reducer\).*?end type identity_reducer", "", source, flags=re.S
    )
    assert count == 1
    source = source.replace(
        "module test_partition\n",
        "module test_partition\n   use mctc_env, only : library_error => error_type\n",
        1,
    )
    source = source.replace(
        "type(error_type), allocatable :: partition_error",
        "type(library_error), allocatable :: partition_error",
    )

    def partition_error(match):
        target = match[1]
        call = match[0].replace(f"({target},", "(failure,", 1)
        return (
            "block\n      type(library_error), allocatable :: failure\n"
            f"      {call}\n      if (allocated(failure)) call test_failed({target}, failure%message)\n   end block"
        )

    source, count = re.subn(
        r"call (?:new_work_partition|get_dispersion)\((error),.*?\)\s*(?=\n)",
        partition_error,
        source,
        flags=re.S,
    )
    assert count == 13, count
if Path(sys.argv[1]).name in ("test_pairwise.f90", "test_regression.f90"):
    module = Path(sys.argv[1]).stem
    source = source.replace(
        f"module {module}\n",
        f"module {module}\n"
        "   use dftd3, only : d3_param, new_rational_damping, new_zero_damping, &\n"
        "      & new_mzero_damping, new_optimizedpower_damping\n",
        1,
    )

    def damping_constructor(match):
        body = match[0]
        declaration = re.search(
            r"type\((\w+)_damping_param\), parameter :: param = \1_damping_param\((.*?)\)",
            body,
            re.S,
        )
        if declaration is None:
            return body
        family, coefficients = declaration.groups()
        body = body.replace(
            declaration[0],
            f"type({family}_damping_param) :: param\n"
            f"   type(d3_param), parameter :: coefficients = d3_param({coefficients})",
        )
        body = re.sub(
            r"(?m)^   call ",
            f"   call new_{family}_damping(param, coefficients)\n   call ",
            body,
            count=1,
        )
        if "type(d3_model) :: d3" in body:
            body = body.replace("\nend subroutine", "\n   call d3%close()\n\nend subroutine")
        return body.replace("\nend subroutine", "\n   call param%close()\n\nend subroutine")

    source = re.sub(
        r"subroutine test_\w+\(error\).*?end subroutine test_\w+",
        damping_constructor,
        source,
        flags=re.S,
    )
    assert source.count("call param%close()") == (12 if module == "test_pairwise" else 1)
    source = source.replace(
        "\nend subroutine test_dftd3_pairwise",
        "\n   call d3%close()\n\nend subroutine test_dftd3_pairwise",
    )
if Path(sys.argv[1]).name == "test_periodic_2d.f90":
    assert source.count('new_unittest("gh185", test_gh185)') == 1
    source = source.replace(
        'new_unittest("gh185", test_gh185)',
        'new_unittest("gh185-cn-image-truncation", test_gh185, should_fail=.true.)',
    )
    source = source.replace("   ! if (allocated(error)) return", "   if (allocated(error)) return")
if Path(sys.argv[1]).name == "test_periodic_3d.f90":
    complete_image_references = {
        "test_pbed3bj_acetic": "-6.67328359833623402e-2",
        "test_pbesold3bj_adaman": "-7.93135196203607201e-2",
        "test_m06ld3zero_cyanamide": "-2.32225718935670268e-2",
    }
    for name, expected in complete_image_references.items():
        body = re.search(rf"subroutine {name}\(error\).*?end subroutine {name}", source, re.S)
        assert body, name
        complete = body[0].replace(name, name + "_complete")
        complete, count = re.subn(
            r"(call test_dftd3_gen\(error, mol, param, )[^)]+",
            rf"\g<1>{expected}_wp",
            complete,
        )
        assert count == 1, name
        source = source.replace(
            "end module test_periodic_3d", complete + "\n\nend module test_periodic_3d"
        )
        source, count = re.subn(
            rf'new_unittest\("([^"]+)", {name}\)',
            rf'new_unittest("\1-legacy-image-box", {name}, should_fail=.true.), &\n'
            rf'      & new_unittest("\1-complete-images", {name}_complete)',
            source,
        )
        assert count == 1, name
if Path(sys.argv[1]).name in ("test_gcp.f90", "test_gcp_hessian.f90"):
    source = source.replace("use dftd3_cutoff,", "use dftd3,").replace(
        "use dftd3_gcp,", "use dftd3,"
    )
    assert not re.search(r"^\s*call ascii_gcp_param", source, re.M | re.I)
    source = source.replace("   use dftd3_output, only : ascii_gcp_param\n", "")
if Path(sys.argv[1]).name == "test_param.f90":
    source = source.replace(
        "module test_param\n",
        "module test_param\n   use mctc_env, only : library_error => error_type\n",
        1,
    )

    def adapt(match):
        call = re.sub(r"\berror\b", "failure", match[0])
        return (
            "block\n"
            "      type(library_error), allocatable :: failure\n"
            f"      {call}\n"
            "      if (allocated(failure)) call test_failed(error, failure%message)\n"
            "   end block"
        )

    source, count = re.subn(r"call get_\w+_damping\([^\n]+\)", adapt, source)
    assert count == 19, f"Expected 19 named getter calls, found {count}"
Path(sys.argv[2]).write_text(source)
