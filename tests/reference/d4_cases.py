"""Extract pinned upstream D4/D4S tests without altering their values."""

import re
import sys
from pathlib import Path

source = Path(sys.argv[1]).read_text()
lines = iter(Path(sys.argv[2]).read_text().splitlines())
geometries = []
lattices = {}
for index in range(21):
    header = next(lines)
    geometries.append([header, *(next(lines) for _ in range(int(header.split()[0])))])
    if index >= 17:
        lattices[index] = next(lines)
assert next(lines, None) is None, "Unexpected exported geometry"
count = 0
for name, body in re.findall(
    r"\nsubroutine (test_\w+)\(error\)(.*?)\nend subroutine", source, re.S | re.I
):
    if not re.fullmatch(
        r"test_\w+_mb\d+|test_\w+_amf3|test_actinides_d4s?|test_smooth_cutoff", name
    ):
        continue
    record = re.search(r'get_structure\(mol, "MB16-43", "(\d+)"\)', body)
    if record:
        geometry = geometries[int(record[1]) - 1]
    else:
        numbers = re.search(r"num\(nat\)\s*=\s*\[(.*?)\]", body, re.S)
        coordinates = re.search(r"reshape\(\[(.*?)\]", body, re.S)
        assert numbers and coordinates, name
        numbers = re.findall(r"\d+", numbers[1])
        coordinates = re.findall(r"([+\-\d.Ee]+)_wp", coordinates[1])
        assert len(coordinates) == 3 * len(numbers), name
        geometry = [f"{len(numbers)} 0"]
        geometry.extend(
            " ".join([number, *coordinates[3 * atom : 3 * atom + 3]])
            for atom, number in enumerate(numbers)
        )
    constructor = re.search(r"param\s*=\s*rational_damping_param\((.*?)\)", body, re.S)
    values = dict(re.findall(r"(s6|s8|s9|a1|a2|alp)\s*=\s*([+\-\d.Ee]+)_wp", constructor[1]))
    assert len(values) == 6, name
    energy = re.search(r"call test_dftd4_gen\(.*?,\s*([+\-\d.Ee]+)_wp\)", body)
    if energy is None and "call test_dftd4_gen(" in body:
        energy = re.search(r":: ref\s*=\s*([+\-\d.Ee]+)_wp", body)
        assert energy, name
    kinds = []
    if energy:
        kinds.append("energy")
    if "call test_numgrad(" in body:
        kinds.append("gradient")
    if "call test_numsigma(" in body:
        kinds.append("sigma")
    if name == "test_smooth_cutoff":
        kinds.append("pairwise")
    assert kinds, name
    for kind in kinds:
        print(
            name, "d4s" if "type(d4s_model)" in body else "d4", kind, energy[1] if energy else "0"
        )
        print(*(values[key] for key in ("s6", "s8", "s9", "a1", "a2", "alp")))
        controls = {"ga": "3.0", "gc": "2.0", "wf": "6.0"}
        controls.update(re.findall(r"\b(ga|gc|wf)\s*=\s*([+\-\d.Ee]+)_wp", body))
        print(*(controls[key] for key in ("ga", "gc", "wf")))
        cutoff = {"cn": "30.0", "disp2": "60.0", "disp3": "40.0", "width2": "0.0", "width3": "0.0"}
        if name == "test_smooth_cutoff":
            smooth = re.search(r"smooth = realspace_cutoff\((.*?)\)", body, re.S)
            cutoff.update(
                re.findall(r"(cn|disp2|disp3|width2|width3)\s*=\s*([+\-\d.Ee]+)_wp", smooth[1])
            )
        print(*(cutoff[key] for key in ("cn", "disp2", "disp3", "width2", "width3")))
        print(*geometry, sep="\n")
    count += 1
assert count == 41, f"Expected 41 upstream cases, found {count}; update the coverage inventory"

parameters = Path(sys.argv[3]).read_text()
arrays = dict(re.findall(r":: (\w+)\(\*\)\s*=\s*\[(.*?)\]", parameters, re.S))
names = re.findall(r"'([^']+)'", arrays["names"])
aliases = re.findall(r"'([^']+)'", arrays["libxc_names"])
assert len(names) == len(aliases) == 67
for name, alias in zip(names, aliases, strict=True):
    print("alias", name, alias)
methods = re.findall(r"'([^']+)'", arrays["func"])
energies = re.findall(r"([+\-\d.Ee]+)_wp", arrays["ref"])
assert len(methods) == 118 and len(energies) == 119
for name, expected in zip(methods, energies[: len(methods)], strict=True):
    print(name, "d4", "named", expected)
    print(*geometries[16], sep="\n")

crystals = {
    name: index
    for index, name in enumerate(("ammonia", "acetic", "adaman", "anthracene"), start=17)
}
for filename, suite, expected_count in [
    (sys.argv[4], "pairwise", 6),
    (sys.argv[5], "periodic", 10),
]:
    source = Path(filename).read_text()
    cases = re.findall(
        r"\nsubroutine (test_\w+)\(error\)(.*?)\nend subroutine", source, re.S | re.I
    )
    registered = re.findall(r'new_unittest\("[^"]+",\s*(test_\w+)\)', source)
    assert len(cases) == len(registered) == expected_count, suite
    assert {name for name, _ in cases} == set(registered), suite
    assert re.search(r"thr\s*=\s*100\*epsilon\(1.0_wp\)", source)
    assert re.search(r"realspace_cutoff\(cn=30_wp, disp2=60.0_wp, disp3=15.0_wp\)", source)
    if suite == "pairwise":
        assert "sum(energy2) + sum(energy3), thr=thr" in source
    else:
        assert "thr2 = sqrt(epsilon(1.0_wp))" in source
        assert "thr3 = 100*sqrt(epsilon(1.0_wp))" in source
        assert re.findall(r"step = ([+\-\d.Ee]+)_wp", source) == ["1.0e-6", "1.0e-7"]
    for name, body in cases:
        record = re.search(r'get_structure\(mol, "(MB16-43|X23)", "(\w+)"\)', body)
        assert record, name
        periodic = record[1] == "X23"
        index = crystals[record[2]] if periodic else int(record[2]) - 1
        constructor = re.search(r"param\s*=\s*rational_damping_param\((.*?)\)", body, re.S)
        values = dict(re.findall(r"(s6|s8|s9|a1|a2|alp)\s*=\s*([+\-\d.Ee]+)_wp", constructor[1]))
        assert len(values) == 6, name
        expected = re.search(r"call test_dftd4_gen\(.*?,\s*([+\-\d.Ee]+)_wp\)", body)
        calls = re.findall(
            r"call (test_dftd4_pairwise|test_dftd4_gen|test_numgrad|test_numsigma)\(", body
        )
        assert len(calls) == 1, name
        kind = {
            "test_dftd4_pairwise": "upstream-pairwise",
            "test_dftd4_gen": "energy",
            "test_numgrad": "gradient",
            "test_numsigma": "sigma",
        }[calls[0]]
        assert (kind == "upstream-pairwise") == (suite == "pairwise"), name
        assert kind != "energy" or expected, name
        print(
            name,
            "d4s" if "type(d4s_model)" in body else "d4",
            ("periodic-" if periodic else "") + kind,
            expected[1] if expected else "0",
        )
        print(*(values[key] for key in ("s6", "s8", "s9", "a1", "a2", "alp")))
        print("3.0 2.0 6.0")
        print("30.0 60.0 15.0 0.0 0.0")
        print(*geometries[index], sep="\n")
        if periodic:
            print(lattices[index])
