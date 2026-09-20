import math
import sys
from pathlib import Path

import disprs
from disprs import D3, D4, version

if "--bundled" in sys.argv:
    assert Path(disprs._lib._name).parent == Path(disprs.__file__).parent

assert version() == "0.1.0"
selection_numbers = [1, 28, 1]
selection_positions = [0, 0, 0, 3, 0.2, 0, 0.3, 1.4, 0.1]
for kind in (D3.D3, D3.D3S):
    with D3(selection_numbers, selection_positions, "blyp", model=kind) as selected:
        with D3(selection_numbers, selection_positions, "blyp", d3s=bool(kind)) as legacy:
            assert selected.dispersion(True) == legacy.dispersion(True)
            assert selected.properties() == legacy.properties()
for options in ({"model": 2}, {"model": D3.D3, "d3s": True}, {"model": D3.D3S, "d3s": False}):
    try:
        D3(selection_numbers, selection_positions, "blyp", **options)
    except ValueError:
        pass
    else:
        raise AssertionError("invalid/conflicting D3 model selector accepted")

for damping in (D3.ZERO, D3.RATIONAL):
    with D3(
        [1, 28, 1],
        [0, 0, 0, 3, 0.2, 0, 0.3, 1.4, 0.1],
        "blyp",
        damping=damping,
        atm=True,
        d3s=True,
    ) as smooth:
        original = smooth.dispersion(gradient=True)
        assert abs(sum(map(sum, smooth.pairwise())) - original[0]) < 1e-13
        assert all(math.isfinite(value) for value in smooth.hessian()[1])
        try:
            smooth.set_ewald(mesh=-1)
        except RuntimeError as error:
            assert "D3S Fourier" in str(error)
        else:
            raise AssertionError("D3S Fourier summation accepted")
        assert smooth.dispersion(gradient=True) == original
try:
    D3([95], [0, 0, 0], "pbe", d3s=True)
except RuntimeError as error:
    assert "1 through 94" in str(error)
else:
    raise AssertionError("unsupported D3S element accepted")

charge_numbers = [6, 8, 1]
charge_positions = [1.0, 1.0, 1.0, 3.0, 2.0, 1.0, 0.5, 3.0, 1.5]
charge_lattice = [12.0, 0.0, 0.0, 1.0, 13.0, 0.0, 0.5, 0.2, 14.0]
for model_type, kind in ((D3, D3.D3), (D3, D3.D3S), (D4, D4.D4), (D4, D4.D4S)):
    with model_type(charge_numbers, charge_positions, "pbe", model=kind, atm=True) as model:
        model.set_realspace_cutoff(
            coordination=8.0, dispersion2=9.0, dispersion3=7.0, width2=1.0, width3=1.0
        )
        energy, gradient, virial = model.dispersion(gradient=True)
        scalar = model.dispersion()
        assert abs(scalar[0] - energy) < 1e-12, (model_type.__name__, kind, scalar[0], energy)
        assert scalar[1:] == (None, None)
        assert len(gradient) == 9 and len(virial) == 9
        pair2, pair3 = model.pairwise()
        assert len(pair2) == len(pair3) == 9
        assert abs(sum(pair2) + sum(pair3) - energy) < 1e-12
        hessian_energy, hessian = model.hessian()
        assert abs(hessian_energy - energy) < 1e-12 and len(hessian) == 81
        for row in range(9):
            for column in range(9):
                assert abs(hessian[row + 9 * column] - hessian[column + 9 * row]) < 1e-10
        response = model.property_response()
        assert {"coordination", "c6"} <= response.keys()
        for name, count in (("coordination", 3), ("c6", 9)):
            assert tuple(map(len, response[name])) == (count, 9 * count, 9 * count)
        for indices in ([-1], [3]):
            try:
                model.set_ghosts(indices)
            except ValueError:
                pass
            else:
                raise AssertionError("invalid ghost index accepted")
        assert model.dispersion() == scalar
        partial = []
        for part in range(2):
            model.set_work_partition(part=part, parts=2)
            partial.append(model.dispersion(gradient=True))
        assert abs(sum(value[0] for value in partial) - energy) < 1e-12
        for output, full in ((1, gradient), (2, virial)):
            assert (
                max(
                    abs(left + right - total)
                    for left, right, total in zip(partial[0][output], partial[1][output], full)
                )
                < 1e-10
            )
        model.set_work_partition()
        model.update(charge_positions)
        assert model.dispersion() == scalar

for kind in (D4.D4, D4.D4S):
    with D4(
        charge_numbers,
        charge_positions,
        "pbe",
        model=kind,
        atm=False,
        lattice=charge_lattice,
        periodic=[True] * 3,
    ) as model:
        properties = model.properties()
        model.set_ewald(mesh=-1, tolerance=1e-10)
        direct = model.dispersion(gradient=True)
        for invalid in ({"tolerance": math.nan}, {"kcut": math.inf}, {"mesh": 2147483647}):
            try:
                model.set_ewald(**invalid)
            except RuntimeError:
                pass
            else:
                raise AssertionError("invalid Fourier controls accepted")
            assert model.dispersion(gradient=True) == direct
        model.set_ewald(mesh=32, tolerance=1e-10)
        mesh = model.dispersion(gradient=True)
        assert abs(mesh[0] - direct[0]) < 5e-8
        for actual, expected in zip(mesh[1:], direct[1:]):
            assert max(abs(left - right) for left, right in zip(actual, expected)) < 5e-8
        assert model.properties() == properties
        for unsupported in (model.hessian, model.pairwise):
            try:
                unsupported()
            except RuntimeError:
                pass
            else:
                raise AssertionError("unsupported Fourier output accepted")
        model.set_realspace_cutoff(width2=1.0)
        try:
            model.dispersion()
        except RuntimeError:
            pass
        else:
            raise AssertionError("Fourier pair smoothing accepted")
for dimensions, atm in ((0, False), (1, False), (2, False), (3, True)):
    with D4(
        charge_numbers,
        charge_positions,
        "pbe",
        atm=atm,
        lattice=charge_lattice,
        periodic=[axis < dimensions for axis in range(3)],
    ) as model:
        model.set_ewald(mesh=-1)
        try:
            model.dispersion()
        except RuntimeError:
            pass
        else:
            raise AssertionError("unsupported Fourier periodicity/ATM accepted")
for model_type, options in [
    (D3, {}),
    (D3, {"d3s": True}),
    (D4, {}),
    (D4, {"model": D4.D4S}),
    (D4, {"charge_model": "eeqbc"}),
    (D4, {"ga": 2.0, "gc": 1.0, "wf": 4.0, "charge": 0.25}),
    (D4, {"model": D4.D4S, "ga": 2.0, "gc": 1.0, "charge_model": "eeqbc", "charge": 0.25}),
    (D4, {"fixed_charges": [0.2, -0.4, 0.1]}),
    (D4, {"model": D4.D4S, "fixed_charges": [0.2, -0.4, 0.1], "charge_model": "eeqbc"}),
]:
    for dimensions in range(4):
        with model_type(
            charge_numbers,
            charge_positions,
            "pbe",
            lattice=charge_lattice,
            periodic=[axis < dimensions for axis in range(3)],
            **options,
        ) as model:
            model.set_realspace_cutoff(coordination=7.0)
            responses = model.property_response()
            if "fixed_charges" in options:
                assert responses["charges"][0] == options["fixed_charges"]
                assert responses["charges"][1:] == ([0.0] * 27, [0.0] * 27)
            for expected, (values, cartesian, strain) in zip(
                model.properties(), responses.values()
            ):
                assert (
                    max(abs(actual - reference) for actual, reference in zip(values, expected))
                    < 1e-10
                )
                assert len(cartesian) == 9 * len(values) and len(strain) == 9 * len(values)
                for prop in range(len(values)):
                    for axis in range(3):
                        assert (
                            abs(sum(cartesian[9 * prop + 3 * atom + axis] for atom in range(3)))
                            < 1e-10
                        )
            model.set_ghosts([1])
            model.set_work_partition(1, 3)
            assert model.property_response() == responses
            step = 1e-5
            for component in range(18):
                changed = []
                for sign in (1, -1):
                    xyz = charge_positions.copy()
                    cell = charge_lattice.copy()
                    if component < 9:
                        xyz[component] += sign * step
                    else:
                        row, column = (component - 9) % 3, (component - 9) // 3
                        for atom in range(3):
                            xyz[3 * atom + row] += sign * step * charge_positions[3 * atom + column]
                            cell[3 * atom + row] += sign * step * charge_lattice[3 * atom + column]
                    model.update(xyz, lattice=cell)
                    changed.append(model.properties())
                for response, plus, minus in zip(responses.values(), *changed):
                    derivative = response[1 if component < 9 else 2]
                    for prop, (left, right) in enumerate(zip(plus, minus)):
                        assert (
                            abs(derivative[9 * prop + component % 9] - (left - right) / (2 * step))
                            < 3e-5
                        )
            model.update(charge_positions, lattice=charge_lattice)
            model.set_realspace_cutoff(coordination=0.0)
            assert all(
                value == 0.0
                for array in model.property_response()["coordination"][1:]
                for value in array
            )
            if model_type is D4:
                for derivatives in responses["charges"][1:]:
                    for component in range(9):
                        assert (
                            abs(sum(derivatives[9 * atom + component] for atom in range(3))) < 1e-10
                        )
        try:
            model.property_response()
        except RuntimeError:
            pass
        else:
            raise AssertionError("closed property response accepted")

for kind in (D4.D4, D4.D4S):
    for charge_model in ("eeq", "eeqbc"):
        for dimensions in range(4):
            with D4(
                charge_numbers,
                charge_positions,
                "pbe",
                model=kind,
                charge_model=charge_model,
                lattice=charge_lattice,
                periodic=[axis < dimensions for axis in range(3)],
            ) as model:
                model.set_realspace_cutoff(coordination=7.0, dispersion2=8.0, dispersion3=6.0)
                reference = model.properties()
                original_energy = model.dispersion()[0]
                model.set_fixed_charges(reference[1])
                assert model.properties() == reference
                assert abs(model.dispersion()[0] - original_energy) < 1e-13
                supplied = [0.2, -0.4, 0.1]
                model.set_fixed_charges(supplied)
                supplied[0] = 999.0
                for invalid in ([0.0], [], [math.nan, 0.0, 0.0], [math.inf, 0.0, 0.0]):
                    try:
                        model.set_fixed_charges(invalid)
                    except (ValueError, RuntimeError):
                        pass
                    else:
                        raise AssertionError("invalid fixed charges accepted")
                    assert model.properties()[1] == [0.2, -0.4, 0.1]
                energy, gradient, virial = model.dispersion(gradient=True)
                assert abs(sum(map(sum, model.pairwise())) - energy) < 1e-13
                step = 1e-5
                for component in range(18):
                    energies = []
                    for sign in (1, -1):
                        xyz, cell = charge_positions.copy(), charge_lattice.copy()
                        if component < 9:
                            xyz[component] += sign * step
                        else:
                            row, column = (component - 9) % 3, (component - 9) // 3
                            for atom in range(3):
                                xyz[3 * atom + row] += (
                                    sign * step * charge_positions[3 * atom + column]
                                )
                                cell[3 * atom + row] += (
                                    sign * step * charge_lattice[3 * atom + column]
                                )
                        model.update(xyz, lattice=cell)
                        energies.append(model.dispersion()[0])
                        assert model.properties()[1] == [0.2, -0.4, 0.1]
                    derivative = gradient[component] if component < 9 else virial[component - 9]
                    assert abs((energies[0] - energies[1]) / (2 * step) - derivative) < 2e-8
                model.update(charge_positions, lattice=charge_lattice)
                _, hessian = model.hessian()
                moved = charge_positions.copy()
                moved[3] += 1e-4
                model.update(moved)
                plus = model.dispersion(gradient=True)[1]
                moved[3] -= 2e-4
                model.update(moved)
                minus = model.dispersion(gradient=True)[1]
                assert (
                    max(
                        abs((left - right) / 2e-4 - hessian[27 + index])
                        for index, (left, right) in enumerate(zip(plus, minus))
                    )
                    < 2e-8
                )
                model.update(charge_positions)
                model.set_fixed_charges()
                assert model.properties() == reference

for dimensions in range(4):
    with D3(
        charge_numbers,
        charge_positions,
        "pbe",
        lattice=charge_lattice,
        periodic=[axis < dimensions for axis in range(3)],
    ) as model:
        model.set_realspace_cutoff(coordination=7.0)
        reference_cn, reference_c6 = model.properties()
        assert len(reference_cn) == 3 and len(reference_c6) == 9
        assert all(math.isfinite(value) and value > 0 for value in reference_cn + reference_c6)
        assert all(
            reference_c6[3 * first + second] == reference_c6[3 * second + first]
            for first in range(3)
            for second in range(3)
        )
        model.set_ghosts([1])
        model.set_work_partition(1, 3)
        assert model.properties() == (reference_cn, reference_c6)
        moved = charge_positions.copy()
        moved[3] += 0.2
        model.update(moved)
        assert model.properties()[0] != reference_cn
        model.update(charge_positions)
        assert model.properties() == (reference_cn, reference_c6)
        if dimensions:
            moved = charge_positions.copy()
            for axis in range(3):
                moved[axis] += 3 * charge_lattice[axis]
            model.update(moved)
            translated_cn, translated_c6 = model.properties()
            assert (
                max(
                    abs(actual - expected)
                    for actual, expected in zip(
                        translated_cn + translated_c6, reference_cn + reference_c6
                    )
                )
                < 1e-11
            )
            model.update(charge_positions)
        model.set_realspace_cutoff(coordination=0.0)
        assert model.properties()[0] == [0.0] * 3
        assert model.properties()[1] != reference_c6
    try:
        model.properties()
    except RuntimeError:
        pass
    else:
        raise AssertionError("closed D3 property query accepted")

for charge_model in ("eeq", "eeqbc"):
    for dimensions in range(4):
        options = dict(
            charge=0.25,
            charge_model=charge_model,
            lattice=charge_lattice,
            periodic=[axis < dimensions for axis in range(3)],
        )
        charges, cartesian, strain = disprs.get_charges(
            charge_numbers, charge_positions, cartesian=True, strain=True, **options
        )
        assert abs(sum(charges) - 0.25) < 1e-12
        assert len(cartesian) == len(strain) == 27
        assert all(math.isfinite(value) for value in [*charges, *cartesian, *strain])
        with D4(charge_numbers, charge_positions, "pbe", **options) as model:
            assert (
                max(
                    abs(actual - expected)
                    for actual, expected in zip(charges, model.properties()[1])
                )
                < 1e-12
            )
        assert disprs.get_charges(charge_numbers, charge_positions, **options) == (
            charges,
            None,
            None,
        )
        assert disprs.get_charges(charge_numbers, charge_positions, cartesian=True, **options) == (
            charges,
            cartesian,
            None,
        )
        assert disprs.get_charges(charge_numbers, charge_positions, strain=True, **options) == (
            charges,
            None,
            strain,
        )
        for component in range(9):
            assert abs(sum(cartesian[9 * atom + component] for atom in range(3))) < 1e-11
            assert abs(sum(strain[9 * atom + component] for atom in range(3))) < 1e-11
        shifted = charge_positions.copy()
        shifted[4] += 1e-5
        plus = disprs.get_charges(charge_numbers, shifted, **options)[0]
        shifted[4] -= 2e-5
        minus = disprs.get_charges(charge_numbers, shifted, **options)[0]
        assert all(
            abs(cartesian[9 * atom + 4] - (plus[atom] - minus[atom]) / 2e-5) < 1e-7
            for atom in range(3)
        )
        strained_charges = []
        for step in (1e-5, -1e-5):
            shifted = charge_positions.copy()
            cell = charge_lattice.copy()
            for atom in range(3):
                shifted[3 * atom + 1] += step * charge_positions[3 * atom]
                cell[3 * atom + 1] += step * charge_lattice[3 * atom]
            strained_charges.append(
                disprs.get_charges(charge_numbers, shifted, **{**options, "lattice": cell})[0]
            )
        assert all(
            abs(
                strain[9 * atom + 1]
                - (strained_charges[0][atom] - strained_charges[1][atom]) / 2e-5
            )
            < 1e-7
            for atom in range(3)
        )
for options in ({"charge_model": "unknown"}, {"charge": float("nan")}, {"positions": [0.0]}):
    try:
        disprs.get_charges(**{"numbers": charge_numbers, "positions": charge_positions, **options})
    except (ValueError, RuntimeError):
        pass
    else:
        raise AssertionError("invalid charge query accepted")

for element in (112, 118):
    with D4([element, 8], charge_positions[:6], "pbe", charge=0.25) as model:
        model.update(charge_positions[:6])
        assert math.isfinite(model.dispersion()[0])
        actual = disprs.get_charges([element, 8], charge_positions[:6], charge=0.25)[0]
        assert (
            max(abs(value - expected) for value, expected in zip(actual, model.properties()[1]))
            < 1e-12
        )
    try:
        disprs.get_charges([element, 8], charge_positions[:6], charge_model="eeqbc")
    except RuntimeError:
        pass
    else:
        raise AssertionError("unsupported EEQBC element accepted")

selection_numbers = [6, 8, 7, 1]
selection_positions = [0.2, 0.3, 0.4, 2.6, 0.7, 0.8, 1.2, 2.8, 0.6, 0.7, 1.1, 2.3]
for factory, options in [(D3, {})] + [
    (D4, dict(model=kind, charge_model=charge_model, charge=0.5))
    for kind in (D4.D4, D4.D4S)
    for charge_model in ("eeq", "eeqbc")
]:
    for dimensions in range(4):
        with factory(
            selection_numbers,
            selection_positions,
            "pbe",
            atm=True,
            ghosts=[3],
            lattice=[5.3, 0.0, 0.0, 0.6, 5.7, 0.0, 0.3, 0.4, 6.2],
            periodic=[axis < dimensions for axis in range(3)],
            **options,
        ) as model:
            model.set_realspace_cutoff(7.0, 6.0, 6.0, width2=1.0, width3=1.0)
            if factory is D4:
                model.set_charge_cutoff(10.0)
                properties = model.properties()
            energy, gradient, virial = model.dispersion(gradient=True)
            if factory is D4:
                hessian_energy, hessian = model.hessian()
                assert abs(hessian_energy - energy) < 1e-13
                assert len(hessian) == 144 and all(math.isfinite(value) for value in hessian)
                assert all(
                    abs(hessian[column * 12 + row] - hessian[row * 12 + column]) < 1e-12
                    for row in range(12)
                    for column in range(12)
                )
                assert all(
                    abs(sum(hessian[(atom * 3 + axis) * 12 + row] for atom in range(4))) < 1e-12
                    for row in range(12)
                    for axis in range(3)
                )
                total_hessian = [0.0] * len(hessian)
                shifted = selection_positions.copy()
                shifted[4] += 2e-4
                model.update(shifted)
                plus_gradient = model.dispersion(gradient=True)[1]
                shifted[4] -= 4e-4
                model.update(shifted)
                minus_gradient = model.dispersion(gradient=True)[1]
                model.update(selection_positions)
                assert all(
                    abs(hessian[12 * 4 + row] - (plus_gradient[row] - minus_gradient[row]) / 4e-4)
                    < 2e-7
                    for row in range(12)
                )
            pairs = model.pairwise()
            expected = [energy, *gradient, *virial, *pairs[0], *pairs[1]]
            totals = [0.0] * len(expected)
            for part in range(3):
                model.set_work_partition(part, 3)
                value, grad, stress = model.dispersion(gradient=True)
                pair2, pair3 = model.pairwise()
                assert abs(value - sum(pair2) - sum(pair3)) < 1e-13
                assert all(
                    matrix[3 * 4 + atom] == matrix[atom * 4 + 3] == 0.0
                    for matrix in (pair2, pair3)
                    for atom in range(4)
                )
                actual = [value, *grad, *stress, *pair2, *pair3]
                totals = [total + item for total, item in zip(totals, actual, strict=True)]
                if factory is D4:
                    assert model.properties() == properties
                    _, local_hessian = model.hessian()
                    total_hessian = [
                        total + value for total, value in zip(total_hessian, local_hessian)
                    ]
            assert all(
                abs(total - value) < 1e-12 for total, value in zip(totals, expected, strict=True)
            )
            if factory is D4:
                assert (
                    max(abs(total - value) for total, value in zip(total_hessian, hessian)) < 1e-9
                )
            model.set_work_partition()
            model.update(selection_positions)
            assert abs(model.dispersion()[0] - energy) < 1e-13
            model.set_ghosts([0, 1, 2])
            value, grad, stress = model.dispersion(gradient=True)
            assert value == 0.0 and all(item == 0.0 for item in [*grad, *stress])

numbers = [6, 6]
positions = [0.0, 0.0, 0.0, 6.0, 0.0, 0.0]

with D3(numbers, positions, "pbe") as model:
    energy, gradient, _ = model.dispersion(gradient=True)
    assert math.isfinite(energy) and energy < 0.0
    assert abs(gradient[0] + gradient[3]) < 1e-12
    energy, gradient, virial = model.counterpoise("pbeh3c", gradient=True)
    assert abs(energy - 2.581420786404301e-5) < 1e-15
    assert abs(gradient[0] - 5.324722896549186e-5) < 1e-14
    assert abs(virial[0] + 3.1948337379295114e-4) < 1e-13
    energy, hessian = model.counterpoise_hessian("pbeh3c")
    assert abs(energy - 2.581420786404301e-5) < 1e-15
    assert abs(hessian[0] - 9.51416946673085e-5) < 1e-13

with D3(numbers, positions, "pbe", ghosts=[0]) as model:
    energy, gradient, _ = model.dispersion(gradient=True)
    assert energy == 0.0
    assert all(value == 0.0 for value in gradient)

with D3(numbers, [0.0, 0.0, 0.0, 5.5, 0.0, 0.0], "pbe") as model:
    model.set_realspace_cutoff(dispersion2=6.0, width2=2.0)
    energy, gradient, _ = model.dispersion(gradient=True)
    assert abs(energy - -6.405046062946415e-5) < 1e-15
    assert abs(gradient[0] - -3.4189782676699354e-4) < 1e-14

with D3(
    [6, 8],
    [0.4, 0.8, 1.2, 4.4, 2.8, 3.6],
    "pbe",
    lattice=[8.0, 0.0, 0.0, 1.0, 9.0, 0.0, 0.5, 0.7, 10.0],
    periodic=[True, True, True],
) as model:
    energy, gradient, virial = model.counterpoise(
        "pbeh3c", gradient=True, cutoff=12.0, srb_cutoff=12.0
    )
    assert abs(energy - 0.000733389638788683) < 1e-13
    assert abs(gradient[0] - 6.214839320259188e-7) < 1e-13
    assert abs(virial[0] + 0.003486778083367054) < 1e-12
    model.update(
        [0.4, 0.8, 1.2, 4.4, 2.8, 3.6],
        lattice=[8.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 8.0],
    )
    energy, _, _ = model.counterpoise("pbeh3c", cutoff=12.0, srb_cutoff=12.0)
    assert abs(energy - 0.0007503293440007774) < 1e-13
    energy, hessian = model.counterpoise_hessian("pbeh3c", cutoff=12.0, srb_cutoff=12.0)
    assert abs(energy - 0.0007503293440007774) < 1e-13
    assert abs(hessian[0] - 0.0006505276577134137) < 1e-12
    assert abs(hessian[3] + 0.0006505276577134137) < 1e-12

with D4(numbers, positions, "pbe", atm=False) as model:
    energy, gradient, _ = model.dispersion(gradient=True)
    assert math.isfinite(energy) and energy < 0.0
    assert abs(gradient[0] + gradient[3]) < 1e-12
    coordination, charges, c6, polarizabilities = model.properties()
    assert len(coordination) == len(charges) == len(polarizabilities) == 2
    assert len(c6) == 4

with D3(numbers, positions, "pbe") as model:
    assert model.counterpoise("b973c", srb=False)[0] == 0.0
    assert model.counterpoise("b973c", srb=True)[0] != 0.0
    assert model.counterpoise("pbeh3c", eta=0.0, base=False, srb=False)[0] == 0.0
    value, gradient, _ = model.counterpoise("r2scan3c", eta=1.1, base=False, gradient=True)
    assert abs(value - model.counterpoise_hessian("r2scan3c", eta=1.1, base=False)[0]) < 1e-15
    step = 1e-5
    model.update([*positions[:3], positions[3] + step, *positions[4:]])
    plus = model.counterpoise("r2scan3c", eta=1.1, base=False)[0]
    model.update([*positions[:3], positions[3] - step, *positions[4:]])
    minus = model.counterpoise("r2scan3c", eta=1.1, base=False)[0]
    assert abs(gradient[3] - (plus - minus) / (2 * step)) < 1e-10
    model.update(positions)
    values = model.counterpoise_parameters("pbeh3c")
    assert values["zeff"] == numbers and len(values["rvdw"]) == 4
    original = model.counterpoise("pbeh3c")[0]
    assert abs(model.counterpoise("pbeh3c", parameters=values)[0] - original) < 1e-15
    assert (
        abs(
            model.counterpoise("pbeh3c", parameters={"sigma": 2 * values["sigma"]})[0]
            - 2 * original
        )
        < 1e-15
    )
    custom = {"dmp_scal": 2.0, "dmp_exp": 4.0, "slater": [1.5, 1.8], "rvdw": [5.0] * 4}
    value, gradient, _ = model.counterpoise("pbeh3c", parameters=custom, gradient=True)
    hessian_value, hessian = model.counterpoise_hessian("pbeh3c", parameters=custom)
    assert abs(value - hessian_value) < 1e-15
    model.update([*positions[:3], positions[3] + step, *positions[4:]])
    plus, plus_gradient, _ = model.counterpoise("pbeh3c", parameters=custom, gradient=True)
    model.update([*positions[:3], positions[3] - step, *positions[4:]])
    minus, minus_gradient, _ = model.counterpoise("pbeh3c", parameters=custom, gradient=True)
    assert abs(gradient[3] - (plus - minus) / (2 * step)) < 1e-10
    assert abs(hessian[3 + 6 * 3] - (plus_gradient[3] - minus_gradient[3]) / (2 * step)) < 1e-9

for kind in (D4.D4, D4.D4S):
    with D4(
        numbers, positions, "pbe", model=kind, ga=2.0, gc=1.0, wf=5.0 if kind == D4.D4 else 6.0
    ) as model:
        original = model.dispersion()[0]
        model.set_realspace_cutoff(dispersion2=8.0, dispersion3=8.0, width2=4.0, width3=4.0)
        actual = model.dispersion(gradient=True)[0]
        assert abs(actual - original) > 1e-10
        pair2, pair3 = model.pairwise()
        assert abs(sum(pair2) + sum(pair3) - actual) < 1e-14

water_numbers = [8, 1, 1]
water_positions = [0.0, 0.0, 0.0, 0.0, 1.4, 1.0, 0.0, -1.4, 1.0]
for periodic in ([True, False, False], [True, False, True], [True] * 3):
    lattice = [9.0, 0.2, 0.0, 0.4, 10.0, 0.1, 0.2, 0.3, 11.0]
    for api in (D3, D4):
        with api(
            water_numbers, water_positions, "pbe", lattice=lattice, periodic=periodic
        ) as model:
            model.set_realspace_cutoff(dispersion2=12.0, dispersion3=9.0, coordination=10.0)
            if api is D4:
                model.set_charge_cutoff(80.0)
            energy, gradient, _ = model.dispersion(gradient=True)
            pair2, pair3 = model.pairwise()
            assert abs(sum(pair2) + sum(pair3) - energy) < 1e-13
            if api is D3:
                hessian_energy, hessian = model.hessian()
                assert abs(hessian_energy - energy) < 1e-13
                shifted = water_positions.copy()
                shifted[4] += 1e-5
                model.update(shifted)
                plus = model.dispersion(gradient=True)[1]
                shifted[4] -= 2e-5
                model.update(shifted)
                minus = model.dispersion(gradient=True)[1]
                for row in range(9):
                    assert abs(hessian[9 * row + 4] - (plus[row] - minus[row]) / 2e-5) < 1e-8

for periodic in ([True, False, False], [True, False, True]):
    with D4(
        [6, 8],
        [0.4, 0.8, 1.2, 3.1, 2.2, 1.7],
        "pbe",
        charge=0.5,
        lattice=[9.0, 0.2, 0.0, 0.4, 10.0, 0.1, 0.2, 0.3, 11.0],
        periodic=periodic,
    ) as model:
        charges = []
        for cutoff in (8.0, 12.0, 30.0, 60.0):
            model.set_charge_cutoff(cutoff)
            charges.append(model.properties()[1][0])
        differences = [abs(right - left) for left, right in zip(charges, charges[1:])]
        assert max(differences) < 1e-12, differences
        model.set_charge_cutoff(60.0)
        original = model.dispersion(gradient=True)
        model.update([8.9, 1.1, 1.4, 11.6, 2.5, 1.9])
        translated = model.dispersion(gradient=True)
        assert abs(original[0] - translated[0]) < 1e-13
        assert max(abs(left - right) for left, right in zip(original[1], translated[1])) < 1e-12

for kind in (D4.D4, D4.D4S):
    with D4(
        water_numbers, water_positions, "pbe", model=kind, charge=0.5, charge_model="eeqbc"
    ) as model:
        energy, gradient, virial = model.dispersion(gradient=True)
        assert abs(sum(model.properties()[1]) - 0.5) < 1e-13
        pair2, pair3 = model.pairwise()
        assert abs(sum(pair2) + sum(pair3) - energy) < 1e-14
        shifted = water_positions.copy()
        shifted[4] += 1e-5
        model.update(shifted)
        plus = model.dispersion()[0]
        shifted[4] -= 2e-5
        model.update(shifted)
        minus = model.dispersion()[0]
        assert abs(gradient[4] - (plus - minus) / 2e-5) < 1e-9
    with D4(water_numbers, water_positions, "pbe", model=kind, charge=0.5) as model:
        assert abs(model.dispersion()[0] - energy) > 1e-8
try:
    D4(water_numbers, water_positions, "pbe", charge_model="unknown")
except ValueError:
    pass
else:
    raise AssertionError("unknown charge model accepted")

for api in (D3, D4):
    for atm in (False, True):
        with api(water_numbers, water_positions, "pbe", atm=atm) as model:
            for scale in (1.0, 1.2):
                model.update([value * scale for value in water_positions])
                pair2, pair3 = model.pairwise()
                assert len(pair2) == len(pair3) == 9
                assert abs(sum(pair2) + sum(pair3) - model.dispersion()[0]) < 1e-14
                for row in range(3):
                    for column in range(3):
                        assert pair2[3 * row + column] == pair2[3 * column + row]
                        assert pair3[3 * row + column] == pair3[3 * column + row]
                assert any(value != 0.0 for value in pair3) == atm

for damping, values in enumerate(
    (
        {"s8": 0.722, "rs6": 1.217},
        {"s8": 0.7875, "a1": 0.4289, "a2": 4.4407},
        {"s8": 0.0, "rs6": 2.340218, "bet": 0.129434},
        {"s8": 0.358940, "a1": 0.012092, "a2": 5.938951},
        {"s6": 0.91826, "s8": 0.0, "a1": 0.200, "a2": 4.750, "bet": 6.0},
        {"a1": 0.24},
        {"s8": 1.0, "a1": 200770.0},
    )
):
    with D3(water_numbers, water_positions, damping=damping, atm=True, parameters=values) as model:
        energy, gradient, virial = model.dispersion(gradient=True)
        assert all(math.isfinite(value) for value in [energy, *gradient, *virial])
        pair2, pair3 = model.pairwise()
        assert abs(sum(pair2) + sum(pair3) - energy) < 1e-14
        if damping != D3.Z:
            with D3(water_numbers, water_positions, "pbe", damping=damping, atm=True) as named:
                assert model.dispersion(gradient=True) == named.dispersion(gradient=True)
            loaded = D3.damping_parameters("pbe", damping, atm=True)
            with D3(
                water_numbers, water_positions, damping=damping, atm=True, parameters=loaded
            ) as loaded_model:
                assert loaded_model.dispersion(gradient=True) == model.dispersion(gradient=True)

assert D3.damping_parameters("pbe")["s9"] == 0.0
for method in ("pbe", "r2scan-3c", "dftb(3ob)", "lc-dftb"):
    for atm in (False, True):
        loaded = D4.damping_parameters(method, s9=None if atm else 0.0)
        with D4(water_numbers, water_positions, method, atm=atm) as named:
            expected = named.dispersion(gradient=True)
        with D4(water_numbers, water_positions, parameters=loaded, atm=atm) as explicit:
            assert explicit.dispersion(gradient=True) == expected
assert D4.damping_parameters("dftb(3ob)", s9=0.0)["s8"] == 0.4727337
assert D4.damping_parameters("pbe", s9=0.5)["s9"] == 0.5
for factory in (D3, D4):
    try:
        factory.damping_parameters("unknown")
    except RuntimeError:
        pass
    else:
        raise AssertionError("unknown named parameters accepted")
try:
    D4.damping_parameters("pbe", s9=float("nan"))
except RuntimeError:
    pass
else:
    raise AssertionError("nonfinite ATM scaling accepted")

for kind in (D4.D4, D4.D4S):
    for atm in (False, True):
        with D4(water_numbers, water_positions, "pbe", model=kind, atm=atm) as named:
            expected = named.dispersion(gradient=True)
        with D4(
            water_numbers,
            water_positions,
            model=kind,
            atm=atm,
            parameters={"s8": 0.95948085, "a1": 0.38574991, "a2": 4.80688534},
        ) as explicit:
            assert explicit.dispersion(gradient=True) == expected

for parameters in ({}, {"s8": float("nan"), "a1": 0.4, "a2": 4.0}):
    try:
        D4(numbers, positions, parameters=parameters)
    except (ValueError, RuntimeError):
        pass
    else:
        raise AssertionError("invalid explicit D4 parameters accepted")

with D4(
    numbers,
    positions,
    "pbe",
    model=D4.D4S,
    atm=False,
    lattice=[20.0, 0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 20.0],
    periodic=[True, True, True],
) as model:
    energy, gradient, virial = model.dispersion(gradient=True)
    assert math.isfinite(energy) and energy < 0.0
    assert all(math.isfinite(value) for value in gradient + virial)
