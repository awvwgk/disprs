import importlib.util
import sys
import unittest


@unittest.skipUnless(importlib.util.find_spec("ase"), "ASE is not installed")
class TestASE(unittest.TestCase):
    def test_energy_forces_and_updates(self):
        import numpy as np
        from ase import Atoms
        from ase.units import Bohr, Hartree
        from disprs import D3, D4
        from disprs.ase_ext import Disprs

        atoms = Atoms("CO", positions=[[0, 0, 0], [2.0, 0, 0]])
        atoms.calc = Disprs(method="pbe")
        for level, model_type in (("d3", D3), ("d3s", D3), ("d4", D4), ("d4s", D4)):
            atoms.calc.set(level=level)
            options = dict(model=D4.D4S) if level == "d4s" else {}
            if level == "d3s":
                options = dict(d3s=True)
            with model_type(
                [6, 8], (atoms.positions / Bohr).ravel(), "pbe", atm=True, **options
            ) as model:
                energy, gradient, _ = model.dispersion(gradient=True)
            self.assertAlmostEqual(atoms.get_potential_energy(), energy * Hartree)
            np.testing.assert_allclose(
                atoms.get_forces(), -np.reshape(gradient, (2, 3)) * Hartree / Bohr
            )
        before = atoms.get_potential_energy()
        atoms.positions[1, 0] += 0.2
        self.assertNotEqual(before, atoms.get_potential_energy())

    def test_periodic_stress_and_charge(self):
        import numpy as np
        from ase import Atoms
        from ase.units import Bohr, Hartree
        from disprs import D4
        from disprs.ase_ext import Disprs

        atoms = Atoms("CO", positions=[[0, 0, 0], [2, 0.3, 0]], cell=[8, 9, 10], pbc=True)
        atoms.set_initial_charges([1, 0])
        atoms.calc = Disprs(atm=False)
        with D4(
            [6, 8],
            (atoms.positions / Bohr).ravel(),
            "pbe",
            charge=1,
            atm=False,
            lattice=(atoms.cell.array / Bohr).ravel(),
            periodic=[True] * 3,
        ) as model:
            energy, _, virial = model.dispersion(gradient=True)
        self.assertAlmostEqual(atoms.get_potential_energy(), energy * Hartree)
        np.testing.assert_allclose(
            atoms.get_stress(),
            np.asarray(virial)[[0, 4, 8, 5, 2, 1]] * Hartree / atoms.get_volume(),
        )
        stress = atoms.get_stress()[0]
        volume = atoms.get_volume()
        delta = 1e-5
        energies = []
        for strain in (-delta, delta):
            displaced = atoms.copy()
            displaced.calc = Disprs(atm=False)
            cell = atoms.cell.array.copy()
            cell[:, 0] *= 1 + strain
            displaced.set_cell(cell, scale_atoms=True)
            energies.append(displaced.get_potential_energy())
        self.assertAlmostEqual(stress, (energies[1] - energies[0]) / (2 * delta * volume), places=8)

    def test_partial_periodicity(self):
        import numpy as np
        from ase import Atoms
        from disprs.ase_ext import Disprs

        for level in ("d3", "d3s", "d4", "d4s"):
            for periodic in ([True, False, False], [True, False, True]):
                atoms = Atoms(
                    "CO", positions=[[0.2, 0.3, 0.1], [2, 0.4, 0.2]], cell=[10] * 3, pbc=periodic
                )
                atoms.calc = Disprs(level=level)
                energy = atoms.get_potential_energy()
                forces = atoms.get_forces()
                cell = atoms.cell.array.copy()
                cell[~np.asarray(periodic)] = 0
                atoms.set_cell(cell)
                self.assertAlmostEqual(atoms.get_potential_energy(), energy, places=12)
                np.testing.assert_allclose(atoms.get_forces(), forces, atol=1e-12)
                atoms.positions[0] += atoms.cell[0]
                self.assertAlmostEqual(atoms.get_potential_energy(), energy, places=12)


@unittest.skipUnless(importlib.util.find_spec("pyscf"), "PySCF is not installed")
class TestPySCF(unittest.TestCase):
    def test_scf_energy_gradient_and_reset(self):
        import numpy as np
        from disprs.pyscf_ext import Dispersion, energy
        from pyscf import gto, lib, scf

        lib.num_threads(1)
        mol = gto.M(atom="H 0 0 0; H 0 0 1.4", unit="Bohr", basis="sto-3g", verbose=0)
        plain = scf.RHF(mol).run()
        for level in ("d3", "d3s", "d4", "d4s"):
            corrected = energy(plain, level=level).run()
            correction, gradient = Dispersion(mol, level=level).kernel()
            self.assertAlmostEqual(corrected.e_tot - plain.e_tot, correction, places=10)
            np.testing.assert_allclose(
                corrected.nuc_grad_method().kernel() - plain.nuc_grad_method().kernel(),
                gradient,
                atol=1e-9,
            )
            np.testing.assert_allclose(
                corrected.Gradients().grad_nuc(atmlst=[1]) - plain.Gradients().grad_nuc(atmlst=[1]),
                gradient[[1]],
                atol=1e-12,
            )
            repeated = energy(corrected, level=level)
            self.assertAlmostEqual(repeated.energy_nuc(), corrected.energy_nuc())
            moved = mol.copy().set_geom_("H 0 0 0; H 0 0 1.8", unit="Bohr")
            corrected.reset(moved)
            self.assertAlmostEqual(
                corrected.energy_nuc() - moved.energy_nuc(),
                Dispersion(moved, level=level).kernel()[0],
            )
        self.assertFalse(hasattr(plain, "with_disprs"))
        self.assertNotIn("disprs", plain.scf_summary)


@unittest.skipUnless(importlib.util.find_spec("qcelemental"), "QCElemental is not installed")
class TestQCSchema(unittest.TestCase):
    @unittest.skipIf(sys.version_info >= (3, 14), "QCSchema v1 requires Python <3.14")
    def test_v1(self):
        from disprs.qcschema_ext import run_qcschema
        from qcelemental.models import v1

        data = v1.AtomicInput(
            molecule=dict(symbols=["H", "H"], geometry=[0, 0, 0, 0, 0, 1.4]),
            driver="energy",
            model=dict(method="hf"),
        )
        self.assertIsInstance(run_qcschema(data), v1.AtomicResult)
        invalid = data.dict()
        invalid["keywords"] = {"level": "invalid"}
        self.assertIsInstance(run_qcschema(invalid), v1.FailedOperation)

    def test_ghosts_charge_and_levels(self):
        import numpy as np
        from disprs import D3, D4
        from disprs.qcschema_ext import run_qcschema
        from qcelemental.models import v2

        molecule = v2.Molecule(
            symbols=["H", "H", "He"],
            geometry=[0, 0, 0, 0, 0, 1.4, 0, 4, 0],
            real=[True, True, False],
            molecular_charge=1,
            molecular_multiplicity=2,
        )
        for level in ("d3", "d3s", "d4", "d4s"):
            result = run_qcschema(
                v2.AtomicInput(
                    molecule=molecule,
                    specification=dict(
                        driver="gradient", model=dict(method="pbe"), keywords=dict(level=level)
                    ),
                )
            )
            self.assertTrue(result.success)
            options = (
                dict(d3s=level == "d3s")
                if level in ("d3", "d3s")
                else dict(charge=1, model=D4.D4S if level == "d4s" else D4.D4)
            )
            with (D3 if level in ("d3", "d3s") else D4)(
                [1, 1], [0, 0, 0, 0, 0, 1.4], "pbe", atm=True, **options
            ) as model:
                _, gradient, _ = model.dispersion(gradient=True)
            np.testing.assert_allclose(result.return_result[:2], np.reshape(gradient, (2, 3)))
            np.testing.assert_array_equal(result.return_result[2], 0)

    def test_results_and_failures(self):
        import numpy as np
        from disprs import D4
        from disprs.qcschema_ext import run_qcschema
        from qcelemental.models.v2 import AtomicInput, AtomicResult, FailedOperation

        data = dict(
            molecule=dict(symbols=["H", "H"], geometry=[0, 0, 0, 0, 0, 1.4]),
            specification=dict(driver="energy", model=dict(method="pbe"), extras={"keep": True}),
        )
        with D4([1, 1], [0, 0, 0, 0, 0, 1.4], "pbe") as model:
            expected, gradient, _ = model.dispersion(gradient=True)
        result = run_qcschema(data)
        self.assertIsInstance(result, AtomicResult)
        self.assertTrue(result.success)
        self.assertAlmostEqual(result.return_result, expected)
        self.assertTrue(result.input_data.specification.extras["keep"])
        data["specification"]["driver"] = "gradient"
        result = run_qcschema(AtomicInput(**data))
        np.testing.assert_allclose(result.return_result, np.reshape(gradient, (2, 3)))
        data["specification"]["driver"] = "hessian"
        self.assertIsInstance(run_qcschema(data), FailedOperation)
        data["specification"]["driver"] = "energy"
        data["specification"]["keywords"] = {"ignored_typo": True}
        self.assertFalse(run_qcschema(data).success)


if __name__ == "__main__":
    unittest.main()
