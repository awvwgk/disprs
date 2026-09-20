"""ASE dispersion calculator; combine with another calculator using SumCalculator."""

from ase.calculators.calculator import Calculator, all_changes
from ase.units import Bohr, Hartree

from . import D3
from ._framework import evaluate


class Disprs(Calculator):
    implemented_properties = ["energy", "free_energy", "forces", "stress"]
    default_parameters = dict(method="pbe", level="d4", atm=True, damping=D3.RATIONAL, charge=None)

    def set(self, **kwargs):
        changed = super().set(**kwargs)
        if changed:
            self.reset()
        return changed

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        atoms = self.atoms
        charge = self.parameters.charge
        if charge is None:
            charge = float(atoms.get_initial_charges().sum())
        energy, gradient, virial = evaluate(
            atoms.numbers.tolist(),
            atoms.positions / Bohr,
            self.parameters.method,
            level=self.parameters.level,
            charge=charge,
            atm=self.parameters.atm,
            damping=self.parameters.damping,
            lattice=atoms.cell.array / Bohr,
            periodic=atoms.pbc,
        )
        self.results = dict(
            energy=energy * Hartree, free_energy=energy * Hartree, forces=-gradient * Hartree / Bohr
        )
        if atoms.pbc.all():
            self.results["stress"] = (virial * Hartree / atoms.get_volume()).flat[
                [0, 4, 8, 5, 2, 1]
            ]
