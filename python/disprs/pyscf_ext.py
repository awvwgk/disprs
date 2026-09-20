"""Molecular PySCF SCF energy and nuclear-gradient corrections."""

import numpy as np
from pyscf import lib
from pyscf.scf.hf import SCF

from . import D3
from ._framework import evaluate


class Dispersion(lib.StreamObject):
    def __init__(self, mol, method="hf", level="d4", atm=True, damping=D3.RATIONAL):
        self.mol = mol
        self.method = method
        self.level = level
        self.atm = atm
        self.damping = damping

    def kernel(self, mol=None):
        mol = self.mol if mol is None else mol
        if hasattr(mol, "pbc_intor"):
            raise ValueError("The PySCF adapter supports molecular calculations only")
        numbers = mol.atom_charges()
        if np.any(numbers <= 0):
            raise ValueError("Ghost atoms are not supported by the PySCF adapter")
        energy, gradient, _ = evaluate(
            numbers.tolist(),
            mol.atom_coords(),
            self.method,
            level=self.level,
            charge=mol.charge,
            atm=self.atm,
            damping=self.damping,
        )
        return energy, gradient


class _DisprsSCF:
    def energy_nuc(self):
        energy = super().energy_nuc()
        correction, _ = self.with_disprs.kernel(self.mol)
        self.scf_summary["disprs"] = correction
        return energy + correction

    def reset(self, mol=None):
        result = super().reset(mol)
        self.with_disprs.mol = self.mol
        return result

    def nuc_grad_method(self):
        gradient = super().nuc_grad_method()
        return lib.set_class(gradient, (_DisprsGrad, gradient.__class__))

    Gradients = nuc_grad_method


class _DisprsGrad:
    def grad_nuc(self, mol=None, atmlst=None):
        mol = self.mol if mol is None else mol
        gradient = super().grad_nuc(mol, atmlst)
        _, correction = self.base.with_disprs.kernel(mol)
        return gradient + (correction if atmlst is None else correction[atmlst])


def energy(mf, *, method=None, level="d4", atm=True, damping=D3.RATIONAL):
    """Return a dispersion-corrected SCF copy without modifying the input method."""
    if not isinstance(mf, SCF):
        raise TypeError("Expected a PySCF SCF method")
    if (
        getattr(mf, "disp", None)
        or getattr(mf, "with_dftd3", None)
        or getattr(mf, "with_dftd4", None)
    ):
        raise ValueError("Remove the existing dispersion correction before adding disprs")
    method = method or getattr(mf, "xc", "hf")
    xc = getattr(mf, "xc", "").lower()
    if any(suffix in value for suffix in ("-d3", "-d4") for value in (method.lower(), xc)):
        raise ValueError("Use an XC method without a dispersion suffix")
    corrected = mf.copy()
    corrected.scf_summary = dict(mf.scf_summary)
    corrected.with_disprs = Dispersion(mf.mol, method, level, atm, damping)
    corrected._keys = set(mf._keys) | {"with_disprs"}
    if not isinstance(corrected, _DisprsSCF):
        lib.set_class(corrected, (_DisprsSCF, mf.__class__))
    return corrected
