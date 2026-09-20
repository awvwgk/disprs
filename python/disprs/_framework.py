import numpy as np

from . import D3, D4


def evaluate(
    numbers,
    positions,
    method,
    *,
    level="d4",
    charge=0.0,
    atm=True,
    damping=D3.RATIONAL,
    lattice=None,
    periodic=None,
):
    """Evaluate a correction in atomic units for the optional framework adapters."""
    level = level.lower()
    if level not in ("d3", "d3s", "d4", "d4s"):
        raise ValueError("level must be 'd3', 'd3s', 'd4', or 'd4s'")
    periodic = np.asarray(periodic if periodic is not None else [False] * 3, dtype=bool)
    positions = np.asarray(positions, dtype=float).reshape(-1).tolist()
    lattice = None if lattice is None else np.asarray(lattice, dtype=float).reshape(-1).tolist()
    options = dict(atm=atm, lattice=lattice, periodic=periodic.tolist())
    if level in ("d3", "d3s"):
        model = D3(
            numbers,
            positions,
            method,
            damping=damping,
            model=D3.D3S if level == "d3s" else D3.D3,
            **options,
        )
    else:
        model = D4(
            numbers,
            positions,
            method,
            charge=charge,
            model=D4.D4S if level == "d4s" else D4.D4,
            **options,
        )
    with model:
        energy, gradient, virial = model.dispersion(gradient=True)
    return energy, np.asarray(gradient).reshape(-1, 3), np.asarray(virial).reshape(3, 3)
