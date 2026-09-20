"""QCSchema v1/v2 energy and gradient driver using QCElemental models."""

import sys

import numpy as np
from qcelemental.models import v1, v2

from . import version
from ._framework import evaluate


def run_qcschema(input_data):
    """Return an AtomicResult or FailedOperation; geometry is in Bohr."""
    original = input_data
    schema = (
        v2
        if (
            isinstance(input_data, v2.AtomicInput)
            or isinstance(input_data, dict)
            and (input_data.get("schema_version") == 2 or "specification" in input_data)
        )
        else v1
    )
    try:
        data = (
            input_data
            if isinstance(input_data, schema.AtomicInput)
            else schema.AtomicInput(**input_data)
        )
        specification = data.specification if schema is v2 else data
        if specification.driver not in ("energy", "gradient"):
            raise ValueError("Only energy and gradient drivers are supported")
        keywords = dict(specification.keywords)
        unknown = keywords.keys() - {"level", "atm", "damping"}
        if unknown:
            raise ValueError(f"Unsupported keywords: {sorted(unknown)}")
        molecule = data.molecule
        real = np.asarray(molecule.real, dtype=bool)
        if not real.any():
            raise ValueError("At least one real atom is required")
        energy, gradient, _ = evaluate(
            molecule.atomic_numbers[real].tolist(),
            molecule.geometry[real],
            specification.model.method,
            charge=molecule.molecular_charge,
            **keywords,
        )
        full_gradient = np.zeros_like(molecule.geometry)
        full_gradient[real] = gradient
        result = dict(input_data=data, molecule=molecule) if schema is v2 else data.dict()
        result.update(
            success=True,
            return_result=energy if specification.driver == "energy" else full_gradient,
            properties={"return_energy": energy},
            provenance={
                "creator": "disprs",
                "version": version(),
                "routine": "disprs.qcschema_ext.run_qcschema",
            },
        )
        return schema.AtomicResult(**result)
    except (ValueError, TypeError, RuntimeError) as error:
        failure_schema = v2 if sys.version_info >= (3, 14) else schema
        return failure_schema.FailedOperation(
            input_data=original,
            error={"error_type": "input_error", "error_message": str(error)},
        )
