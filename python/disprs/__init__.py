import ctypes
import ctypes.util
import math
import os
from pathlib import Path

_bundled = next(
    (
        path
        for name in ("libdisprs.so", "libdisprs.dylib", "disprs.dll")
        if (path := Path(__file__).with_name(name)).is_file()
    ),
    None,
)
_path = os.environ.get("DISPRS_LIBRARY") or _bundled or ctypes.util.find_library("disprs")
if not _path:
    raise ImportError("libdisprs was not found; set DISPRS_LIBRARY to its path")
_lib = ctypes.CDLL(_path)
_lib.disprs_get_version.restype = ctypes.c_char_p

_handle = ctypes.c_void_p
_int_p = ctypes.POINTER(ctypes.c_int)
_double_p = ctypes.POINTER(ctypes.c_double)
_handle_p = ctypes.POINTER(_handle)


def _prototype(name, arguments, result=None):
    function = getattr(_lib, name)
    function.argtypes = arguments
    function.restype = result


for _prefix in ("d3", "d4"):
    _prototype(
        f"disprs_{_prefix}_get_pairwise_dispersion",
        [_handle, _handle, _handle, _handle, _double_p, _double_p],
    )
    _prototype(f"disprs_{_prefix}_new_error", [], _handle)
    _prototype(f"disprs_{_prefix}_check_error", [_handle], ctypes.c_int)
    _prototype(f"disprs_{_prefix}_get_error", [_handle, ctypes.c_void_p, _int_p])
    _prototype(f"disprs_{_prefix}_delete_error", [_handle_p])
    _prototype(
        f"disprs_{_prefix}_set_model_ghost_index",
        [_handle, _handle, _int_p, ctypes.c_int],
    )
    _prototype(
        f"disprs_{_prefix}_set_model_work_partition",
        [_handle, _handle, ctypes.c_int, ctypes.c_int],
    )
    _prototype(
        f"disprs_{_prefix}_set_model_ewald",
        [_handle, _handle, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int],
    )

_prototype(
    "disprs_d3_new_structure",
    [_handle, ctypes.c_int, _int_p, _double_p, ctypes.c_void_p, ctypes.c_void_p],
    _handle,
)
_prototype("disprs_d3_update_structure", [_handle, _handle, _double_p, ctypes.c_void_p])
_prototype("disprs_d3_delete_structure", [_handle_p])
_prototype("disprs_d3_new_model", [_handle, _handle], _handle)
_prototype("disprs_d3_new_smooth_model", [_handle, _handle], _handle)
_prototype("disprs_d3_new_model_kind", [_handle, _handle, ctypes.c_int], _handle)
_prototype("disprs_d3_delete_model", [_handle_p])
_prototype("disprs_d3_get_properties", [_handle, _handle, _handle, _double_p, _double_p])
_prototype("disprs_d3_get_property_response", [_handle] * 3 + [_double_p] * 6)
_prototype(
    "disprs_d3_set_model_realspace_cutoff_smooth",
    [_handle, _handle] + [ctypes.c_double] * 5,
)
_prototype(
    "disprs_d3_load_param",
    [_handle, ctypes.c_int, ctypes.c_char_p, ctypes.c_bool],
    _handle,
)
_prototype("disprs_d3_delete_param", [_handle_p])
_prototype(
    "disprs_d3_get_named_parameters",
    [_handle, ctypes.c_int, ctypes.c_char_p, _double_p],
)
_d3_parameter_keys = (
    ("zero", "s6 s8 s9 rs6 rs8 alp"),
    ("rational", "s6 s8 s9 a1 a2 alp"),
    ("mzero", "s6 s8 s9 rs6 rs8 alp bet"),
    ("mrational", "s6 s8 s9 a1 a2 alp"),
    ("optimizedpower", "s6 s8 s9 a1 a2 alp bet"),
    ("cso", "s6 s9 a1 a2 a3 a4 alp"),
    ("z", "s6 s8 s9 a1 alp"),
)
for _name, _keys in _d3_parameter_keys:
    _prototype(
        f"disprs_d3_new_{_name}_damping",
        [_handle] + [ctypes.c_double] * len(_keys.split()),
        _handle,
    )
_prototype(
    "disprs_d3_get_dispersion",
    [_handle, _handle, _handle, _handle, _double_p, ctypes.c_void_p, ctypes.c_void_p],
)
_prototype(
    "disprs_d3_get_dispersion_hessian",
    [_handle, _handle, _handle, _handle, _double_p, _double_p],
)
_prototype(
    "disprs_d3_load_gcp",
    [_handle, _handle, ctypes.c_char_p, ctypes.c_char_p],
    _handle,
)
_prototype("disprs_d3_delete_gcp", [_handle_p])
_prototype(
    "disprs_d3_set_gcp_controls", [_handle, _handle, _double_p, ctypes.c_void_p, ctypes.c_void_p]
)
for _operation in ("get", "set"):
    _prototype(
        f"disprs_d3_{_operation}_gcp_parameters",
        [_handle, _handle, ctypes.c_int] + [ctypes.c_void_p] * 8,
    )
_prototype(
    "disprs_d3_set_gcp_realspace_cutoff",
    [_handle, _handle, ctypes.c_double, ctypes.c_double],
)
_prototype(
    "disprs_d3_get_counterpoise",
    [_handle, _handle, _handle, _double_p, ctypes.c_void_p, ctypes.c_void_p],
)
_prototype(
    "disprs_d3_get_counterpoise_hessian",
    [_handle, _handle, _handle, _double_p, _double_p],
)

_prototype(
    "disprs_d4_new_structure",
    [_handle, ctypes.c_int, _int_p, _double_p, _double_p, ctypes.c_void_p, ctypes.c_void_p],
    _handle,
)
_prototype("disprs_d4_update_structure", [_handle, _handle, _double_p, ctypes.c_void_p])
_prototype("disprs_d4_delete_structure", [_handle_p])
_prototype("disprs_d4_new_model", [_handle, _handle, ctypes.c_int], _handle)
_prototype(
    "disprs_d4_new_custom_model", [_handle, _handle, ctypes.c_int] + [ctypes.c_double] * 3, _handle
)
_prototype("disprs_d4_set_realspace_cutoff", [_handle, _handle] + [ctypes.c_double] * 5)
_prototype("disprs_d4_set_charge_model", [_handle, _handle, ctypes.c_int])
_prototype("disprs_d4_set_fixed_charges", [_handle, _handle, _double_p, ctypes.c_int])
_prototype("disprs_d4_set_charge_cutoff", [_handle, _handle, ctypes.c_double])
_prototype("disprs_d4_delete_model", [_handle_p])
_prototype("disprs_d4_load_param", [_handle, ctypes.c_char_p, ctypes.c_bool], _handle)
_prototype("disprs_d4_new_rational_damping", [_handle] + [ctypes.c_double] * 6, _handle)
_prototype("disprs_d4_delete_param", [_handle_p])
_prototype(
    "disprs_d4_get_named_parameters_s9",
    [_handle, ctypes.c_char_p, _double_p, _double_p],
)
_prototype(
    "disprs_d4_get_dispersion",
    [_handle, _handle, _handle, _handle, _double_p, ctypes.c_void_p, ctypes.c_void_p],
)
_prototype(
    "disprs_d4_get_properties",
    [_handle, _handle, _handle, _double_p, _double_p, _double_p, _double_p],
)
_prototype("disprs_d4_get_property_response", [_handle] * 3 + [_double_p] * 12)
_prototype(
    "disprs_d4_get_dispersion_hessian",
    [_handle, _handle, _handle, _handle, _double_p],
)
_prototype(
    "disprs_get_charges",
    [_handle, _handle, ctypes.c_int, _double_p, _double_p, _double_p],
)


def version():
    return _lib.disprs_get_version().decode()


def _array(kind, values, size=None):
    values = tuple(values)
    if size is not None and len(values) != size:
        raise ValueError(f"expected {size} values, got {len(values)}")
    return (kind * len(values))(*values)


def _check(prefix, error):
    if not getattr(_lib, f"disprs_{prefix}_check_error")(error):
        return
    size = ctypes.c_int(512)
    message = ctypes.create_string_buffer(size.value)
    getattr(_lib, f"disprs_{prefix}_get_error")(error, message, ctypes.byref(size))
    raise RuntimeError(message.value.decode())


def _delete(name, handle):
    if handle:
        value = ctypes.c_void_p(handle)
        getattr(_lib, name)(ctypes.byref(value))


def get_charges(
    numbers,
    positions,
    charge=0.0,
    charge_model="eeq",
    lattice=None,
    periodic=None,
    *,
    cartesian=False,
    strain=False,
):
    """Return charges and optional flat Fortran-order (3,N,N)/(3,3,N) responses."""
    if charge_model not in ("eeq", "eeqbc"):
        raise ValueError("charge_model must be 'eeq' or 'eeqbc'")
    numbers = _array(ctypes.c_int, numbers)
    natoms = len(numbers)
    positions = _array(ctypes.c_double, positions, 3 * natoms)
    total_charge = ctypes.c_double(charge)
    lattice = _array(ctypes.c_double, lattice, 9) if lattice is not None else None
    periodic = _array(ctypes.c_bool, periodic, 3) if periodic is not None else None
    charges = (ctypes.c_double * natoms)()
    dqdr = (ctypes.c_double * (3 * natoms * natoms))() if cartesian else None
    dqdstrain = (ctypes.c_double * (9 * natoms))() if strain else None
    error = _lib.disprs_d4_new_error()
    structure = None
    try:
        structure = _lib.disprs_d4_new_structure(
            error, natoms, numbers, positions, ctypes.byref(total_charge), lattice, periodic
        )
        _check("d4", error)
        _lib.disprs_get_charges(
            error, structure, int(charge_model == "eeqbc"), charges, dqdr, dqdstrain
        )
        _check("d4", error)
        return (
            list(charges),
            list(dqdr) if dqdr is not None else None,
            list(dqdstrain) if dqdstrain is not None else None,
        )
    finally:
        _delete("disprs_d4_delete_structure", structure)
        _delete("disprs_d4_delete_error", error)


class D3:
    """Native D3/D3S; select model=D3.D3S (elements 1-94) for smooth C6."""

    D3 = 0
    D3S = 1
    ZERO = 0
    RATIONAL = 1
    MODIFIED_ZERO = 2
    MODIFIED_RATIONAL = 3
    OPTIMIZED_POWER = 4
    CSO = 5
    Z = 6

    @staticmethod
    def damping_parameters(method, damping=RATIONAL, *, atm=False):
        """Return named damping values accepted by the parameters constructor option."""
        if damping not in range(len(_d3_parameter_keys)):
            raise ValueError("invalid D3 damping function")
        method = method.encode()
        values = (ctypes.c_double * 9)()
        error = _lib.disprs_d3_new_error()
        try:
            _lib.disprs_d3_get_named_parameters(error, damping, method, values)
            _check("d3", error)
            result = dict(zip("s6 s8 s9 rs6 rs8 a1 a2 alp bet".split(), values))
            result["s9"] = float(bool(atm))
            if damping == D3.CSO:
                result.update(a3=result["rs6"], a4=result["rs8"])
            return {key: result[key] for key in _d3_parameter_keys[damping][1].split()}
        finally:
            _delete("disprs_d3_delete_error", error)

    def __init__(
        self,
        numbers,
        positions,
        method=None,
        damping=RATIONAL,
        atm=False,
        ghosts=(),
        lattice=None,
        periodic=None,
        *,
        parameters=None,
        model=None,
        d3s=None,
    ):
        if model is None:
            model = self.D3S if d3s else self.D3
        elif d3s is not None and model != int(bool(d3s)):
            raise ValueError("model and d3s select different D3 models")
        if model not in (self.D3, self.D3S):
            raise ValueError("model must be D3.D3 or D3.D3S")
        if (method is None) == (parameters is None):
            raise ValueError("specify exactly one of method or parameters")
        if damping not in range(len(_d3_parameter_keys)):
            raise ValueError("invalid D3 damping function")
        if parameters is not None:
            name, keys = _d3_parameter_keys[damping]
            keys = keys.split()
            defaults = {"s6": 1.0, "s9": 1.0, "alp": 14.0, "rs8": 1.0}
            if damping == self.CSO:
                defaults.update(a2=2.5, a3=0.0, a4=6.25)
            values = {key: value for key, value in defaults.items() if key in keys}
            values.update(parameters)
            if set(values) != set(keys):
                raise ValueError(f"parameters for {name} require {', '.join(keys)}")
            values = [float(values[key]) for key in keys]
            if not all(math.isfinite(value) for value in values):
                raise ValueError("D3 damping parameters must be finite")
            if not atm:
                values[keys.index("s9")] = 0.0
        self.natoms = len(numbers)
        self.error = _lib.disprs_d3_new_error()
        self.structure = self.model = self.param = None
        numbers = _array(ctypes.c_int, numbers)
        positions = _array(ctypes.c_double, positions, 3 * self.natoms)
        lattice = _array(ctypes.c_double, lattice, 9) if lattice is not None else None
        periodic = _array(ctypes.c_bool, periodic, 3) if periodic is not None else None
        self.structure = _lib.disprs_d3_new_structure(
            self.error, self.natoms, numbers, positions, lattice, periodic
        )
        _check("d3", self.error)
        try:
            self.model = _lib.disprs_d3_new_model_kind(self.error, self.structure, model)
            _check("d3", self.error)
            self.set_ghosts(ghosts)
            if parameters is None:
                method = ctypes.create_string_buffer(method.encode())
                self.param = _lib.disprs_d3_load_param(self.error, damping, method, bool(atm))
            else:
                self.param = getattr(_lib, f"disprs_d3_new_{name}_damping")(self.error, *values)
            _check("d3", self.error)
            if not self.param:
                raise ValueError("invalid D3 damping function")
        except Exception:
            self.close()
            raise

    def update(self, positions, lattice=None):
        positions = _array(ctypes.c_double, positions, 3 * self.natoms)
        lattice = _array(ctypes.c_double, lattice, 9) if lattice is not None else None
        _lib.disprs_d3_update_structure(self.error, self.structure, positions, lattice)
        _check("d3", self.error)

    def set_ghosts(self, indices):
        """Cumulatively exclude zero-based atoms from dispersion, not CN."""
        indices = tuple(indices)
        if any(index < 0 or index >= self.natoms for index in indices):
            raise ValueError("D3 ghost index is out of range")
        indices = _array(ctypes.c_int, indices)
        _lib.disprs_d3_set_model_ghost_index(self.error, self.model, indices, len(indices))
        _check("d3", self.error)

    def set_work_partition(self, part=0, parts=1):
        """Select a zero-based work partition; sum dispersion outputs across parts."""
        _lib.disprs_d3_set_model_work_partition(self.error, self.model, part, parts)
        _check("d3", self.error)

    def set_realspace_cutoff(
        self,
        dispersion2=60.0,
        dispersion3=40.0,
        coordination=40.0,
        width2=0.0,
        width3=0.0,
    ):
        _lib.disprs_d3_set_model_realspace_cutoff_smooth(
            self.error,
            self.model,
            dispersion2,
            dispersion3,
            coordination,
            width2,
            width3,
        )
        _check("d3", self.error)

    def dispersion(self, gradient=False):
        energy = ctypes.c_double()
        grad = (ctypes.c_double * (3 * self.natoms))() if gradient else None
        virial = (ctypes.c_double * 9)() if gradient else None
        _lib.disprs_d3_get_dispersion(
            self.error, self.structure, self.model, self.param, ctypes.byref(energy), grad, virial
        )
        _check("d3", self.error)
        return energy.value, list(grad) if grad else None, list(virial) if virial else None

    def pairwise(self):
        """Return flat two-body and ATM matrices; their sums give the energy."""
        pair2 = (ctypes.c_double * (self.natoms * self.natoms))()
        pair3 = (ctypes.c_double * (self.natoms * self.natoms))()
        _lib.disprs_d3_get_pairwise_dispersion(
            self.error, self.structure, self.model, self.param, pair2, pair3
        )
        _check("d3", self.error)
        return list(pair2), list(pair3)

    def properties(self):
        """Return full-system CN and flat undamped C6, including ghost atoms."""
        coordination = (ctypes.c_double * self.natoms)()
        c6 = (ctypes.c_double * (self.natoms * self.natoms))()
        _lib.disprs_d3_get_properties(self.error, self.structure, self.model, coordination, c6)
        _check("d3", self.error)
        return list(coordination), list(c6)

    def property_response(self):
        """Return {name: (values, Cartesian, strain)} as flat Fortran-order arrays."""
        sizes = (self.natoms, self.natoms**2)
        values = [(ctypes.c_double * size)() for size in sizes]
        responses = [
            ((ctypes.c_double * (3 * self.natoms * size))(), (ctypes.c_double * (9 * size))())
            for size in sizes
        ]
        _lib.disprs_d3_get_property_response(
            self.error,
            self.structure,
            self.model,
            *values,
            *(array for response in responses for array in response),
        )
        _check("d3", self.error)
        return {
            name: (list(value), list(response[0]), list(response[1]))
            for name, value, response in zip(("coordination", "c6"), values, responses)
        }

    def set_ewald(self, *, rank=0, tolerance=1e-4, kcut=0.0, mesh=0):
        """Enable periodic D3 Fourier summation; mesh=-1 selects direct Ewald."""
        _lib.disprs_d3_set_model_ewald(self.error, self.model, rank, tolerance, kcut, mesh)
        _check("d3", self.error)

    def counterpoise_parameters(self, method=None, basis=None, *, eta=None):
        method = ctypes.create_string_buffer(method.encode()) if method else None
        basis = ctypes.create_string_buffer(basis.encode()) if basis else None
        gcp = _lib.disprs_d3_load_gcp(self.error, self.structure, method, basis)
        _check("d3", self.error)
        try:
            if eta is not None:
                _lib.disprs_d3_set_gcp_controls(
                    self.error, gcp, ctypes.byref(ctypes.c_double(eta)), None, None
                )
                _check("d3", self.error)
            return _gcp_parameters(self.error, gcp, self.natoms)
        finally:
            _delete("disprs_d3_delete_gcp", gcp)

    def counterpoise(
        self,
        method,
        basis=None,
        gradient=False,
        cutoff=60.0,
        srb_cutoff=60.0,
        *,
        eta=None,
        base=None,
        srb=None,
        parameters=None,
    ):
        method = ctypes.create_string_buffer(method.encode())
        basis = ctypes.create_string_buffer(basis.encode()) if basis else None
        gcp = _lib.disprs_d3_load_gcp(self.error, self.structure, method, basis)
        _check("d3", self.error)
        try:
            _lib.disprs_d3_set_gcp_controls(
                self.error,
                gcp,
                ctypes.byref(ctypes.c_double(eta)) if eta is not None else None,
                ctypes.byref(ctypes.c_bool(base)) if base is not None else None,
                ctypes.byref(ctypes.c_bool(srb)) if srb is not None else None,
            )
            _check("d3", self.error)
            if parameters is not None:
                _gcp_parameters(self.error, gcp, self.natoms, parameters)
            _lib.disprs_d3_set_gcp_realspace_cutoff(self.error, gcp, cutoff, srb_cutoff)
            energy = ctypes.c_double()
            grad = (ctypes.c_double * (3 * self.natoms))() if gradient else None
            virial = (ctypes.c_double * 9)() if gradient else None
            _lib.disprs_d3_get_counterpoise(
                self.error, self.structure, gcp, ctypes.byref(energy), grad, virial
            )
            _check("d3", self.error)
            return energy.value, list(grad) if grad else None, list(virial) if virial else None
        finally:
            _delete("disprs_d3_delete_gcp", gcp)

    def hessian(self):
        """Return energy and the flat Cartesian Hessian in atomic units."""
        energy = ctypes.c_double()
        hessian = (ctypes.c_double * (3 * self.natoms) ** 2)()
        _lib.disprs_d3_get_dispersion_hessian(
            self.error, self.structure, self.model, self.param, ctypes.byref(energy), hessian
        )
        _check("d3", self.error)
        return energy.value, list(hessian)

    def counterpoise_hessian(
        self,
        method,
        basis=None,
        cutoff=60.0,
        srb_cutoff=60.0,
        *,
        eta=None,
        base=None,
        srb=None,
        parameters=None,
    ):
        method = ctypes.create_string_buffer(method.encode())
        basis = ctypes.create_string_buffer(basis.encode()) if basis else None
        gcp = _lib.disprs_d3_load_gcp(self.error, self.structure, method, basis)
        _check("d3", self.error)
        try:
            _lib.disprs_d3_set_gcp_controls(
                self.error,
                gcp,
                ctypes.byref(ctypes.c_double(eta)) if eta is not None else None,
                ctypes.byref(ctypes.c_bool(base)) if base is not None else None,
                ctypes.byref(ctypes.c_bool(srb)) if srb is not None else None,
            )
            _check("d3", self.error)
            if parameters is not None:
                _gcp_parameters(self.error, gcp, self.natoms, parameters)
            _lib.disprs_d3_set_gcp_realspace_cutoff(self.error, gcp, cutoff, srb_cutoff)
            energy = ctypes.c_double()
            hessian = (ctypes.c_double * (3 * self.natoms) ** 2)()
            _lib.disprs_d3_get_counterpoise_hessian(
                self.error, self.structure, gcp, ctypes.byref(energy), hessian
            )
            _check("d3", self.error)
            return energy.value, list(hessian)
        finally:
            _delete("disprs_d3_delete_gcp", gcp)

    def close(self):
        _delete("disprs_d3_delete_param", self.param)
        _delete("disprs_d3_delete_model", self.model)
        _delete("disprs_d3_delete_structure", self.structure)
        _delete("disprs_d3_delete_error", self.error)
        self.error = self.structure = self.model = self.param = None

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


def _gcp_parameters(error, gcp, natoms, parameters=None):
    scalar_names = ("sigma", "alpha", "beta", "dmp_scal", "dmp_exp", "rscal", "qscal")
    flag_names = ("damp", "base", "srb")
    array_names = ("zeff", "emiss", "xv", "slater", "rvdw", "rvdw_srb")
    scalar_buffer = (ctypes.c_double * 7)()
    flag_buffer = (ctypes.c_bool * 3)()
    arrays = [
        (kind * count)()
        for kind, count in [
            (ctypes.c_int, natoms),
            *[(ctypes.c_double, natoms)] * 3,
            *[(ctypes.c_double, natoms * natoms)] * 2,
        ]
    ]
    _lib.disprs_d3_get_gcp_parameters(error, gcp, natoms, scalar_buffer, flag_buffer, *arrays)
    _check("d3", error)
    values = dict(zip(scalar_names, scalar_buffer, strict=True))
    values.update(zip(flag_names, flag_buffer, strict=True))
    values.update((name, list(array)) for name, array in zip(array_names, arrays, strict=True))
    if parameters is not None:
        unknown = set(parameters) - values.keys()
        if unknown:
            raise ValueError(f"unknown gCP parameters: {sorted(unknown)}")
        values.update(parameters)
        scalar_buffer = _array(ctypes.c_double, [values[name] for name in scalar_names], 7)
        flag_buffer = _array(ctypes.c_bool, [values[name] for name in flag_names], 3)
        arrays = [
            _array(
                ctypes.c_int if name == "zeff" else ctypes.c_double,
                values[name],
                natoms * natoms if name in ("rvdw", "rvdw_srb") else natoms,
            )
            for name in array_names
        ]
        _lib.disprs_d3_set_gcp_parameters(error, gcp, natoms, scalar_buffer, flag_buffer, *arrays)
        _check("d3", error)
    return values


class D4:
    D4 = 0
    D4S = 1

    @staticmethod
    def damping_parameters(method, *, s9=None):
        """Return named rational damping values, optionally selecting an ATM scaling."""
        method = method.encode()
        scale = ctypes.c_double(s9) if s9 is not None else None
        values = (ctypes.c_double * 6)()
        error = _lib.disprs_d4_new_error()
        try:
            _lib.disprs_d4_get_named_parameters_s9(
                error, method, ctypes.byref(scale) if scale is not None else None, values
            )
            _check("d4", error)
            return dict(zip("s6 s8 s9 a1 a2 alp".split(), values))
        finally:
            _delete("disprs_d4_delete_error", error)

    def __init__(
        self,
        numbers,
        positions,
        method=None,
        charge=0.0,
        model=D4,
        atm=True,
        lattice=None,
        periodic=None,
        *,
        parameters=None,
        ga=3.0,
        gc=2.0,
        wf=6.0,
        charge_model="eeq",
        ghosts=(),
        fixed_charges=None,
    ):
        if charge_model not in ("eeq", "eeqbc"):
            raise ValueError("charge_model must be 'eeq' or 'eeqbc'")
        if (method is None) == (parameters is None):
            raise ValueError("specify exactly one of method or parameters")
        if parameters is not None:
            values = {"s6": 1.0, "s9": 1.0, "alp": 16.0, **parameters}
            if set(values) != {"s6", "s8", "s9", "a1", "a2", "alp"}:
                raise ValueError("parameters require s8, a1, a2; optional keys are s6, s9, alp")
            values = [float(values[key]) for key in ("s6", "s8", "s9", "a1", "a2", "alp")]
            if not atm:
                values[2] = 0.0
        self.natoms = len(numbers)
        self.error = _lib.disprs_d4_new_error()
        self.structure = self.model = self.param = None
        numbers = _array(ctypes.c_int, numbers)
        positions = _array(ctypes.c_double, positions, 3 * self.natoms)
        charge = ctypes.c_double(charge)
        lattice = _array(ctypes.c_double, lattice, 9) if lattice is not None else None
        periodic = _array(ctypes.c_bool, periodic, 3) if periodic is not None else None
        self.structure = _lib.disprs_d4_new_structure(
            self.error,
            self.natoms,
            numbers,
            positions,
            ctypes.byref(charge),
            lattice,
            periodic,
        )
        _check("d4", self.error)
        try:
            self.model = _lib.disprs_d4_new_custom_model(
                self.error, self.structure, model, ga, gc, wf
            )
            _check("d4", self.error)
            self.set_ghosts(ghosts)
            _lib.disprs_d4_set_charge_model(self.error, self.model, int(charge_model == "eeqbc"))
            _check("d4", self.error)
            if fixed_charges is not None:
                self.set_fixed_charges(fixed_charges)
            if parameters is None:
                method = ctypes.create_string_buffer(method.encode())
                self.param = _lib.disprs_d4_load_param(self.error, method, bool(atm))
            else:
                self.param = _lib.disprs_d4_new_rational_damping(self.error, *values)
            _check("d4", self.error)
        except Exception:
            self.close()
            raise

    def update(self, positions, lattice=None):
        positions = _array(ctypes.c_double, positions, 3 * self.natoms)
        lattice = _array(ctypes.c_double, lattice, 9) if lattice is not None else None
        _lib.disprs_d4_update_structure(self.error, self.structure, positions, lattice)
        _check("d4", self.error)

    def set_ghosts(self, indices):
        """Cumulatively exclude zero-based atoms from dispersion, not CN or charges."""
        indices = tuple(indices)
        if any(index < 0 or index >= self.natoms for index in indices):
            raise ValueError("D4 ghost index is out of range")
        indices = _array(ctypes.c_int, indices)
        _lib.disprs_d4_set_model_ghost_index(self.error, self.model, indices, len(indices))
        _check("d4", self.error)

    def set_work_partition(self, part=0, parts=1):
        """Select a zero-based work partition; properties remain full-system values."""
        _lib.disprs_d4_set_model_work_partition(self.error, self.model, part, parts)
        _check("d4", self.error)

    def set_charge_cutoff(self, cutoff=60.0):
        """Set the Ewald real-space splitting hint (bohr) for 1D/2D EEQ."""
        _lib.disprs_d4_set_charge_cutoff(self.error, self.model, cutoff)
        _check("d4", self.error)

    def set_ewald(self, *, rank=0, tolerance=1e-4, kcut=0.0, mesh=0):
        """Enable 3D two-body Fourier dispersion; mesh<0 is direct, >=0 is SPME.

        Rank/tolerance apply to each species-pair reference C6 block. ATM,
        Hessians, pair matrices and real-space pair smoothing are unsupported.
        """
        _lib.disprs_d4_set_model_ewald(self.error, self.model, rank, tolerance, kcut, mesh)
        _check("d4", self.error)

    def set_fixed_charges(self, charges=None):
        """Copy finite atomic charges without rescaling; None restores EEQ/EEQBC."""
        values = _array(ctypes.c_double, charges, self.natoms) if charges is not None else None
        _lib.disprs_d4_set_fixed_charges(
            self.error, self.model, values, self.natoms if values is not None else 0
        )
        _check("d4", self.error)

    def set_realspace_cutoff(
        self, dispersion2=60.0, dispersion3=40.0, coordination=30.0, *, width2=0.0, width3=0.0
    ):
        _lib.disprs_d4_set_realspace_cutoff(
            self.error, self.model, coordination, dispersion2, dispersion3, width2, width3
        )
        _check("d4", self.error)

    def dispersion(self, gradient=False):
        energy = ctypes.c_double()
        grad = (ctypes.c_double * (3 * self.natoms))() if gradient else None
        virial = (ctypes.c_double * 9)() if gradient else None
        _lib.disprs_d4_get_dispersion(
            self.error, self.structure, self.model, self.param, ctypes.byref(energy), grad, virial
        )
        _check("d4", self.error)
        return energy.value, list(grad) if grad else None, list(virial) if virial else None

    def pairwise(self):
        """Return flat two-body and ATM matrices; their sums give the energy."""
        pair2 = (ctypes.c_double * (self.natoms * self.natoms))()
        pair3 = (ctypes.c_double * (self.natoms * self.natoms))()
        _lib.disprs_d4_get_pairwise_dispersion(
            self.error, self.structure, self.model, self.param, pair2, pair3
        )
        _check("d4", self.error)
        return list(pair2), list(pair3)

    def hessian(self):
        """Return energy and a flat analytical Cartesian Hessian at fixed cell (column-major)."""
        energy = self.dispersion()[0]
        hessian = (ctypes.c_double * (3 * self.natoms) ** 2)()
        _lib.disprs_d4_get_dispersion_hessian(
            self.error, self.structure, self.model, self.param, hessian
        )
        _check("d4", self.error)
        return energy, list(hessian)

    def properties(self):
        coordination = (ctypes.c_double * self.natoms)()
        charges = (ctypes.c_double * self.natoms)()
        c6 = (ctypes.c_double * (self.natoms * self.natoms))()
        polarizabilities = (ctypes.c_double * self.natoms)()
        _lib.disprs_d4_get_properties(
            self.error, self.structure, self.model, coordination, charges, c6, polarizabilities
        )
        _check("d4", self.error)
        return list(coordination), list(charges), list(c6), list(polarizabilities)

    def property_response(self):
        """Return {name: (values, Cartesian, strain)}, including charge relaxation."""
        names = ("coordination", "charges", "c6", "polarizabilities")
        sizes = (self.natoms, self.natoms, self.natoms**2, self.natoms)
        values = [(ctypes.c_double * size)() for size in sizes]
        responses = [
            ((ctypes.c_double * (3 * self.natoms * size))(), (ctypes.c_double * (9 * size))())
            for size in sizes
        ]
        _lib.disprs_d4_get_property_response(
            self.error,
            self.structure,
            self.model,
            *values,
            *(array for response in responses for array in response),
        )
        _check("d4", self.error)
        return {
            name: (list(value), list(response[0]), list(response[1]))
            for name, value, response in zip(names, values, responses)
        }

    def close(self):
        _delete("disprs_d4_delete_param", self.param)
        _delete("disprs_d4_delete_model", self.model)
        _delete("disprs_d4_delete_structure", self.structure)
        _delete("disprs_d4_delete_error", self.error)
        self.error = self.structure = self.model = self.param = None

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()
