use crate::{Geometry, Query, Result, Run, NAMES};
#[cfg(not(feature = "static"))]
use libloading::Library;
use serde_json::{json, Map, Value};
#[cfg(not(feature = "static"))]
use std::path::PathBuf;
use std::{
    ffi::{c_char, c_void, CStr, CString},
    ptr::{null, null_mut},
};

type Handle = *mut c_void;

macro_rules! api {
    ($($name:ident($($argument:ty),*) $(-> $return:ty)?;)+) => {
        pub struct Api {
            $($name: unsafe extern "C" fn($($argument),*) $(-> $return)?,)+
            #[cfg(not(feature = "static"))]
            _library: Library,
        }
        impl Api {
            #[cfg(not(feature = "static"))]
            unsafe fn from_library(library: Library) -> Result<Self> {
                Ok(Self {
                    $($name: *library.get(concat!("disprs_", stringify!($name), "\0").as_bytes())?,)+
                    _library: library,
                })
            }

            #[cfg(feature = "static")]
            pub fn load() -> Result<Self> {
                #[link(name = "disprs", kind = "static")]
                extern "C" {
                    $(#[link_name = concat!("disprs_", stringify!($name))]
                    fn $name($(_: $argument),*) $(-> $return)?;)+
                }
                Ok(Self { $($name,)+ })
            }
        }
    };
}

api! {
    d3_new_error() -> Handle;
    d3_check_error(Handle) -> i32;
    d3_get_error(Handle, *mut c_char, *const i32);
    d3_delete_error(*mut Handle);
    d3_new_structure(Handle, i32, *const i32, *const f64, *const f64, *const bool) -> Handle;
    d3_delete_structure(*mut Handle);
    d3_new_model_kind(Handle, Handle, i32) -> Handle;
    d3_delete_model(*mut Handle);
    d3_load_param(Handle, i32, *const c_char, bool) -> Handle;
    d3_get_named_parameters(Handle, i32, *const c_char, *mut f64);
    d3_new_zero_damping(Handle, f64, f64, f64, f64, f64, f64) -> Handle;
    d3_new_rational_damping(Handle, f64, f64, f64, f64, f64, f64) -> Handle;
    d3_new_mzero_damping(Handle, f64, f64, f64, f64, f64, f64, f64) -> Handle;
    d3_new_mrational_damping(Handle, f64, f64, f64, f64, f64, f64) -> Handle;
    d3_new_optimizedpower_damping(Handle, f64, f64, f64, f64, f64, f64, f64) -> Handle;
    d3_new_cso_damping(Handle, f64, f64, f64, f64, f64, f64, f64) -> Handle;
    d3_new_z_damping(Handle, f64, f64, f64, f64, f64) -> Handle;
    d3_delete_param(*mut Handle);
    d3_set_model_ewald(Handle, Handle, i32, f64, f64, i32);
    d3_set_model_ghost_index(Handle, Handle, *const i32, i32);
    d3_get_dispersion(Handle, Handle, Handle, Handle, *mut f64, *mut f64, *mut f64);
    d3_get_dispersion_hessian(Handle, Handle, Handle, Handle, *mut f64, *mut f64);
    d3_get_pairwise_dispersion(Handle, Handle, Handle, Handle, *mut f64, *mut f64);
    d3_get_properties(Handle, Handle, Handle, *mut f64, *mut f64);
    d3_load_gcp(Handle, Handle, *const c_char, *const c_char) -> Handle;
    d3_delete_gcp(*mut Handle);
    d3_get_gcp_parameters(Handle, Handle, i32, *mut f64, *mut bool, *mut i32, *mut f64, *mut f64, *mut f64, *mut f64, *mut f64);
    d3_set_gcp_controls(Handle, Handle, *const f64, *const bool, *const bool);
    d3_get_counterpoise(Handle, Handle, Handle, *mut f64, *mut f64, *mut f64);
    d3_get_counterpoise_hessian(Handle, Handle, Handle, *mut f64, *mut f64);
    d4_new_error() -> Handle;
    d4_check_error(Handle) -> i32;
    d4_get_error(Handle, *mut c_char, *const i32);
    d4_delete_error(*mut Handle);
    d4_new_structure(Handle, i32, *const i32, *const f64, *const f64, *const f64, *const bool) -> Handle;
    d4_delete_structure(*mut Handle);
    d4_new_model(Handle, Handle, i32) -> Handle;
    d4_delete_model(*mut Handle);
    d4_set_charge_model(Handle, Handle, i32);
    d4_set_model_ewald(Handle, Handle, i32, f64, f64, i32);
    d4_set_model_ghost_index(Handle, Handle, *const i32, i32);
    d4_load_param(Handle, *const c_char, bool) -> Handle;
    d4_get_named_parameters_s9(Handle, *const c_char, *const f64, *mut f64);
    d4_new_rational_damping(Handle, f64, f64, f64, f64, f64, f64) -> Handle;
    d4_delete_param(*mut Handle);
    d4_get_dispersion(Handle, Handle, Handle, Handle, *mut f64, *mut f64, *mut f64);
    d4_get_pairwise_dispersion(Handle, Handle, Handle, Handle, *mut f64, *mut f64);
    d4_get_properties(Handle, Handle, Handle, *mut f64, *mut f64, *mut f64, *mut f64);
    d4_get_dispersion_hessian(Handle, Handle, Handle, Handle, *mut f64);
}

#[cfg(not(feature = "static"))]
impl Api {
    pub fn load() -> Result<Self> {
        let filename = libloading::library_filename("disprs");
        let candidates = if let Some(path) = std::env::var_os("DISPRS_LIBRARY") {
            vec![PathBuf::from(path)]
        } else {
            let mut paths = Vec::new();
            if let Some(path) = crate::build_config::LIBRARY {
                paths.push(PathBuf::from(path));
            }
            if let Ok(executable) = std::env::current_exe() {
                if let Some(prefix) = executable.parent().and_then(|directory| directory.parent()) {
                    paths.push(prefix.join("lib").join(&filename));
                    paths.push(prefix.join("lib64").join(&filename));
                }
            }
            paths.push(PathBuf::from(filename));
            paths
        };
        let mut errors = Vec::new();
        for path in candidates {
            match unsafe { Library::new(&path) } {
                Ok(library) => {
                    return unsafe { Self::from_library(library) }.map_err(|error| {
                        format!("incompatible disprs library {}: {error}", path.display()).into()
                    })
                }
                Err(error) => errors.push(format!("{}: {error}", path.display())),
            }
        }
        Err(format!(
            "cannot load the installed disprs library; set DISPRS_LIBRARY to its path\n{}",
            errors.join("\n")
        )
        .into())
    }
}

struct Calculation<'a> {
    api: &'a Api,
    d4: bool,
    error: Handle,
    structure: Handle,
    model: Handle,
    param: Handle,
    gcp: Handle,
}

impl Drop for Calculation<'_> {
    fn drop(&mut self) {
        unsafe {
            (self.api.d3_delete_gcp)(&mut self.gcp);
            if self.d4 {
                (self.api.d4_delete_param)(&mut self.param);
                (self.api.d4_delete_model)(&mut self.model);
                (self.api.d4_delete_structure)(&mut self.structure);
                (self.api.d4_delete_error)(&mut self.error);
            } else {
                (self.api.d3_delete_param)(&mut self.param);
                (self.api.d3_delete_model)(&mut self.model);
                (self.api.d3_delete_structure)(&mut self.structure);
                (self.api.d3_delete_error)(&mut self.error);
            }
        }
    }
}

const D3_KEYS: [&str; 7] = [
    "s6 s8 s9 rs6 rs8 alp",
    "s6 s8 s9 a1 a2 alp",
    "s6 s8 s9 rs6 rs8 alp bet",
    "s6 s8 s9 a1 a2 alp",
    "s6 s8 s9 a1 a2 alp bet",
    "s6 s9 a1 a2 a3 a4 alp",
    "s6 s8 s9 a1 alp",
];
const D4_KEYS: &str = "s6 s8 s9 a1 a2 alp";

fn finite(values: &[f64]) -> Result<()> {
    if values.iter().any(|value| !value.is_finite()) {
        return Err("nonfinite numerical result".into());
    }
    Ok(())
}

fn fields(keys: &str, values: &[f64]) -> Result<Value> {
    finite(values)?;
    Ok(Value::Object(
        keys.split_whitespace()
            .zip(values)
            .map(|(key, value)| (key.to_owned(), json!(value)))
            .collect(),
    ))
}

impl<'a> Calculation<'a> {
    fn new(api: &'a Api, d4: bool) -> Result<Self> {
        let error = unsafe {
            if d4 {
                (api.d4_new_error)()
            } else {
                (api.d3_new_error)()
            }
        };
        if error.is_null() {
            return Err("could not allocate native error handle".into());
        }
        Ok(Self {
            api,
            d4,
            error,
            structure: null_mut(),
            model: null_mut(),
            param: null_mut(),
            gcp: null_mut(),
        })
    }

    fn check(&self) -> Result<()> {
        unsafe {
            let failed = if self.d4 {
                (self.api.d4_check_error)(self.error)
            } else {
                (self.api.d3_check_error)(self.error)
            };
            if failed != 0 {
                let mut buffer = [0 as c_char; 512];
                if self.d4 {
                    (self.api.d4_get_error)(self.error, buffer.as_mut_ptr(), &512)
                } else {
                    (self.api.d3_get_error)(self.error, buffer.as_mut_ptr(), &512)
                }
                buffer[511] = 0;
                return Err(CStr::from_ptr(buffer.as_ptr())
                    .to_string_lossy()
                    .into_owned()
                    .into());
            }
        }
        Ok(())
    }

    fn structure(&mut self, geometry: &Geometry) -> Result<()> {
        let count = i32::try_from(geometry.numbers.len())?;
        let lattice = geometry
            .lattice
            .as_ref()
            .map_or(null(), |cell| cell.as_ptr());
        unsafe {
            self.structure = if self.d4 {
                (self.api.d4_new_structure)(
                    self.error,
                    count,
                    geometry.numbers.as_ptr(),
                    geometry.positions.as_ptr(),
                    &geometry.charge,
                    lattice,
                    geometry.periodic.as_ptr(),
                )
            } else {
                (self.api.d3_new_structure)(
                    self.error,
                    count,
                    geometry.numbers.as_ptr(),
                    geometry.positions.as_ptr(),
                    lattice,
                    geometry.periodic.as_ptr(),
                )
            };
        }
        self.check()
    }

    fn named(&self, name: &str, method: &str, atm: Option<bool>) -> Result<Value> {
        let method = CString::new(method)?;
        if self.d4 {
            let mut values = [f64::NAN; 6];
            let scale = atm.map(|value| if value { 1.0 } else { 0.0 });
            unsafe {
                (self.api.d4_get_named_parameters_s9)(
                    self.error,
                    method.as_ptr(),
                    scale.as_ref().map_or(null(), |value| value),
                    values.as_mut_ptr(),
                );
            }
            self.check()?;
            fields(D4_KEYS, &values)
        } else {
            let damping = NAMES
                .iter()
                .position(|candidate| *candidate == name.replacen("d3s(", "d3(", 1))
                .ok_or("invalid D3 selector")?;
            let mut values = [f64::NAN; 9];
            unsafe {
                (self.api.d3_get_named_parameters)(
                    self.error,
                    damping as i32,
                    method.as_ptr(),
                    values.as_mut_ptr(),
                );
            }
            self.check()?;
            values[2] = if atm.unwrap_or(false) { 1.0 } else { 0.0 };
            let mut all = fields("s6 s8 s9 rs6 rs8 a1 a2 alp bet", &values)?;
            all["a3"] = all["rs6"].clone();
            all["a4"] = all["rs8"].clone();
            Ok(Value::Object(
                D3_KEYS[damping]
                    .split_whitespace()
                    .map(|key| (key.to_owned(), all[key].clone()))
                    .collect(),
            ))
        }
    }

    fn load_gcp(
        &mut self,
        name: &str,
        method: &str,
        basis: Option<&str>,
        count: usize,
    ) -> Result<Value> {
        let method = CString::new(method.to_ascii_lowercase().replace('-', ""))?;
        let basis = basis.map(CString::new).transpose()?;
        unsafe {
            self.gcp = (self.api.d3_load_gcp)(
                self.error,
                self.structure,
                method.as_ptr(),
                basis.as_ref().map_or(null(), |value| value.as_ptr()),
            );
        }
        self.check()?;
        let mut scalars = [f64::NAN; 7];
        let mut flags = [false; 3];
        let mut slater = vec![f64::NAN; count];
        unsafe {
            (self.api.d3_get_gcp_parameters)(
                self.error,
                self.gcp,
                count as i32,
                scalars.as_mut_ptr(),
                flags.as_mut_ptr(),
                null_mut(),
                null_mut(),
                null_mut(),
                slater.as_mut_ptr(),
                null_mut(),
                null_mut(),
            );
        }
        self.check()?;
        if name != "srb" && (scalars[0] == 0.0 || !slater.iter().all(|value| *value > 0.0)) {
            return Err("no active gCP fit for this method/basis".into());
        }
        if name != "gcp" && !(flags[1] || flags[2]) {
            return Err("no active SRB fit for this method/basis".into());
        }
        let mut result = fields("sigma alpha beta dmp_scal dmp_exp rscal qscal", &scalars)?;
        for (key, value) in ["damp", "base", "srb"].into_iter().zip(flags) {
            result[key] = json!(value);
        }
        Ok(result)
    }

    fn damping(
        &mut self,
        name: &str,
        options: &Run,
        explicit: Option<&Value>,
        atm: bool,
    ) -> Result<()> {
        let damping = NAMES
            .iter()
            .position(|candidate| *candidate == name.replacen("d3s(", "d3(", 1))
            .ok_or("invalid selector")?;
        if let Some(explicit) = explicit {
            let keys = if self.d4 { D4_KEYS } else { D3_KEYS[damping] };
            let mut defaults =
                json!({"s6": 1.0, "s9": 1.0, "alp": if self.d4 {16.0} else {14.0}, "rs8": 1.0});
            if name == "d3(cso)" {
                defaults["a2"] = json!(2.5);
                defaults["a3"] = json!(0.0);
                defaults["a4"] = json!(6.25);
            }
            let object = explicit
                .as_object()
                .ok_or("parameters must be a JSON object")?;
            if object
                .keys()
                .any(|key| !keys.split_whitespace().any(|allowed| key == allowed))
            {
                return Err(format!("parameters for {name} require {keys}").into());
            }
            let mut values = Vec::new();
            for key in keys.split_whitespace() {
                let value = object
                    .get(key)
                    .unwrap_or(&defaults[key])
                    .as_f64()
                    .ok_or_else(|| format!("parameters for {name} require numerical {key}"))?;
                if !value.is_finite() {
                    return Err("damping parameters must be finite".into());
                }
                values.push(if key == "s9" && !atm { 0.0 } else { value });
            }
            let values = values.as_slice();
            unsafe {
                self.param = if self.d4 {
                    (self.api.d4_new_rational_damping)(
                        self.error, values[0], values[1], values[2], values[3], values[4],
                        values[5],
                    )
                } else {
                    match damping {
                        0 => (self.api.d3_new_zero_damping)(
                            self.error, values[0], values[1], values[2], values[3], values[4],
                            values[5],
                        ),
                        1 => (self.api.d3_new_rational_damping)(
                            self.error, values[0], values[1], values[2], values[3], values[4],
                            values[5],
                        ),
                        2 => (self.api.d3_new_mzero_damping)(
                            self.error, values[0], values[1], values[2], values[3], values[4],
                            values[5], values[6],
                        ),
                        3 => (self.api.d3_new_mrational_damping)(
                            self.error, values[0], values[1], values[2], values[3], values[4],
                            values[5],
                        ),
                        4 => (self.api.d3_new_optimizedpower_damping)(
                            self.error, values[0], values[1], values[2], values[3], values[4],
                            values[5], values[6],
                        ),
                        5 => (self.api.d3_new_cso_damping)(
                            self.error, values[0], values[1], values[2], values[3], values[4],
                            values[5], values[6],
                        ),
                        6 => (self.api.d3_new_z_damping)(
                            self.error, values[0], values[1], values[2], values[3], values[4],
                        ),
                        _ => return Err("invalid D3 selector".into()),
                    }
                };
            }
        } else {
            let method = CString::new(options.method.as_deref().ok_or("--method is required")?)?;
            unsafe {
                self.param = if self.d4 {
                    (self.api.d4_load_param)(self.error, method.as_ptr(), atm)
                } else {
                    (self.api.d3_load_param)(self.error, damping as i32, method.as_ptr(), atm)
                };
            }
        }
        self.check()
    }
}

pub fn query(api: &Api, options: &Query) -> Result<Value> {
    let selection = &options.selection;
    if let Some(name) = selection.name.as_deref() {
        if NAMES[13..].contains(&name) && selection.atm().is_some() {
            return Err("gCP/SRB does not support ATM".into());
        }
        if !NAMES[13..].contains(&name) && selection.basis.is_some() {
            return Err("--basis only applies to gCP/SRB".into());
        }
    }
    let mut available = Map::new();
    let mut unavailable = Map::new();
    for name in NAMES.into_iter().filter(|name| {
        selection
            .name
            .as_deref()
            .is_none_or(|selected| *name == selected)
    }) {
        let mut calculation = Calculation::new(api, name.starts_with("d4"))?;
        let values = if NAMES[13..].contains(&name) {
            calculation.structure(&Geometry {
                numbers: vec![1],
                ghosts: vec![],
                positions: vec![0.0; 3],
                lattice: None,
                periodic: [false; 3],
                charge: 0.0,
            })?;
            calculation.load_gcp(name, &options.method, selection.basis.as_deref(), 1)
        } else {
            calculation.named(name, &options.method, selection.atm())
        };
        match values {
            Ok(values) => {
                available.insert(name.to_owned(), values);
            }
            Err(error) => {
                unavailable.insert(name.to_owned(), json!(error.to_string()));
            }
        }
    }
    if available.is_empty() {
        return Err(format!("no parameters available: {}", Value::Object(unavailable)).into());
    }
    Ok(
        json!({"method": options.method, "basis": selection.basis, "available": available, "unavailable": unavailable,
        "gcp_scope": "global scalars/flags; element arrays depend on the structure"}),
    )
}

fn output_array(result: &mut Value, key: &str, values: Vec<f64>) -> Result<()> {
    finite(&values)?;
    result[key] = json!(values);
    Ok(())
}

pub fn run(api: &Api, options: &Run, geometry: &Geometry) -> Result<Value> {
    let name = options.name()?;
    let d4 = name.starts_with("d4");
    let gcp = NAMES[13..].contains(&name);
    let count = geometry.numbers.len();
    let dimension = count.checked_mul(3).ok_or("too many atoms")?;
    let pairs = count.checked_mul(count).ok_or("too many atoms")?;
    let atm = options.selection.atm().unwrap_or(d4);
    let explicit: Option<Value> = options
        .parameters
        .as_ref()
        .map(|path| -> Result<Value> {
            let molecule = crate::qcschema::parse(&std::fs::read_to_string(path)?)?;
            molecule.pointer("/extras/disprs/parameters").filter(|value| value.is_object()).cloned()
                .ok_or("--parameters requires a QCSchema Molecule with an extras.disprs.parameters object".into())
        })
        .transpose()?;
    let mut calculation = Calculation::new(api, d4)?;
    calculation.structure(geometry)?;
    let error = calculation.error;
    let structure = calculation.structure;
    unsafe {
        if gcp {
            calculation.load_gcp(
                name,
                options.method.as_deref().ok_or("gCP requires --method")?,
                options.selection.basis.as_deref(),
                count,
            )?;
            if name == "gcp" {
                (api.d3_set_gcp_controls)(error, calculation.gcp, null(), &false, &false);
            } else if name == "srb" {
                (api.d3_set_gcp_controls)(error, calculation.gcp, &0.0, null(), null());
            }
            calculation.check()?;
        } else {
            calculation.model = if d4 {
                (api.d4_new_model)(error, structure, i32::from(name.starts_with("d4s")))
            } else {
                (api.d3_new_model_kind)(error, structure, i32::from(name.starts_with("d3s(")))
            };
            calculation.check()?;
            if !geometry.ghosts.is_empty() {
                let set_ghosts = if d4 {
                    api.d4_set_model_ghost_index
                } else {
                    api.d3_set_model_ghost_index
                };
                set_ghosts(
                    error,
                    calculation.model,
                    geometry.ghosts.as_ptr(),
                    geometry.ghosts.len() as i32,
                );
                calculation.check()?;
            }
            if d4 {
                (api.d4_set_charge_model)(
                    error,
                    calculation.model,
                    i32::from(name.contains("eeqbc")),
                );
                calculation.check()?;
            }
            calculation.damping(name, options, explicit.as_ref(), atm)?;
            if options.fft {
                let set_ewald = if d4 {
                    api.d4_set_model_ewald
                } else {
                    api.d3_set_model_ewald
                };
                set_ewald(
                    error,
                    calculation.model,
                    options.rank.unwrap_or(0),
                    options.tolerance.unwrap_or(1e-4),
                    0.0,
                    options.mesh.unwrap_or(0),
                );
                calculation.check()?;
            }
        }
    }
    let model = calculation.model;
    let param = calculation.param;
    let mut energy = f64::NAN;
    let mut gradient = if options.gradient {
        vec![f64::NAN; dimension]
    } else {
        Vec::new()
    };
    let mut virial = if options.gradient {
        vec![f64::NAN; 9]
    } else {
        Vec::new()
    };
    let gradient_ptr = if options.gradient {
        gradient.as_mut_ptr()
    } else {
        null_mut()
    };
    let virial_ptr = if options.gradient {
        virial.as_mut_ptr()
    } else {
        null_mut()
    };
    unsafe {
        if gcp {
            (api.d3_get_counterpoise)(
                error,
                structure,
                calculation.gcp,
                &mut energy,
                gradient_ptr,
                virial_ptr,
            );
        } else if d4 {
            (api.d4_get_dispersion)(
                error,
                structure,
                model,
                param,
                &mut energy,
                gradient_ptr,
                virial_ptr,
            );
        } else {
            (api.d3_get_dispersion)(
                error,
                structure,
                model,
                param,
                &mut energy,
                gradient_ptr,
                virial_ptr,
            );
        }
    }
    calculation.check()?;
    finite(&[energy])?;
    let mut result = json!({"name": name, "method": options.method, "natoms": count, "energy": energy,
        "units": {"energy": "hartree", "gradient": "hartree/bohr", "virial": "hartree", "hessian": "hartree/bohr^2"},
        "summation": if options.fft { "fft" } else { "realspace" }, "periodic": geometry.periodic,
        "charge": geometry.charge, "atm": atm, "basis": options.selection.basis});
    if let Some(explicit) = explicit {
        result["parameters"] = explicit;
    }
    if options.fft {
        result["fft"] = json!({"mesh": options.mesh.unwrap_or(0), "rank": options.rank.unwrap_or(0), "tolerance": options.tolerance.unwrap_or(1e-4)});
    }
    if options.gradient {
        output_array(&mut result, "gradient", gradient)?;
        output_array(&mut result, "virial", virial)?;
    }
    if options.hessian {
        let mut hessian = vec![
            f64::NAN;
            dimension
                .checked_mul(dimension)
                .ok_or("too many Hessian components")?
        ];
        unsafe {
            if gcp {
                (api.d3_get_counterpoise_hessian)(
                    error,
                    structure,
                    calculation.gcp,
                    &mut energy,
                    hessian.as_mut_ptr(),
                );
            } else if d4 {
                (api.d4_get_dispersion_hessian)(
                    error,
                    structure,
                    model,
                    param,
                    hessian.as_mut_ptr(),
                );
            } else {
                (api.d3_get_dispersion_hessian)(
                    error,
                    structure,
                    model,
                    param,
                    &mut energy,
                    hessian.as_mut_ptr(),
                );
            }
        }
        calculation.check()?;
        output_array(&mut result, "hessian", hessian)?;
    }
    if options.pairwise {
        let mut pair2 = vec![f64::NAN; pairs];
        let mut pair3 = vec![f64::NAN; pairs];
        unsafe {
            if d4 {
                (api.d4_get_pairwise_dispersion)(
                    error,
                    structure,
                    model,
                    param,
                    pair2.as_mut_ptr(),
                    pair3.as_mut_ptr(),
                );
            } else {
                (api.d3_get_pairwise_dispersion)(
                    error,
                    structure,
                    model,
                    param,
                    pair2.as_mut_ptr(),
                    pair3.as_mut_ptr(),
                );
            }
        }
        calculation.check()?;
        output_array(&mut result, "pair2", pair2)?;
        output_array(&mut result, "pair3", pair3)?;
    }
    if options.properties {
        let mut coordination = vec![f64::NAN; count];
        let mut c6 = vec![f64::NAN; pairs];
        let mut charges = vec![f64::NAN; count];
        let mut polarizabilities = vec![f64::NAN; count];
        unsafe {
            if d4 {
                (api.d4_get_properties)(
                    error,
                    structure,
                    model,
                    coordination.as_mut_ptr(),
                    charges.as_mut_ptr(),
                    c6.as_mut_ptr(),
                    polarizabilities.as_mut_ptr(),
                );
            } else {
                (api.d3_get_properties)(
                    error,
                    structure,
                    model,
                    coordination.as_mut_ptr(),
                    c6.as_mut_ptr(),
                );
            }
        }
        calculation.check()?;
        let mut properties = json!({});
        output_array(&mut properties, "coordination", coordination)?;
        output_array(&mut properties, "c6", c6)?;
        if d4 {
            output_array(&mut properties, "charges", charges)?;
            output_array(&mut properties, "polarizabilities", polarizabilities)?;
        }
        result["properties"] = properties;
    }
    Ok(result)
}
