use clap::{Args, CommandFactory, Parser, Subcommand, ValueHint};
use serde_json::{json, Value};
mod qcschema;
use std::{
    fs,
    io::{self, Read, Write},
    path::{Path, PathBuf},
};

mod disprs_ext;

#[cfg(meson)]
mod build_config {
    include!("build_config.rs");
}
#[cfg(not(meson))]
mod build_config {
    pub const VERSION: &str = env!("CARGO_PKG_VERSION");
    #[cfg(not(feature = "static"))]
    pub const LIBRARY: Option<&str> = option_env!("DISPRS_DEFAULT_LIBRARY");
}

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const BOHR: f64 = 0.529177210903;
const NAMES: [&str; 16] = [
    "d3(zero)",
    "d3(bj)",
    "d3(mzero)",
    "d3(mbj)",
    "d3(op)",
    "d3(cso)",
    "d3(z)",
    "d4(eeq)",
    "d4(eeqbc)",
    "d4s(eeq)",
    "d4s(eeqbc)",
    "d3s(zero)",
    "d3s(bj)",
    "gcp",
    "srb",
    "gcp+srb",
];
const ELEMENTS: &str = "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og";

fn basis_name(value: &str) -> std::result::Result<String, String> {
    let normalized = value.to_ascii_lowercase().replace('-', "");
    let normalized = if normalized == "def2svp" {
        "svp".to_owned()
    } else {
        normalized
    };
    if "sv sv(p) def2sv(p) sv_p def2sv_p svx svp minis 631gd 631gs tz def2tzvp deftzvp def1tzvp ccdz ccpvdz accdz augccpvdz accpvdz pobtz pobtzvp minix hf3c gcore 2g twog fitg dzp dz msvp def2msvp lanl pbeh3c hse3c mtzvp def2mtzvp mtzvpp def2mtzvpp r2scan3c"
        .split_whitespace().any(|basis| basis == normalized) {
        Ok(normalized)
    } else {
        Err(format!("unknown gCP/SRB basis: {value}"))
    }
}

#[derive(Debug, Parser)]
#[command(
    name = "disprs",
    version = build_config::VERSION,
    disable_help_subcommand = true,
    about = "Native dispersion and counterpoise corrections",
    after_long_help = "Uses the disprs C API without Python or an upstream numerical fallback. Results are JSON on stdout; errors go to stderr with exit status 2 (success: 0).

Static builds include the numerical library and ignore DISPRS_LIBRARY. In dynamic builds, DISPRS_LIBRARY selects an explicit shared library. Otherwise discovery tries the Meson-configured path, lib/lib64 under the executable's prefix, then the platform loader's search path (LD_LIBRARY_PATH on Linux). DISPRS_NUM_THREADS controls native workers. Help, version and completion require no library.

Examples:
    disprs run water.xyz --name 'd3(op)' --method pbe --gradient
    disprs run cell.json --name 'd4s(eeqbc)' --method b3lyp
    disprs parameters pbe
    disprs completion ./cli-assets
    disprs help ./cli-assets --format roff
    disprs help ./cli-docs --format rst"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Evaluate a dispersion or counterpoise correction
    #[command(
        after_long_help = "Quote selectors containing parentheses. --name is required. D3/D3S default to no ATM; D4/D4S default to ATM. D3S supports elements 1-94, reuses D3 damping parameters, and requires real-space summation. Named Z fits are unavailable; use --parameters for d3(z).

gcp selects only fitted overlap correction; srb selects fitted short-range terms (HF-3c base or B97-3c SRB); gcp+srb requires both. Unavailable terms are errors, not zero substitutes. Ordinary gCP functional names require --basis.

Input: XYZ records contain an element symbol or atomic number and three coordinates, defaulting to angstrom. The comment line is ignored, including extended XYZ metadata. Every JSON input must be a top-level QCSchema Molecule with schema_name=qcschema_molecule, schema_version=2 or 3, title-case symbols, and flat geometry in bohr. Optional molecular_charge defaults to zero; real=false marks ghost atoms (D3/D4 only). Unknown fields and malformed metadata are rejected, even with validated=true. Legacy numbers/positions/units/charge and JSON cell fields are not accepted. Use --lattice and --periodic for cells (bohr for JSON); a cell defaults to full periodicity. --charge overrides molecular_charge only after validation. JSON cannot use --units angstrom.

--parameters also reads a QCSchema Molecule; damping values belong in its extras.disprs.parameters object. Its molecular data is validated but does not replace the positional geometry input.
Explicit parameter keys in extras.disprs.parameters:
    D3 zero/mzero: s6, s8, s9, rs6, rs8, alp (plus bet for mzero).
    D3 BJ/modified BJ/OP: s6, s8, s9, a1, a2, alp (plus bet for OP).
    D3 CSO: s6, s9, a1, a2, a3, a4, alp. D3 Z: s6, s8, s9, a1, alp.
    D4/D4S: s6, s8, s9, a1, a2, alp.
Defaults: s6=1, s9=1, alp=14 (16 for D4), rs8=1; CSO also defaults a2=2.5, a3=0, a4=6.25. Other keys are required. Disabled ATM overrides s9 with zero. Unknown keys are rejected. Not applicable to gCP/SRB.

FFT requires d3(zero), d3(bj) or D4/D4S, full 3D periodicity and disabled ATM (--no-atm for D4/D4S); Hessians and pairwise output are unavailable. Mesh, rank and tolerance require --fft. D4 rank applies per species-pair reference block. This is reciprocal-space dispersion, not an FFT of the C6 coefficient matrix. Unsupported combinations are rejected, never silently downgraded. Check convergence for your structure.

Output units: energy/strain virial in hartree, gradient in hartree/bohr, Hessian in hartree/bohr squared; C6 and polarizabilities in atomic units, CN dimensionless, charge in elementary-charge units. Gradients are atom-major xyz; matrices/virials are flattened column-major; Hessian gradient components are contiguous per displacement. Symmetric pair matrices sum to their corresponding dispersion energy."
    )]
    Run(Box<Run>),
    /// Query installed-library parameterizations for a method
    #[command(
        after_long_help = "Queries all selectors unless --name limits the query. Returns available fits and reasons for unavailable variants. Unknown methods with no available fits are errors. gCP/SRB entries contain global scalars and flags; element arrays depend on geometry. No fit table is copied into the CLI."
    )]
    Parameters(Query),
    /// Generate Bash, Zsh and Fish completions without loading the library
    #[command(
        after_long_help = "Writes disprs (Bash), _disprs (Zsh), and disprs.fish (Fish). Existing generated files are overwritten. Source the Bash/Fish file, or add the Zsh directory to fpath before compinit."
    )]
    Completion {
        /// Directory for generated files (created if necessary)
        #[arg(value_hint = ValueHint::DirPath)]
        output_dir: PathBuf,
    },
    /// Generate command manuals without loading the library
    #[command(
        after_long_help = "Writes one manual per command, including disprs itself. Roff uses .1 files for man(1); reStructuredText uses .rst files for Sphinx or Docutils. Existing generated files are overwritten. Use --help on any command for terminal help."
    )]
    Help {
        /// Directory for generated manuals (created if necessary)
        #[arg(value_hint = ValueHint::DirPath)]
        output_dir: PathBuf,
        /// Manual output format
        #[arg(long, value_parser = ["roff", "rst"], default_value = "roff")]
        format: String,
    },
}

#[derive(Debug, Args)]
struct Selection {
    /// Complete correction selector; required for run, optional query filter
    #[arg(long, value_parser = NAMES)]
    name: Option<String>,
    /// Native gCP/SRB basis identifier, for example svp or def2tzvp
    #[arg(long, value_parser = basis_name)]
    basis: Option<String>,
    /// Enable ATM three-body dispersion (not applicable to gCP/SRB)
    #[arg(long, conflicts_with = "no_atm")]
    atm: bool,
    /// Disable ATM three-body dispersion (not applicable to gCP/SRB)
    #[arg(long)]
    no_atm: bool,
}

impl Selection {
    fn atm(&self) -> Option<bool> {
        if self.atm {
            Some(true)
        } else if self.no_atm {
            Some(false)
        } else {
            None
        }
    }
}

#[derive(Debug, Args)]
struct Query {
    /// Functional or composite method, for example pbe, b3lyp or hf-3c
    #[arg(value_hint = ValueHint::Other)]
    method: String,
    #[command(flatten)]
    selection: Selection,
}

#[derive(Debug, Args)]
struct Run {
    /// Single-frame XYZ or JSON geometry; - reads stdin
    #[arg(value_hint = ValueHint::FilePath)]
    input: String,
    #[command(flatten)]
    selection: Selection,
    /// Named functional or composite method, for example pbe or hf-3c
    #[arg(
        short,
        long,
        required_unless_present = "parameters",
        conflicts_with = "parameters",
        value_hint = ValueHint::Other
    )]
    method: Option<String>,
    /// QCSchema Molecule file with extras.disprs.parameters instead of a named method
    #[arg(long)]
    parameters: Option<PathBuf>,
    /// Override input format detection (.json selects JSON; stdin defaults to XYZ)
    #[arg(long, value_parser = ["xyz", "json"])]
    format: Option<String>,
    /// XYZ coordinate/cell units (default: angstrom); JSON requires bohr
    #[arg(long, value_parser = ["angstrom", "bohr"])]
    units: Option<String>,
    /// Three flattened lattice vectors, in input units
    #[arg(long, num_args = 9, allow_negative_numbers = true)]
    lattice: Option<Vec<f64>>,
    /// Active cell-vector directions; requires a cell
    #[arg(long, value_parser = ["none", "x", "y", "z", "xy", "xz", "yz", "xyz"])]
    periodic: Option<String>,
    /// Finite total charge, D4/D4S only (default: molecular_charge or zero)
    #[arg(long, allow_negative_numbers = true)]
    charge: Option<f64>,
    /// Return Cartesian gradients (not forces) and strain virial
    #[arg(long)]
    gradient: bool,
    /// Return analytical Cartesian Hessian at fixed cell for D3/gCP/D4/D4S
    #[arg(long)]
    hessian: bool,
    /// Return two-body and ATM energy matrices (D3/D4/D4S real space only)
    #[arg(long)]
    pairwise: bool,
    /// Return CN and C6; D4/D4S also return charges and polarizabilities
    #[arg(long)]
    properties: bool,
    /// D3 zero/BJ or D4/D4S two-body 3D Fourier dispersion with low-rank C6
    #[arg(long)]
    fft: bool,
    /// Mesh size per axis: power of two, zero for automatic, -1 for direct Ewald
    #[arg(long, allow_negative_numbers = true)]
    mesh: Option<i32>,
    /// C6 factorization rank: zero for automatic
    #[arg(long, allow_negative_numbers = true)]
    rank: Option<i32>,
    /// Positive finite C6 factorization tolerance (default: 1e-4)
    #[arg(long, allow_negative_numbers = true)]
    tolerance: Option<f64>,
}

impl Run {
    fn name(&self) -> Result<&str> {
        self.selection
            .name
            .as_deref()
            .ok_or("run requires --name".into())
    }

    fn validate(&self, geometry: &Geometry) -> Result<()> {
        let name = self.name()?;
        let gcp = NAMES[13..].contains(&name);
        if gcp && !geometry.ghosts.is_empty() {
            return Err("QCSchema ghost atoms are not supported for gCP/SRB".into());
        }
        if !gcp && self.selection.basis.is_some() {
            return Err("--basis is only valid for gcp, srb and gcp+srb".into());
        }
        if gcp
            && (self.parameters.is_some()
                || self.selection.atm().is_some()
                || self.pairwise
                || self.properties)
        {
            return Err(
                "gCP/SRB does not support damping parameters, ATM, pairwise or properties".into(),
            );
        }
        if !name.starts_with("d4") && (geometry.charge != 0.0 || self.charge.is_some()) {
            return Err("charge is only supported for D4/D4S".into());
        }
        if !self.fft && (self.mesh.is_some() || self.rank.is_some() || self.tolerance.is_some()) {
            return Err("--mesh, --rank and --tolerance require --fft".into());
        }
        if self.fft {
            if !(name.starts_with("d4") || ["d3(zero)", "d3(bj)"].contains(&name))
                || !geometry.periodic.iter().all(|flag| *flag)
            {
                return Err(
                    "--fft requires d3(zero), d3(bj) or D4/D4S with full 3D periodicity".into(),
                );
            }
            if self.selection.atm().unwrap_or(name.starts_with("d4"))
                || self.hessian
                || self.pairwise
            {
                return Err("--fft does not support ATM, Hessians or pairwise output".into());
            }
            let mesh = self.mesh.unwrap_or(0);
            if mesh < -1 || (mesh > 0 && !(mesh as u32).is_power_of_two()) {
                return Err(
                    "--mesh must be -1 (direct Ewald), zero or a positive power of two".into(),
                );
            }
            if self.rank.unwrap_or(0) < 0 {
                return Err("--rank must fit a nonnegative 32-bit integer".into());
            }
            let tolerance = self.tolerance.unwrap_or(1e-4);
            if !tolerance.is_finite() || tolerance <= 0.0 {
                return Err("--tolerance must be finite and positive".into());
            }
        }
        Ok(())
    }
}

#[derive(Debug)]
struct Geometry {
    numbers: Vec<i32>,
    ghosts: Vec<i32>,
    positions: Vec<f64>,
    lattice: Option<Vec<f64>>,
    periodic: [bool; 3],
    charge: f64,
}

impl Geometry {
    fn read(options: &Run) -> Result<Self> {
        let text = if options.input == "-" {
            let mut text = String::new();
            io::stdin().read_to_string(&mut text)?;
            text
        } else {
            fs::read_to_string(&options.input)?
        };
        Self::parse(options, &text)
    }

    fn parse(options: &Run, text: &str) -> Result<Self> {
        let format = options
            .format
            .as_deref()
            .unwrap_or(if options.input.ends_with(".json") {
                "json"
            } else {
                "xyz"
            });
        let data: Value = if format == "json" {
            qcschema::parse(text)?
        } else {
            let mut lines = text.lines();
            let count: usize = lines
                .next()
                .ok_or("XYZ requires an atom count")?
                .trim()
                .parse()?;
            lines.next().ok_or("XYZ requires a comment line")?;
            let mut numbers = Vec::new();
            let mut positions = Vec::new();
            for line in lines {
                let fields: Vec<_> = line.split_whitespace().collect();
                if fields.len() != 4 {
                    return Err("XYZ atom records require an element and three coordinates".into());
                }
                let number = match ELEMENTS
                    .split_whitespace()
                    .position(|symbol| symbol.eq_ignore_ascii_case(fields[0]))
                {
                    Some(index) => index as i32 + 1,
                    None => fields[0].parse()?,
                };
                numbers.push(number);
                for field in &fields[1..] {
                    let value: f64 = field.parse()?;
                    if !value.is_finite() {
                        return Err("positions must be finite".into());
                    }
                    positions.push(value);
                }
            }
            if count == 0 || numbers.len() != count {
                return Err(
                    "XYZ must contain exactly one frame with the declared atom count".into(),
                );
            }
            json!({"numbers": numbers, "positions": positions})
        };
        let numbers: Vec<i32> = if format == "json" {
            data["symbols"]
                .as_array()
                .unwrap()
                .iter()
                .map(|symbol| {
                    ELEMENTS
                        .split_whitespace()
                        .position(|element| Some(element) == symbol.as_str())
                        .unwrap() as i32
                        + 1
                })
                .collect()
        } else {
            serde_json::from_value(data["numbers"].clone())?
        };
        if numbers.is_empty()
            || numbers.len() > i32::MAX as usize
            || numbers.iter().any(|number| !(1..=118).contains(number))
        {
            return Err("numbers must be a nonempty list of atomic numbers 1 through 118".into());
        }
        if format == "json"
            && options
                .units
                .as_deref()
                .is_some_and(|units| units != "bohr")
        {
            return Err("QCSchema geometry and CLI lattice vectors must be in bohr; --units angstrom is not allowed".into());
        }
        let units =
            options
                .units
                .as_deref()
                .unwrap_or(if format == "json" { "bohr" } else { "angstrom" });
        let scale = match units {
            "bohr" => 1.0,
            "angstrom" => 1.0 / BOHR,
            _ => return Err("units must be angstrom or bohr".into()),
        };
        let mut positions: Vec<f64> = serde_json::from_value(
            data[if format == "json" {
                "geometry"
            } else {
                "positions"
            }]
            .clone(),
        )?;
        let mut lattice = options.lattice.clone();
        if positions.len() != 3 * numbers.len()
            || lattice.as_ref().is_some_and(|cell| cell.len() != 9)
        {
            return Err("positions require 3N values and lattice requires 9 values".into());
        }
        for value in positions.iter_mut().chain(lattice.iter_mut().flatten()) {
            *value *= scale;
            if !value.is_finite() {
                return Err("coordinates and cell must be finite".into());
            }
        }
        let periodic = if let Some(axes) = &options.periodic {
            ['x', 'y', 'z'].map(|axis| axes.contains(axis))
        } else {
            [lattice.is_some(); 3]
        };
        if periodic.iter().any(|flag| *flag) && lattice.is_none() {
            return Err("periodicity requires lattice vectors".into());
        }
        let charge = match options.charge {
            Some(value) => value,
            None => match data.get("molecular_charge") {
                Some(value) => value.as_f64().ok_or("charge must be finite")?,
                None => 0.0,
            },
        };
        if !charge.is_finite() {
            return Err("charge must be finite".into());
        }
        Ok(Self {
            numbers,
            ghosts: data
                .get("real")
                .and_then(Value::as_array)
                .map_or_else(Vec::new, |real| {
                    real.iter()
                        .enumerate()
                        .filter(|(_, real)| **real == false)
                        .map(|(atom, _)| atom as i32)
                        .collect()
                }),
            positions,
            lattice,
            periodic,
            charge,
        })
    }
}

fn completion(directory: &Path) -> Result<()> {
    fs::create_dir_all(directory)?;
    for shell in [
        clap_complete::Shell::Bash,
        clap_complete::Shell::Zsh,
        clap_complete::Shell::Fish,
    ] {
        let path = clap_complete::generate_to(shell, &mut Cli::command(), "disprs", directory)?;
        if shell == clap_complete::Shell::Bash {
            fs::rename(&path, directory.join("disprs"))?;
        }
        if shell == clap_complete::Shell::Fish {
            // ponytail: quote selector tokens until clap_complete escapes Fish's second parsing pass.
            let mut text = fs::read_to_string(&path)?;
            for name in NAMES {
                text = text.replace(&format!("{name}\\t"), &format!("'{name}'\\t"));
            }
            fs::write(path, text)?;
        }
    }
    Ok(())
}

fn help(directory: &Path, format: &str) -> Result<()> {
    fs::create_dir_all(directory)?;
    if format == "roff" {
        clap_mangen::generate_to(Cli::command(), directory)?;
        return Ok(());
    }
    let mut root = Cli::command();
    root.build();
    for mut command in std::iter::once(root.clone()).chain(root.get_subcommands().cloned()) {
        let name = if command.get_name() == "disprs" {
            "disprs".to_owned()
        } else {
            format!("disprs-{}", command.get_name())
        };
        let title = format!("{name}(1)");
        let mut text = format!("{title}\n{}\n\n", "=".repeat(title.len()));
        if let Some(about) = command.get_long_about().or(command.get_about()) {
            text.push_str(&format!("{about}\n\n"));
        }
        text.push_str(&format!(
            "Synopsis\n--------\n\n.. code-block:: console\n\n   {}\n\n",
            command
                .render_usage()
                .to_string()
                .trim_start_matches("Usage: ")
                .replace('\n', "\n   ")
        ));
        text.push_str("Arguments And Options\n---------------------\n\n");
        for argument in command
            .get_arguments()
            .filter(|argument| !argument.is_hide_set())
        {
            if let Some(short) = argument
                .get_short()
                .filter(|_| argument.get_long().is_some())
            {
                text.push_str(&format!("``-{short}``, "));
            }
            text.push_str(&format!("``{argument}``\n"));
            if let Some(help) = argument.get_long_help().or(argument.get_help()) {
                for line in help.to_string().lines() {
                    text.push_str(&format!("   {line}\n"));
                }
            }
            let values = argument.get_possible_values();
            if !values.is_empty() {
                let values: Vec<_> = values
                    .iter()
                    .filter(|value| !value.is_hide_set())
                    .map(|value| format!("``{}``", value.get_name()))
                    .collect();
                text.push_str(&format!("\n   Choices: {}.\n", values.join(", ")));
            }
            for default in argument.get_default_values() {
                text.push_str(&format!(
                    "\n   Default: ``{}``.\n",
                    default.to_string_lossy()
                ));
            }
            text.push('\n');
        }
        if command.has_subcommands() {
            text.push_str("Commands\n--------\n\n");
            for subcommand in command.get_subcommands() {
                text.push_str(&format!(
                    "* :doc:`disprs-{}(1) <disprs-{}>`\n",
                    subcommand.get_name(),
                    subcommand.get_name()
                ));
            }
            text.push('\n');
        }
        if let Some(details) = command.get_after_long_help().or(command.get_after_help()) {
            text.push_str("Details\n-------\n\n");
            for line in details.to_string().lines() {
                text.push_str(&format!("| {line}\n"));
            }
        }
        fs::write(directory.join(format!("{name}.rst")), text)?;
    }
    Ok(())
}

fn execute(command: Command) -> Result<Value> {
    match command {
        Command::Run(options) => {
            let geometry = Geometry::read(&options)?;
            options.validate(&geometry)?;
            disprs_ext::run(&disprs_ext::Api::load()?, &options, &geometry)
        }
        Command::Parameters(options) => disprs_ext::query(&disprs_ext::Api::load()?, &options),
        Command::Completion { output_dir } => {
            completion(&output_dir)?;
            Ok(json!({"directory": output_dir}))
        }
        Command::Help { output_dir, format } => {
            help(&output_dir, &format)?;
            Ok(json!({"directory": output_dir}))
        }
    }
}

fn main() {
    let cli = Cli::parse();
    let result = execute(cli.command).and_then(|result| {
        let output = serde_json::to_string_pretty(&result)?;
        writeln!(io::stdout().lock(), "{output}")?;
        Ok(())
    });
    if let Err(error) = result {
        eprintln!("disprs: error: {error}");
        std::process::exit(2);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_options(extra: &[&str]) -> Box<Run> {
        let mut args = vec!["disprs", "run", "-", "--method", "pbe"];
        args.extend(extra);
        match Cli::try_parse_from(args).unwrap().command {
            Command::Run(options) => options,
            _ => unreachable!(),
        }
    }

    #[test]
    fn generated_assets_follow_cli() {
        let directory =
            std::env::temp_dir().join(format!("disprs-cli-assets-{}", std::process::id()));
        completion(&directory).unwrap();
        assert!(!directory.join("disprs.1").exists());
        help(&directory, "roff").unwrap();
        help(&directory, "rst").unwrap();
        let mut command = Cli::command();
        command.build();
        let root_manual = fs::read_to_string(directory.join("disprs.rst")).unwrap();
        for subcommand in command.get_subcommands() {
            let name = subcommand.get_name();
            assert!(root_manual.contains(&format!("* :doc:`disprs-{name}(1) <disprs-{name}>`")));
            for extension in ["1", "rst"] {
                let manual = fs::read_to_string(
                    directory.join(format!("disprs-{}.{extension}", subcommand.get_name())),
                )
                .unwrap()
                .replace("\\-", "-");
                for argument in subcommand.get_arguments() {
                    if let Some(long) = argument.get_long() {
                        assert!(
                            manual.contains(&format!("--{long}")),
                            "missing {long} in {extension}"
                        );
                    }
                    if let Some(short) = argument.get_short() {
                        assert!(
                            manual.contains(&format!("-{short}")),
                            "missing {short} in {extension}"
                        );
                    }
                }
            }
        }
        for filename in ["disprs.1", "disprs.rst", "disprs", "_disprs", "disprs.fish"] {
            let text = fs::read_to_string(directory.join(filename)).unwrap();
            for subcommand in command.get_subcommands() {
                assert!(text.contains(subcommand.get_name()), "{filename}");
            }
        }
        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn help_requires_no_library() {
        assert_eq!(
            Cli::try_parse_from(["disprs", "--help"])
                .unwrap_err()
                .kind(),
            clap::error::ErrorKind::DisplayHelp
        );
    }

    #[test]
    fn geometry_and_fft_validation() {
        let options = run_options(&["--name", "d3(bj)"]);
        let geometry = Geometry::parse(&options, "2\ncarbon pair\nC 0 0 0\nc 3 0 0\n").unwrap();
        assert_eq!(geometry.numbers, [6, 6]);
        assert_eq!(geometry.positions[3], 3.0 * (1.0 / BOHR));
        assert_eq!(geometry.periodic, [false; 3]);
        let text = r#"{"schema_name":"qcschema_molecule","schema_version":2,"symbols":["C","O"],"geometry":[0,0,0,3,1,0]}"#;
        let options = run_options(&[
            "--name",
            "d3(bj)",
            "--format",
            "json",
            "--fft",
            "--lattice",
            "9",
            "0",
            "0",
            "0",
            "10",
            "0",
            "0",
            "0",
            "11",
        ]);
        let geometry = Geometry::parse(&options, text).unwrap();
        assert_eq!(geometry.positions, [0.0, 0.0, 0.0, 3.0, 1.0, 0.0]);
        assert_eq!(geometry.lattice.as_ref().unwrap().len(), 9);
        options.validate(&geometry).unwrap();
        for extra in [
            vec!["--name", "d4(eeq)", "--fft"],
            vec!["--name", "d3(op)", "--fft"],
            vec!["--name", "d3s(bj)", "--fft"],
            vec!["--name", "d3(bj)", "--fft", "--atm"],
            vec!["--name", "d3(bj)", "--fft", "--mesh", "7"],
            vec!["--name", "d3(bj)", "--mesh", "8"],
            vec!["--name", "gcp", "--properties"],
            vec!["--name", "d3(bj)", "--fft", "--rank", "-1"],
            vec!["--name", "d3(bj)", "--fft", "--tolerance", "NaN"],
        ] {
            assert!(
                run_options(&extra).validate(&geometry).is_err(),
                "{extra:?}"
            );
        }
        for invalid in [
            "{}",
            r#"{"numbers":[6],"positions":[0]}"#,
            r#"{"numbers":[true],"positions":[0,0,0]}"#,
            r#"{"numbers":[6],"positions":[0,0,0],"periodic":[1,1,1]}"#,
            r#"{"numbers":[6],"positions":[0,0,0],"charge":null}"#,
        ] {
            assert!(Geometry::parse(&options, invalid).is_err(), "{invalid}");
        }
    }

    #[test]
    fn native_workflows() {
        let api = disprs_ext::Api::load().unwrap();
        let geometry = Geometry {
            numbers: vec![6, 6],
            ghosts: vec![],
            positions: vec![0.0, 0.0, 0.0, 6.0, 0.0, 0.0],
            lattice: None,
            periodic: [false; 3],
            charge: 0.0,
        };
        for name in NAMES.into_iter().filter(|name| *name != "d3(z)") {
            let method = if NAMES[13..].contains(&name) {
                "hf3c"
            } else {
                "pbe"
            };
            let cli = Cli::try_parse_from([
                "disprs",
                "run",
                "-",
                "--name",
                name,
                "--method",
                method,
                "--gradient",
                "--hessian",
            ])
            .unwrap();
            let Command::Run(options) = cli.command else {
                unreachable!()
            };
            let output = disprs_ext::run(&api, &options, &geometry).unwrap();
            assert!(output["energy"].as_f64().unwrap().is_finite(), "{name}");
            assert_eq!(output["gradient"].as_array().unwrap().len(), 6);
            assert_eq!(output["hessian"].as_array().unwrap().len(), 36);
            if name == "d3(bj)" {
                assert!((output["energy"].as_f64().unwrap() + 0.0005341413931338267).abs() < 1e-14);
            }
        }
        for method in ["pbe", "b3lyp"] {
            let Command::Parameters(options) =
                Cli::try_parse_from(["disprs", "parameters", method])
                    .unwrap()
                    .command
            else {
                unreachable!()
            };
            let output = disprs_ext::query(&api, &options).unwrap();
            assert!(output["available"]["d3(bj)"].is_object());
            assert_eq!(
                output["available"]["d3s(bj)"],
                output["available"]["d3(bj)"]
            );
            assert!(output["available"]["d4s(eeqbc)"].is_object());
            assert!(output["unavailable"]["d3(z)"].is_string());
        }
        let options = run_options(&["--name", "d3(bj)", "--pairwise", "--properties"]);
        let output = disprs_ext::run(&api, &options, &geometry).unwrap();
        let sum: f64 = ["pair2", "pair3"]
            .iter()
            .flat_map(|key| output[key].as_array().unwrap())
            .map(|value| value.as_f64().unwrap())
            .sum();
        assert!((sum - output["energy"].as_f64().unwrap()).abs() < 1e-14);
        assert_eq!(output["properties"]["c6"].as_array().unwrap().len(), 4);
        let mut geometry = geometry;
        geometry.lattice = Some(vec![9.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 11.0]);
        geometry.periodic = [true; 3];
        for name in ["d3(zero)", "d3(bj)", "d4(eeq)", "d4s(eeqbc)"] {
            let options = run_options(&[
                "--name",
                name,
                "--fft",
                "--no-atm",
                "--mesh",
                "16",
                "--gradient",
            ]);
            options.validate(&geometry).unwrap();
            let output = disprs_ext::run(&api, &options, &geometry).unwrap();
            assert_eq!(output["summation"], "fft");
        }
        geometry.charge = 0.5;
        for dimensions in 0..=3 {
            geometry.periodic = std::array::from_fn(|axis| axis < dimensions);
            let options = run_options(&["--name", "d4s(eeqbc)", "--properties"]);
            let output = disprs_ext::run(&api, &options, &geometry).unwrap();
            let charge: f64 = output["properties"]["charges"]
                .as_array()
                .unwrap()
                .iter()
                .map(|value| value.as_f64().unwrap())
                .sum();
            assert!((charge - 0.5).abs() < 1e-12);
        }
    }
}
