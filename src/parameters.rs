use std::sync::OnceLock;
use toml::de::{DeTable, DeValue};

const PARAMETERS: &str = include_str!("../assets/parameters.toml");
static PARSED_PARAMETERS: OnceLock<DeTable<'static>> = OnceLock::new();

fn d4_name(method: &str) -> String {
    let lower = method.trim().to_ascii_lowercase();
    let mut parts: Vec<_> = lower.split(':').collect();
    parts.sort_unstable();
    let ordered = parts.join(":");
    let name = match ordered.as_str() {
        "gga_c_am05:gga_x_am05" => "am05",
        "gga_c_lyp:gga_x_b88" => "blyp",
        "gga_c_pbe:gga_x_b88" => "bpbe",
        "gga_c_p86:gga_x_b88" => "bp",
        "gga_c_pw91:gga_x_b88" => "bpw",
        "gga_x_lb" => "lb94",
        "gga_c_lyp:gga_x_mpw91" => "mpwlyp",
        "gga_c_pw91:gga_x_mpw91" => "mpwpw",
        "gga_c_lyp:gga_x_optx" => "olyp",
        "gga_c_pbe:gga_x_optx" => "opbe",
        "gga_c_pbe:gga_x_pbe" => "pbe",
        "gga_c_pbe:gga_x_rpbe" => "rpbe",
        "gga_c_pbe:gga_x_pbe_r" => "revpbe",
        "gga_c_pbe_sol:gga_x_pbe_sol" => "pbesol",
        "gga_c_pbe:gga_x_pw86" => "pw86pbe",
        "gga_c_pbe:gga_x_rpw86" => "rpw86pbe",
        "gga_c_pw91:gga_x_pw91" => "pw91",
        "gga_c_p86:gga_x_pw91" => "pwp",
        "mgga_c_tpss:mgga_x_tpss" => "tpss",
        "mgga_c_revtpss:mgga_x_revtpss" => "revtpss",
        "mgga_c_scan:mgga_x_scan" => "scan",
        "mgga_c_rscan:mgga_x_rscan" => "rscan",
        "mgga_c_r2scan:mgga_x_r2scan" => "r2scan",
        "mgga_c_m06:mgga_x_m06" => "m06",
        "mgga_c_m06_l:mgga_x_m06_l" => "m06l",
        "mgga_c_mn12_sx:mgga_c_mn12_sx" => "mn12sx",
        "gga_c_lyp:gga_x_g96" => "glyp",
        "hyb_gga_xc_b3lyp3" | "hyb_gga_xc_b3lyp5" => "b3lyp",
        "hyb_gga_xc_bhandh" | "hyb_gga_xc_bhandhlyp" => "bhlyp",
        "hyb_gga_xc_b3p86" | "hyb_gga_xc_b3p86_nwchem" => "b3p",
        "hyb_gga_xc_b1pw91" => "b1pw",
        "hyb_gga_xc_b3pw91" => "b3pw",
        "hyb_gga_xc_pbeh" => "pbe0",
        "hyb_mgga_xc_b88b95" => "b1b95",
        "hyb_gga_xc_lc_wpbe08_whs" | "hyb_gga_xc_lc_wpbe_whs" | "hyb_gga_xc_lrc_wpbe" => "lcwpbe",
        "hyb_gga_xc_lc_wpbeh_whs" | "hyb_gga_xc_lrc_wpbeh" => "lcwpbeh",
        "mgga_xc_b97m_v" => "b97m",
        "hyb_mgga_xc_wb97m_v" => "wb97m",
        "hyb_gga_xc_wb97x" => "wb97x_2008",
        "hyb_gga_xc_wb97x_v" => "wb97x",
        _ => {
            let stripped = ["xc_hyb_gga_xc_", "hyb_gga_xc_", "hyb_mgga_xc_", "gga_xc_"]
                .iter()
                .find_map(|prefix| lower.strip_prefix(prefix));
            match stripped {
                Some(suffix)
                    if matches!(
                        suffix,
                        "xlyp"
                            | "b97"
                            | "b97_d"
                            | "r2scanh"
                            | "r2scan0"
                            | "r2scan50"
                            | "b1lyp"
                            | "b3lyp"
                            | "o3lyp"
                            | "mpw1pw"
                            | "mpw1lyp"
                            | "pw6b95"
                            | "tpssh"
                            | "tpss0"
                            | "x3lyp"
                            | "cam_b3lyp"
                            | "cam_qtp_01"
                            | "lc_blyp"
                            | "lc_wpbe"
                            | "b2plyp"
                            | "b2gpplyp"
                            | "mpwb1k"
                            | "mpw1b95"
                            | "hse03"
                            | "hse06"
                            | "hse12"
                            | "hse12s"
                            | "hse_sol"
                            | "revtpssh"
                            | "wb97"
                    ) =>
                {
                    suffix
                }
                _ => &lower,
            }
        }
    };
    let normalized = name
        .replace("omega", "w")
        .replace('\u{03c9}', "w")
        .replace('\u{00b2}', "2")
        .replace(['-', '_', '(', ')'], "");
    match normalized.as_str() {
        "bp86" => "bp".into(),
        "b1p86" => "b1p".into(),
        "b3p86" => "b3p".into(),
        "b1pw91" => "b1pw".into(),
        "b3pw91" => "b3pw".into(),
        "mpwpw91" => "mpwpw".into(),
        "mpw1pw91" => "mpw1pw".into(),
        "pw91p86" => "pwp".into(),
        "revdsdpbepbe" => "revdsdpbe".into(),
        "dftb3ob" => "dftb3".into(),
        "lcdftb" => "dftbob2".into(),
        _ => normalized,
    }
}

pub(crate) fn lookup(method: &str, variant: &str) -> Option<&'static DeTable<'static>> {
    let normalize = |value: &str| {
        let lower = value.trim().to_ascii_lowercase();
        if variant.starts_with("d3.") {
            lower.replace(['-', '_', '/'], "")
        } else {
            d4_name(value)
        }
    };
    let name = normalize(method);
    let name = if variant.starts_with("d3.") {
        match name.as_str() {
            "slaterdiracexchange" => "slaterdirac",
            "bp86" => "bp",
            "b88b95" => "b1b95",
            "b1p86" => "b1p",
            "b3p86" => "b3p",
            "bhandhlyp" => "bhlyp",
            "b3lyp5" => "b3lyp",
            "b3lyp3" | "b3lypg" | "dm21" | "dm21m" | "dm21mc" | "dm21mu"
                if matches!(variant, "d3.bj" | "d3.zero") =>
            {
                "b3lyp"
            }
            "pbeh" => "pbe0",
            "dsdtpsstpss" => "dsdtpss",
            "dsdpbepbe" => "dsdpbe",
            "dodpbepbe" => "dodpbe",
            "mpw1pw91" => "mpw1pw",
            "mpwpw91" => "mpwpw",
            "pw91p86" => "pwp",
            "lcomegahpbe" | "lc\u{03c9}hpbe" => "lcwhpbe",
            "\u{03c4}hcth" => "tauhcth",
            "\u{03c4}hcthhyb" => "tauhcthhyb",
            "skala1.1" => "skala1.0",
            _ => &name,
        }
    } else {
        &name
    };
    let parameters = PARSED_PARAMETERS.get_or_init(|| {
        DeTable::parse(PARAMETERS)
            .expect("invalid embedded parameters.toml")
            .into_inner()
    });
    let (model, damping) = variant.split_once('.')?;
    parameters
        .get("parameter")?
        .get_ref()
        .as_table()?
        .iter()
        .find_map(|(key, values)| {
            (normalize(key.get_ref()) == name)
                .then(|| {
                    values
                        .get_ref()
                        .get(model)?
                        .get_ref()
                        .get(damping)?
                        .get_ref()
                        .as_table()
                })
                .flatten()
        })
}

pub(crate) fn parameter(values: &DeTable<'_>, name: &str, default: f64) -> f64 {
    values
        .get(name)
        .map(|value| match value.get_ref() {
            DeValue::Float(number) => number.as_str().parse().expect("invalid damping float"),
            DeValue::Integer(number) => i64::from_str_radix(number.as_str(), number.radix())
                .expect("invalid damping integer") as f64,
            _ => panic!("embedded damping parameter must be numeric"),
        })
        .unwrap_or(default)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_toml_syntax_and_rejects_invalid_values() {
        let parsed = DeTable::parse(
            r#"
            [parameter."quoted.name".d4.bj] # trailing comment
            "s6" = 1
            s8 = 1_000.5e-2 # numeric separator
            doi = "text, s6=99, # not a comment"
        "#,
        )
        .unwrap();
        let values = parsed.get_ref()["parameter"].get_ref()["quoted.name"].get_ref()["d4"]
            .get_ref()["bj"]
            .get_ref()
            .as_table()
            .unwrap();
        assert_eq!(parameter(values, "s6", 0.0), 1.0);
        assert_eq!(parameter(values, "s8", 0.0), 10.005);
        assert_eq!(parameter(values, "absent", 2.0), 2.0);
        assert!(DeTable::parse("s6 = 1\ns6 = 2").is_err());
        assert!(std::panic::catch_unwind(|| parameter(values, "doi", 0.0)).is_err());
    }

    #[test]
    fn d4_libxc_names_resolve_to_the_same_parameters() {
        for (name, libxc) in [
            ("pbe", "gga_x_pbe:gga_c_pbe"),
            ("blyp", "gga_c_lyp:gga_x_b88"),
            ("b3lyp", "hyb_gga_xc_b3lyp"),
            ("wb97x_2008", "hyb_gga_xc_wb97x"),
        ] {
            let canonical = lookup(name, "d4.bj").unwrap();
            let alias = lookup(libxc, "d4.bj").expect(libxc);
            for key in ["s6", "s8", "s9", "a1", "a2", "alp"] {
                assert_eq!(parameter(canonical, key, 0.0), parameter(alias, key, 0.0));
            }
        }
        assert!(lookup("gga_x_unknown:gga_c_pbe", "d4.bj").is_none());
    }

    #[test]
    fn reads_quoted_names_and_last_numeric_field() {
        let values = lookup("WB97X-3C", "d4.bj").unwrap();
        assert_eq!(parameter(values, "a1", 0.0), 0.2464);
        assert_eq!(parameter(values, "alp", 0.0), 16.0);
        assert_eq!(parameter(values, "absent", 2.0), 2.0);
        assert!(lookup("wb97x-3c", "d3.bj").is_none());
        assert!(lookup("unknown", "d4.bj").is_none());
        assert!(lookup("skala-1.0", "d3.bj").is_some());
    }
}
