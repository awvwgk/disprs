use crate::{Result, ELEMENTS};
use serde_json::Value;

pub fn parse(text: &str) -> Result<Value> {
    let data: Value = serde_json::from_str(text)?;
    validate(&data).map_err(|error| format!("QCSchema Molecule: {error}"))?;
    Ok(data)
}

fn array(value: &Value) -> Result<&Vec<Value>> {
    value.as_array().ok_or("expected an array".into())
}

fn number(value: &Value) -> bool {
    value.as_f64().is_some_and(f64::is_finite)
}

fn integer(value: &Value) -> bool {
    value
        .as_f64()
        .is_some_and(|value| value.is_finite() && value.fract() == 0.0)
}

fn index(value: &Value, count: usize) -> bool {
    integer(value)
        && value
            .as_f64()
            .is_some_and(|value| value >= 0.0 && value < count as f64)
}

fn validate(data: &Value) -> Result<()> {
    let object = data
        .as_object()
        .ok_or("expected a top-level Molecule object")?;
    if data["schema_name"] != "qcschema_molecule"
        || !matches!(data["schema_version"].as_u64(), Some(2 | 3))
    {
        return Err(
            "schema_name must be qcschema_molecule and schema_version must be 2 or 3".into(),
        );
    }
    let symbols = array(&data["symbols"])?;
    let count = symbols.len();
    if count == 0
        || count > i32::MAX as usize
        || symbols.iter().any(|value| {
            !value
                .as_str()
                .is_some_and(|symbol| ELEMENTS.split_whitespace().any(|element| element == symbol))
        })
    {
        return Err(
            "symbols must be a nonempty array of title-case element symbols (Z=1..118)".into(),
        );
    }
    let geometry = array(&data["geometry"])?;
    if geometry.len() != 3 * count || !geometry.iter().all(number) {
        return Err("geometry must contain exactly 3N finite numbers in bohr".into());
    }
    for (key, value) in object {
        if data["schema_version"] == 3
            && value.is_null()
            && matches!(
                key.as_str(),
                "name"
                    | "comment"
                    | "fix_symmetry"
                    | "identifiers"
                    | "masses"
                    | "real"
                    | "atom_labels"
                    | "atomic_numbers"
                    | "mass_numbers"
                    | "connectivity"
                    | "fragments"
                    | "fragment_charges"
                    | "fragment_multiplicities"
            )
        {
            continue;
        }
        let valid = match key.as_str() {
            "schema_name" | "schema_version" | "symbols" | "geometry" | "id" => true,
            "name" | "comment" | "fix_symmetry" => value.is_string(),
            "fix_com" | "fix_orientation" | "validated" => value.is_boolean(),
            "molecular_charge" => number(value),
            "molecular_multiplicity" => {
                number(value)
                    && value.as_f64().unwrap() >= 1.0
                    && (data["schema_version"] == 3 || integer(value))
            }
            "extras" => value.is_object(),
            "real" | "masses" | "atomic_numbers" | "mass_numbers" | "atom_labels" => {
                let values = array(value).map_err(|error| format!("{key}: {error}"))?;
                values.len() == count
                    && values
                        .iter()
                        .enumerate()
                        .all(|(atom, item)| match key.as_str() {
                            "real" => item.is_boolean(),
                            "atom_labels" => item.is_string(),
                            "masses" => number(item) && item.as_f64().unwrap() > 0.0,
                            "mass_numbers" => {
                                integer(item)
                                    && (item.as_f64().unwrap() == -1.0
                                        || item.as_f64().unwrap() > 0.0)
                            }
                            _ => {
                                integer(item)
                                    && ELEMENTS
                                        .split_whitespace()
                                        .position(|symbol| Some(symbol) == symbols[atom].as_str())
                                        .is_some_and(|position| {
                                            item.as_f64() == Some((position + 1) as f64)
                                        })
                            }
                        })
            }
            "fragments" => {
                let mut seen = vec![false; count];
                for fragment in array(value)? {
                    let atoms = array(fragment)?;
                    if atoms.is_empty() {
                        return Err("fragments must not be empty".into());
                    }
                    for atom in atoms {
                        if !index(atom, count) {
                            return Err("fragment atom index out of range".into());
                        }
                        let atom = atom.as_f64().unwrap() as usize;
                        if std::mem::replace(&mut seen[atom], true) {
                            return Err("duplicate fragment atom index".into());
                        }
                    }
                }
                seen.into_iter().all(|present| present)
            }
            "fragment_charges" | "fragment_multiplicities" => {
                let values = array(value)?;
                let fragments = data
                    .get("fragments")
                    .filter(|value| !value.is_null())
                    .map(array)
                    .transpose()?
                    .map_or(1, Vec::len);
                values.len() == fragments
                    && values.iter().all(|item| {
                        number(item)
                            && (key == "fragment_charges"
                                || (item.as_f64().unwrap() >= 1.0
                                    && (data["schema_version"] == 3 || integer(item))))
                    })
            }
            "connectivity" => {
                (data["schema_version"] != 3 || !array(value)?.is_empty())
                    && array(value)?.iter().all(|bond| {
                        bond.as_array().is_some_and(|bond| {
                            bond.len() == 3
                                && index(&bond[0], count)
                                && index(&bond[1], count)
                                && bond[0].as_f64() != bond[1].as_f64()
                                && bond[2]
                                    .as_f64()
                                    .is_some_and(|order| (0.0..=5.0).contains(&order))
                        })
                    })
            }
            "provenance" => value.as_object().is_some_and(|provenance| {
                provenance.get("creator").is_some_and(Value::is_string)
                    && provenance.iter().all(|(key, value)| {
                        !["creator", "version", "routine"].contains(&key.as_str())
                            || value.is_string()
                    })
            }),
            "identifiers" => value.as_object().is_some_and(|identifiers| {
                identifiers.iter().all(|(key, value)| {
                    [
                        "molecule_hash",
                        "molecular_formula",
                        "smiles",
                        "inchi",
                        "inchikey",
                        "canonical_explicit_hydrogen_smiles",
                        "canonical_isomeric_explicit_hydrogen_mapped_smiles",
                        "canonical_isomeric_explicit_hydrogen_smiles",
                        "canonical_isomeric_smiles",
                        "canonical_smiles",
                        "pubchem_cid",
                        "pubchem_sid",
                        "pubchem_conformerid",
                    ]
                    .contains(&key.as_str())
                        && (value.is_string() || (data["schema_version"] == 3 && value.is_null()))
                })
            }),
            _ => {
                return Err(format!("unknown field {key:?}; custom data belongs in extras").into())
            }
        };
        if !valid {
            return Err(format!("invalid {key}").into());
        }
    }
    if let Some(charges) = data
        .get("fragment_charges")
        .filter(|value| !value.is_null())
    {
        let total: f64 = array(charges)?
            .iter()
            .map(|charge| charge.as_f64().unwrap())
            .sum();
        if (total
            - data
                .get("molecular_charge")
                .and_then(Value::as_f64)
                .unwrap_or(0.0))
        .abs()
            > 1e-8
        {
            return Err("fragment_charges must sum to molecular_charge".into());
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn molecule_contract() {
        let valid = json!({"schema_name":"qcschema_molecule", "schema_version":2,
            "symbols":["C", "O"], "geometry":[0,0,0,3,1,0]});
        assert!(parse(&valid.to_string()).is_ok());
        for (key, value) in [
            ("schema_name", json!("qcschema_input")),
            ("schema_version", json!(1)),
            ("symbols", json!(["c", "O"])),
            ("geometry", json!([[0, 0, 0], [3, 1, 0]])),
            ("real", json!([true])),
            ("real", json!([1, 1])),
            ("masses", json!([12, "16"])),
            ("atomic_numbers", json!([6, 7])),
            ("mass_numbers", json!([12.5, 16])),
            ("fragments", json!([[0], [0]])),
            ("fragments", json!([[0]])),
            ("connectivity", json!([[0, 2, 1]])),
            ("connectivity", json!([[0, 1, 6]])),
            ("fragment_charges", json!([1])),
            ("molecular_charge", Value::Null),
            ("molecular_multiplicity", json!(0)),
            ("fix_com", json!(1)),
            ("provenance", json!({"creator":1})),
            ("identifiers", json!({"typo":"x"})),
            ("extras", json!([])),
            ("units", json!("bohr")),
            ("lattice", json!(vec![1; 9])),
        ] {
            let mut invalid = valid.clone();
            invalid[key] = value;
            invalid["validated"] = json!(true);
            assert!(parse(&invalid.to_string()).is_err(), "{invalid}");
        }
        for key in ["schema_name", "schema_version", "symbols", "geometry"] {
            let mut invalid = valid.clone();
            invalid.as_object_mut().unwrap().remove(key);
            assert!(parse(&invalid.to_string()).is_err(), "{key}");
        }
        let mut metadata = valid;
        for (key, value) in [
            ("schema_version", json!(3)),
            ("real", json!([true, false])),
            ("fragments", json!([[0], [1]])),
            ("fragment_charges", json!([0, 0])),
            ("fragment_multiplicities", json!([1, 1])),
            ("connectivity", json!([[0, 1, 1.5]])),
            ("provenance", json!({"creator":"test", "version":"1"})),
            ("extras", json!({"arbitrary":[true, null, 1]})),
        ] {
            metadata[key] = value;
        }
        assert!(parse(&metadata.to_string()).is_ok());
        metadata["connectivity"] = json!([]);
        assert!(parse(&metadata.to_string()).is_err());
        for key in [
            "name",
            "comment",
            "fix_symmetry",
            "identifiers",
            "masses",
            "real",
            "atom_labels",
            "atomic_numbers",
            "mass_numbers",
            "connectivity",
            "fragments",
            "fragment_charges",
            "fragment_multiplicities",
        ] {
            metadata[key] = Value::Null;
        }
        assert!(parse(&metadata.to_string()).is_ok());
    }
}
