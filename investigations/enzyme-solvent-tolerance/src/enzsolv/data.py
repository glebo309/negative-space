from pathlib import Path
import re

import pandas as pd


REQUIRED_COLUMNS = {
    "doi",
    "enzyme_name",
    "enzyme_species",
    "mutation",
    "sequence",
    "property",
    "solvent_name",
    "solvent_volume",
    "measured_value",
    "units",
}


def load_measurements(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False)


def _condition_pair(value: object) -> tuple[float | None, float | None]:
    text = "" if pd.isna(value) else str(value)
    number = r"([-+]?\d+(?:\.\d+)?)"
    incubation = re.search(rf"Incubation:\s*{number}", text, flags=re.IGNORECASE)
    assay = re.search(rf"Assay:\s*{number}", text, flags=re.IGNORECASE)
    if incubation or assay:
        return (
            float(incubation.group(1)) if incubation else None,
            float(assay.group(1)) if assay else None,
        )
    unlabelled = re.search(number, text)
    if unlabelled:
        parsed = float(unlabelled.group(1))
        return parsed, parsed
    return None, None


def add_condition_columns(frame: pd.DataFrame) -> pd.DataFrame:
    parsed = frame.copy()
    temperatures = parsed["temperature"].map(_condition_pair)
    ph_values = parsed["ph"].map(_condition_pair)
    parsed["incubation_temperature_c"] = temperatures.str[0]
    parsed["assay_temperature_c"] = temperatures.str[1]
    parsed["incubation_ph"] = ph_values.str[0]
    parsed["assay_ph"] = ph_values.str[1]
    return parsed


def select_primary_endpoint(measurements: pd.DataFrame) -> pd.DataFrame:
    missing = sorted(REQUIRED_COLUMNS.difference(measurements.columns))
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    selected = measurements.copy()
    selected["solvent_volume"] = pd.to_numeric(
        selected["solvent_volume"], errors="coerce"
    )
    selected["measured_value"] = pd.to_numeric(
        selected["measured_value"], errors="coerce"
    )
    sequence_present = selected["sequence"].fillna("").str.strip().ne("")
    comparable = (
        selected["property"].eq("Stability - Incubation")
        & selected["units"].eq("%")
        & selected["mutation"].eq("WT")
        & sequence_present
        & selected["solvent_volume"].notna()
        & selected["measured_value"].notna()
    )
    return selected.loc[comparable].reset_index(drop=True)


def audit_primary_endpoint(selected: pd.DataFrame) -> dict[str, int]:
    enzyme_records = selected[["enzyme_name", "enzyme_species"]].drop_duplicates()
    return {
        "measurements": len(selected),
        "enzyme_records": len(enzyme_records),
        "unique_sequences": selected["sequence"].nunique(),
        "solvents": selected["solvent_name"].nunique(),
        "publications": selected["doi"].nunique(),
    }
