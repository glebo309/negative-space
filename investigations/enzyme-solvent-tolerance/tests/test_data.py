import numpy as np
import pandas as pd
import pytest

from enzsolv.data import (
    add_condition_columns,
    audit_primary_endpoint,
    select_primary_endpoint,
)


def test_select_primary_endpoint_keeps_comparable_wild_type_rows():
    measurements = pd.DataFrame(
        {
            "doi": ["a", "b", "c", "d", "e"],
            "enzyme_name": ["E1", "E2", "E3", "E4", "E5"],
            "enzyme_species": ["S1", "S2", "S3", "S4", "S5"],
            "mutation": ["WT", "A1V", "WT", "WT", "WT"],
            "sequence": ["AAAA", "BBBB", "CCCC", "", "EEEE"],
            "property": [
                "Stability - Incubation",
                "Stability - Incubation",
                "Activity - Classical",
                "Stability - Incubation",
                "Stability - Incubation",
            ],
            "solvent_name": ["DMSO"] * 5,
            "solvent_volume": [10, 10, 10, 10, "?"],
            "measured_value": [80, 90, 70, 60, 50],
            "units": ["%"] * 5,
        }
    )

    selected = select_primary_endpoint(measurements)

    assert selected["enzyme_name"].tolist() == ["E1"]
    assert selected["solvent_volume"].dtype.kind in "fi"
    assert selected["measured_value"].dtype.kind in "fi"


def test_audit_reports_independent_entities():
    selected = pd.DataFrame(
        {
            "doi": ["a", "a", "b"],
            "enzyme_name": ["E1", "E1", "E2"],
            "enzyme_species": ["S1", "S1", "S2"],
            "sequence": ["AAAA", "AAAA", "BBBB"],
            "solvent_name": ["DMSO", "Methanol", "DMSO"],
        }
    )

    audit = audit_primary_endpoint(selected)

    assert audit == {
        "measurements": 3,
        "enzyme_records": 2,
        "unique_sequences": 2,
        "solvents": 2,
        "publications": 2,
    }


def test_select_primary_endpoint_rejects_missing_schema():
    with pytest.raises(ValueError, match="Missing required columns"):
        select_primary_endpoint(pd.DataFrame({"property": []}))


def test_add_condition_columns_separates_incubation_and_assay_values():
    frame = pd.DataFrame(
        {
            "temperature": [
                "Incubation: 25°C, Assay: 37°C",
                "Incubation: Room temperature, Assay: 60°C",
                "30°C",
            ],
            "ph": [
                "Incubation: 7.5, Assay: 8",
                "Incubation: ?, Assay: 6.5",
                "7",
            ],
        }
    )

    parsed = add_condition_columns(frame)

    assert parsed["incubation_temperature_c"].tolist()[:1] == [25.0]
    assert np.isnan(parsed.loc[1, "incubation_temperature_c"])
    assert parsed.loc[1, "assay_temperature_c"] == 60.0
    assert parsed.loc[2, "incubation_temperature_c"] == 30.0
    assert parsed.loc[2, "assay_temperature_c"] == 30.0
    assert parsed.loc[0, "incubation_ph"] == 7.5
    assert parsed.loc[0, "assay_ph"] == 8.0
    assert np.isnan(parsed.loc[1, "incubation_ph"])
