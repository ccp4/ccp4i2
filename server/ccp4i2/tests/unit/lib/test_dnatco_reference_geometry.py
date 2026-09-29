"""The PDB-NA-RS reference geometry the DNATCO report shows beside a flagged value.

The table is a data file shipped next to the wrapper, so these tests pin both
that it is found and that its keys line up with the names DNATCO's own NAVAL
output uses ("C1'-C2'" in the JSON, "C1'_C2'" in the table).
"""
import json

from ccp4i2.wrappers.dnatco.script import dnatco_data

# The standard nucleotides the reference set covers.
EXPECTED_COMPOUNDS = {"A", "C", "G", "U", "DA", "DC", "DG", "DT", "DU"}


def test_reference_table_ships_next_to_the_wrapper():
    values = dnatco_data.reference_geometry_values()
    assert set(values) == EXPECTED_COMPOUNDS
    for compound, terms in values.items():
        assert terms, f"{compound} has no reference terms"
        for name, value in terms.items():
            assert isinstance(value, float), f"{compound}.{name} is {type(value).__name__}"


def test_lookup_translates_naval_names_to_table_keys():
    # NAVAL writes "C1'-C2'"; the table keys it "C1'_C2'".
    assert dnatco_data.reference_value("A", "C1'-C2'") == 1.528209
    # "r2" marks an atom of the preceding nucleotide, and survives the mapping.
    assert dnatco_data.reference_value("DG", "C3'r2-O3'r2-P") is not None


def test_lookup_is_case_insensitive_and_tolerates_whitespace():
    assert dnatco_data.reference_value("a", " C1'-C2' ") == 1.528209


def test_lookup_returns_none_rather_than_guessing():
    assert dnatco_data.reference_value("PSU", "C1'-C2'") is None   # modified base, not covered
    assert dnatco_data.reference_value("A", "NO-SUCH-TERM") is None
    assert dnatco_data.reference_value("", "C1'-C2'") is None
    assert dnatco_data.reference_value("A", None) is None


def test_concerned_items_carry_their_reference():
    entries = [{
        "compound": "A", "authChain": "A", "authSeqId": 5,
        "details": [
            {"name": "C1'-C2'", "value": 1.61, "naval_tier": "Of Concern",
             "prosco": 0.1, "pGroup": "Rare"},
            {"name": "NOT-A-TERM", "value": 9.9, "naval_tier": "Of Concern",
             "prosco": 0.2, "pGroup": "Rare"},
        ],
    }]
    by_name = {item["name"]: item for item in dnatco_data.concerned_items(entries)}
    assert by_name["C1'-C2'"]["reference"] == 1.528209
    # No reference is None, not a placeholder string: the report formats this
    # number, and a string here would raise instead of drawing the report.
    assert by_name["NOT-A-TERM"]["reference"] is None


def test_a_missing_reference_table_does_not_break_the_report(tmp_path, monkeypatch):
    monkeypatch.setattr(dnatco_data, "REFERENCE_GEOMETRY_FILE", "does-not-exist.json")
    dnatco_data.reference_geometry_values.cache_clear()
    try:
        assert dnatco_data.reference_geometry_values() == {}
        assert dnatco_data.reference_value("A", "C1'-C2'") is None
    finally:
        dnatco_data.reference_geometry_values.cache_clear()
