"""Which model files a PDB-REDO result zip yields (no network).

PDB-REDO is dropping PDB-format output, and many runs already write mmCIF
only. The wrapper used to look for ``_final.pdb`` and ``_besttls.pdb`` alone,
so such a run finished with maps and no model, and said nothing.
"""

from ccp4i2.wrappers.pdb_redo_api.script.pdb_redo_api import choose_models


def test_mmcif_only_run_yields_both_models():
    names = ["1fjl/1fjl_final.cif", "1fjl/1fjl_final.mtz",
             "1fjl/1fjl_besttls.cif", "1fjl/1fjl_besttls.mtz"]
    assert choose_models(names) == {
        "XYZOUT_FINAL": "1fjl/1fjl_final.cif",
        "XYZOUT_BESTTLS": "1fjl/1fjl_besttls.cif",
    }


def test_mmcif_preferred_when_both_formats_written():
    names = ["r/x_final.pdb", "r/x_final.cif", "r/x_besttls.pdb", "r/x_besttls.cif"]
    assert choose_models(names) == {
        "XYZOUT_FINAL": "r/x_final.cif",
        "XYZOUT_BESTTLS": "r/x_besttls.cif",
    }


def test_pdb_format_still_taken_when_it_is_all_there_is():
    names = ["r/x_final.pdb", "r/x_besttls.pdb"]
    assert choose_models(names) == {
        "XYZOUT_FINAL": "r/x_final.pdb",
        "XYZOUT_BESTTLS": "r/x_besttls.pdb",
    }


def test_no_model_is_an_empty_choice():
    assert choose_models(["r/x_final.mtz", "r/x_final.log", "r/data.json"]) == {}
