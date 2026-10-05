"""MrParse's hits become outputs one each, and a hit without a model is skipped.

The loop used to carry the previous hit's paths over when a hit had no
model file, registering a duplicate under the wrong name (or raising
NameError on the first hit); found while drafting the mrparse judgement.
"""
import json

from ccp4i2.core.tasks import get_plugin_class


def test_a_hit_without_a_model_is_skipped(tmp_path):
    plugin = get_plugin_class("mrparse")(workDirectory=str(tmp_path), name="m")
    plugin.hklin = None
    out = tmp_path / "mrparse_0"
    out.mkdir()
    (out / "1abc_A.pdb").write_text("ATOM\n")
    (out / "homologs.json").write_text(json.dumps([
        {"name": "1xyz_B", "pdb_file": None, "seq_ident": 0.99, "ellg": 0.0},
        {"name": "1abc_A", "pdb_file": "1abc_A.pdb", "seq_ident": 0.95, "ellg": 0.0},
    ]))
    (out / "af_models.json").write_text("[]")
    plugin.processOutputFiles()
    outputs = plugin.container.outputData.XYZOUT
    assert len(outputs) == 1
    assert str(outputs[0].annotation) == "PDB hit: 1abc_A"
    assert str(outputs[0].fullPath).endswith("1abc_A.pdb")
