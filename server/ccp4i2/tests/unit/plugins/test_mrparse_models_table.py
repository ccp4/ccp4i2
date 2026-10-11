"""mrparse keeps what it knew about each model in program.xml."""
import json
from xml.etree import ElementTree as ET

from ccp4i2.core.tasks import get_plugin_class


def test_each_model_carries_identity_coverage_ellg(tmp_path):
    # Haiku had to align the sequences itself to estimate an identity the
    # wrapper had read from homologs.json and dropped (Eleanor's case).
    plugin = get_plugin_class("mrparse")(workDirectory=str(tmp_path), name="m")
    plugin.hklin = "hklin.mtz"
    plugin._targets = lambda: [("sp|P11802|CDK4_HUMAN", "M" * 303)]
    out = tmp_path / "mrparse_0"
    out.mkdir()
    for name in ("6p8e_B_19-297.pdb", "1blx_A_10-302.pdb"):
        (out / name).write_text("ATOM\n")
    (out / "homologs.json").write_text(json.dumps([
        {"name": "1blx_A", "pdb_file": "1blx_A_10-302.pdb", "seq_ident": 0.43, "ellg": 120.0,
         "score": 88.0, "rmsd": 1.1, "query_start": 10, "query_stop": 302, "length": 293, "ncopies": 1},
        {"name": "6p8e_B", "pdb_file": "6p8e_B_19-297.pdb", "seq_ident": 0.92, "ellg": 551.7,
         "score": 574.0, "rmsd": 0.485, "query_start": 18, "query_stop": 299, "length": 280, "ncopies": 1,
         "resolution": 2.3},
    ]))
    (out / "af_models.json").write_text("[]")
    plugin.processOutputFiles()
    root = ET.parse(plugin.makeFileName("PROGRAMXML")).getroot()
    models = root.findall(".//Models/Model")
    assert [m.get("name") for m in models] == ["6p8e_B", "1blx_A"]   # by eLLG, reflections given
    first = models[0]
    assert (first.get("file"), first.get("source"), first.get("identity"), first.get("ellg")) == \
        ("6p8e_B_19-297.pdb", "PDB", "0.92", "551.7")
    assert first.get("range") == "18-299" and first.get("coverage") == "0.93"   # 282 of 303
    assert first.get("rmsd") == "0.48" and first.get("resolution") == "2.30"
    assert models[1].get("identity") == "0.43"
    # the judgement's xpaths
    assert root.find(".//Models/Model[1]").get("identity") == "0.92"
    assert root.find(".//Models/Model[2]").get("identity") == "0.43"
