"""The shared drop-directory harvest (cootbridge.harvest) and the moorhen
wrapper's use of it, driven without a database."""

import json
import os

import pytest

django = pytest.importorskip("django")
pytest.importorskip("gemmi", reason="CDataFile stack needs gemmi")

PDB = "ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N\nEND\n"
DICT = """data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
LIG LIG
data_comp_LIG
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
LIG C1 C
"""


@pytest.fixture(scope="module", autouse=True)
def _django_setup():
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "ccp4i2.config.test_settings")
    django.setup()


def _drop(drop_dir, number, text, suffix, meta=None):
    drop_dir.mkdir(parents=True, exist_ok=True)
    path = drop_dir / f"output{number}.{suffix}"
    path.write_text(text)
    if meta is not None:
        path.with_suffix(".meta.json").write_text(json.dumps(meta))
    return path


def test_read_drop_metadata_reads_sidecar_and_tolerates_absence(tmp_path):
    from ccp4i2.cootbridge.harvest import read_drop_metadata

    path = _drop(tmp_path, 0, PDB, "pdb", {"annotation": "first save", "kind": "model"})
    assert read_drop_metadata(path)["annotation"] == "first save"
    assert read_drop_metadata(tmp_path / "output1.pdb") == {}


def test_split_models_and_dictionaries_by_content(tmp_path):
    from ccp4i2.cootbridge.harvest import split_models_and_dictionaries

    model = _drop(tmp_path, 0, PDB, "pdb")
    dictionary = _drop(tmp_path, 1, DICT, "cif")
    models, dicts = split_models_and_dictionaries([model, dictionary])
    assert models == [model]
    assert dicts == [dictionary]


def test_moorhen_harvest_files_outputs_with_sidecar_annotation(tmp_path):
    from ccp4i2.core.CCP4PluginScript import CPluginScript
    from ccp4i2.wrappers.moorhen.script.moorhen import DROP_DIR_NAME, moorhen

    drop_dir = tmp_path / DROP_DIR_NAME
    _drop(drop_dir, 0, PDB, "pdb", {"annotation": "after ligand fit", "kind": "model"})
    _drop(drop_dir, 1, DICT, "cif", {"kind": "dictionary"})

    plugin = moorhen(parent=None, workDirectory=str(tmp_path))
    assert plugin.processOutputFiles() == CPluginScript.SUCCEEDED

    xyzout = plugin.container.outputData.XYZOUT
    assert len(xyzout) == 1
    assert (tmp_path / "XYZOUT_0.pdb").exists()
    metadata = xyzout[0].to_metadata_dict()
    assert metadata.get("annotation") == "after ligand fit"
    assert metadata.get("subType") == 1

    dictout = plugin.container.outputData.DICTOUT
    assert len(dictout) == 1
    assert (tmp_path / "DICTOUT_0.cif").exists()
    assert "Moorhen ligand dictionary" in dictout[0].to_metadata_dict().get("annotation", "")

    program_xml = (tmp_path / "program.xml").read_text()
    assert "<number_output_files>1</number_output_files>" in program_xml


def test_moorhen_harvest_of_empty_session_marks_job_for_deletion(tmp_path):
    from ccp4i2.core.CCP4PluginScript import CPluginScript
    from ccp4i2.wrappers.moorhen.script.moorhen import DROP_DIR_NAME, moorhen

    (tmp_path / DROP_DIR_NAME).mkdir()
    plugin = moorhen(parent=None, workDirectory=str(tmp_path))
    assert plugin.processOutputFiles() == CPluginScript.MARK_TO_DELETE
    assert len(plugin.container.outputData.XYZOUT) == 0


def test_moorhen_start_process_without_database_returns_at_once(tmp_path):
    from ccp4i2.core.CCP4PluginScript import CPluginScript
    from ccp4i2.wrappers.moorhen.script.moorhen import moorhen

    plugin = moorhen(parent=None, workDirectory=str(tmp_path))
    assert plugin.makeCommandAndScript() == CPluginScript.SUCCEEDED
    assert plugin.startProcess() == CPluginScript.SUCCEEDED
