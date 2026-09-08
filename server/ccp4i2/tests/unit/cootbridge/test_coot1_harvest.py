"""Harvest metadata survival for the Coot task wrappers.

Regression guard for the annotation/subType loss caused by truncating an
output CList with ``.set(list[:n])``: CList.set() deep-copies items
through CDataFile.get()/set(), which carries only the path fields, so
the gleaner fell back to the bare param name. The wrappers truncate with
pop() instead; these tests pin that the populated metadata survives.
"""

import os

import pytest

django = pytest.importorskip("django")
pytest.importorskip("gemmi", reason="CDataFile stack needs gemmi")


@pytest.fixture(scope="module", autouse=True)
def _django_setup():
    os.environ.setdefault("DJANGO_SETTINGS_MODULE",
                          "ccp4i2.config.test_settings")
    django.setup()


def _populate_one(xyzout, path):
    from ccp4i2.core.CCP4ModelData import CPdbDataFile

    while len(xyzout) < 1:
        xyzout.append(xyzout.makeItem())
    xyzout[0].setFullPath(str(path))
    xyzout[0].annotation.set("Coot output: model.pdb")
    xyzout[0].subType.set(CPdbDataFile.SUBTYPE_MODEL)
    xyzout[0].contentFlag.set(CPdbDataFile.CONTENT_FLAG_PDB)


def test_pop_truncate_preserves_annotation_and_subtype(tmp_path):
    from ccp4i2.wrappers.coot1.script.coot1 import coot1

    plugin = coot1(parent=None, workDirectory=str(tmp_path))
    xyzout = plugin.container.outputData.XYZOUT
    model = tmp_path / "model.pdb"
    model.write_text("")
    _populate_one(xyzout, model)

    while len(xyzout) > 1:  # the wrappers' truncation
        xyzout.pop()

    metadata = xyzout[0].to_metadata_dict()
    assert metadata.get("annotation") == "Coot output: model.pdb"
    assert metadata.get("subType") == 1


def test_cdatafile_set_dict_preserves_metadata(tmp_path):
    """CDataFile.set(dict) must carry annotation/subType, not just the
    path. This is the core-level fix (cdata_file.py) behind the harvest
    bug: previously the 'fullPath' branch set the path and dropped every
    other key, so CList.set(slice) - and the client's file-parameter
    path - silently lost metadata."""
    from ccp4i2.core.CCP4ModelData import CPdbDataFile

    source = CPdbDataFile()
    source.setFullPath(str(tmp_path / "model.pdb"))
    source.annotation.set("kept through copy")
    source.subType.set(CPdbDataFile.SUBTYPE_MODEL)
    source.contentFlag.set(CPdbDataFile.CONTENT_FLAG_PDB)

    copy = CPdbDataFile()
    copy.set(source.get())  # the get()/set() round trip CList.set() uses

    metadata = copy.to_metadata_dict()
    assert metadata.get("annotation") == "kept through copy"
    assert metadata.get("subType") == 1
    assert metadata.get("contentFlag") == 1


def test_cdatafile_set_pathonly_dict_leaves_other_fields(tmp_path):
    """A dict of only {fullPath} sets just the path - it must not disturb
    metadata already present (the fix applies remaining keys per child,
    never unsetting absent ones)."""
    from ccp4i2.core.CCP4ModelData import CPdbDataFile

    the_file = CPdbDataFile()
    the_file.annotation.set("pre-existing")
    the_file.set({"fullPath": str(tmp_path / "model.pdb")})

    assert str(the_file.annotation) == "pre-existing"
    assert the_file.getFullPath().endswith("model.pdb")
