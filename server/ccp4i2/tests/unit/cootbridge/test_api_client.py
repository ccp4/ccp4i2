"""CCP4-free tests for the cootbridge shared data layer.

The module under test must run on stock Python (it also targets Coot's
embedded interpreters), so nothing here may import Django, gemmi, or
any GUI toolkit.
"""

import json
import os

import pytest

from ccp4i2.cootbridge import api_client


# ---------------------------------------------------------------------------
# BridgeConfig: the environment contract
# ---------------------------------------------------------------------------


def test_config_explicit_api_url_wins_and_trailing_slash_stripped():
    config = api_client.BridgeConfig(
        {"CCP4I2_API_URL": "http://127.0.0.1:43434/", "UVICORN_PORT": "999"})
    assert config.api_url == "http://127.0.0.1:43434"


def test_config_falls_back_to_uvicorn_port_then_default():
    assert (api_client.BridgeConfig({"UVICORN_PORT": "5151"}).api_url
            == "http://127.0.0.1:5151")
    assert api_client.BridgeConfig({}).api_url == "http://127.0.0.1:8000"


def test_config_token_falls_back_to_local_session_token():
    config = api_client.BridgeConfig(
        {"CCP4I2_LOCAL_SESSION_TOKEN": "sekrit"})
    assert config.token == "sekrit"
    config = api_client.BridgeConfig(
        {"CCP4I2_ACCESS_TOKEN": "a", "CCP4I2_LOCAL_SESSION_TOKEN": "b"})
    assert config.token == "a"


def test_config_derives_drop_dir_from_job_directory():
    config = api_client.BridgeConfig(
        {"CCP4I2_JOB_DIRECTORY": os.path.join("proj", "CCP4_JOBS", "job_9")})
    assert config.drop_dir == os.path.join(
        "proj", "CCP4_JOBS", "job_9", "COOT_FILE_DROP")


# ---------------------------------------------------------------------------
# params.xml parsing and the load plan
# ---------------------------------------------------------------------------

PARAMS_XML = """<?xml version='1.0' encoding='ASCII'?>
<ccp4:ccp4i2 xmlns:ccp4="http://www.ccp4.ac.uk/ccp4ns">
  <ccp4i2_header><function>PARAMS</function></ccp4i2_header>
  <ccp4i2_body>
    <inputData>
      <XYZIN_LIST>
        <CPdbDataFile>
          <dbFileId>abc</dbFileId>
          <baseName>model1.pdb</baseName>
          <project>PROJ-UUID</project>
          <relPath>CCP4_JOBS/job_3</relPath>
          <annotation>Refined model</annotation>
        </CPdbDataFile>
        <CPdbDataFile>
          <baseName>model2.cif</baseName>
          <project>PROJ-UUID</project>
          <relPath>CCP4_JOBS/job_4</relPath>
        </CPdbDataFile>
      </XYZIN_LIST>
      <FPHIIN>
        <baseName>FPHIOUT.mtz</baseName>
        <project>PROJ-UUID</project>
        <subType>1</subType>
        <relPath>CCP4_JOBS/job_3</relPath>
        <annotation>2Fo-Fc map</annotation>
      </FPHIIN>
      <DICT>
        <baseName>LIG.cif</baseName>
        <project>PROJ-UUID</project>
        <relPath>CCP4_IMPORTED_FILES</relPath>
      </DICT>
      <UNRELATED><baseName>ignored.txt</baseName></UNRELATED>
    </inputData>
  </ccp4i2_body>
</ccp4:ccp4i2>
"""


def test_parse_input_params_extracts_lists_singles_and_order():
    entries = api_client.parse_input_params(PARAMS_XML)
    kinds = [entry["kind"] for entry in entries]
    names = [entry["base_name"] for entry in entries]
    # Dictionaries first (models may need them), then declared order.
    assert kinds == ["dictionary", "coordinates", "coordinates", "map_2fofc"]
    assert names == ["LIG.cif", "model1.pdb", "model2.cif", "FPHIOUT.mtz"]
    assert entries[1]["annotation"] == "Refined model"
    assert entries[1]["db_file_id"] == "abc"
    assert entries[2]["db_file_id"] == ""


def test_parse_input_params_tolerates_missing_input_data():
    assert api_client.parse_input_params(
        "<ccp4i2><ccp4i2_body/></ccp4i2>") == []


def _project_with_files(tmp_path):
    project = tmp_path / "proj"
    for rel, name in [
        ("CCP4_JOBS/job_3", "model1.pdb"),
        ("CCP4_JOBS/job_4", "model2.cif"),
        ("CCP4_JOBS/job_3", "FPHIOUT.mtz"),
        ("CCP4_IMPORTED_FILES", "LIG.cif"),
    ]:
        target = project.joinpath(*rel.split("/"))
        target.mkdir(parents=True, exist_ok=True)
        (target / name).write_text("x")
    return project


def test_resolve_entry_path_local_project(tmp_path):
    project = _project_with_files(tmp_path)
    config = api_client.BridgeConfig({
        "CCP4I2_PROJECT_UUID": "proj-uuid",
        "CCP4I2_PROJECT_DIRECTORY": str(project),
    })
    entry = {"base_name": "model1.pdb", "rel_path": "CCP4_JOBS/job_3",
             "project": "PROJ-UUID", "db_file_id": ""}
    resolved = api_client.resolve_entry_path(entry, config)
    assert resolved == str(project / "CCP4_JOBS" / "job_3" / "model1.pdb")


def test_resolve_entry_path_returns_none_when_unresolvable(tmp_path):
    config = api_client.BridgeConfig(
        {"CCP4I2_PROJECT_DIRECTORY": str(tmp_path)})
    entry = {"base_name": "nope.pdb", "rel_path": "CCP4_JOBS/job_1",
             "project": "", "db_file_id": ""}
    assert api_client.resolve_entry_path(entry, config) is None


def test_load_plan_offline_from_job_directory(tmp_path):
    project = _project_with_files(tmp_path)
    job_dir = project / "CCP4_JOBS" / "job_9"
    job_dir.mkdir(parents=True)
    (job_dir / "input_params.xml").write_text(PARAMS_XML)
    config = api_client.BridgeConfig({
        "CCP4I2_PROJECT_UUID": "proj-uuid",
        "CCP4I2_PROJECT_DIRECTORY": str(project),
        "CCP4I2_JOB_DIRECTORY": str(job_dir),
    })
    plan = api_client.load_plan(config, client=None)
    assert [item["kind"] for item in plan] == [
        "dictionary", "coordinates", "coordinates", "map_2fofc"]
    assert plan[1]["label"] == "Refined model"
    assert plan[2]["label"] == "model2.cif"  # falls back to base name
    assert all(os.path.exists(item["path"]) for item in plan)


# ---------------------------------------------------------------------------
# Browse model
# ---------------------------------------------------------------------------


def test_classify_file_types_and_map_subtypes():
    assert api_client.classify_file(
        {"type": "chemical/x-pdb"}) == "coordinates"
    assert api_client.classify_file(
        {"type": "application/CCP4-mtz-map", "sub_type": 2}) == "map_fofc"
    assert api_client.classify_file(
        {"type": "application/CCP4-mtz-map", "sub_type": 3}) == "map_anom"
    assert api_client.classify_file(
        {"type": "application/CCP4-mtz-map", "sub_type": None}) == "map_2fofc"
    assert api_client.classify_file(
        {"type": "application/refmac-dictionary"}) == "dictionary"
    assert api_client.classify_file(
        {"type": "application/CCP4-unmerged-mtz"}) is None


JOB_TREE = [
    {
        "id": 11, "uuid": "u11", "number": "2", "title": "Refine",
        "status": 6,
        "files": [
            {"id": 101, "uuid": "f101", "name": "XYZOUT.pdb",
             "type": "chemical/x-pdb", "sub_type": 1,
             "annotation": "Refined model"},
            {"id": 102, "uuid": "f102", "name": "FPHIOUT.mtz",
             "type": "application/CCP4-mtz-map", "sub_type": 1,
             "annotation": ""},
            {"id": 103, "uuid": "f103", "name": "log.txt",
             "type": "text/plain", "sub_type": None, "annotation": ""},
        ],
        "children": [
            {"id": 12, "uuid": "u12", "number": "2.1", "title": "sub",
             "status": 6, "files": [], "children": []},
        ],
    },
    {
        "id": 10, "uuid": "u10", "number": "1", "title": "Import",
        "status": 6,
        "files": [{"id": 100, "uuid": "f100", "name": "LIG.cif",
                   "type": "application/refmac-dictionary", "sub_type": None,
                   "annotation": "Ligand dictionary"}],
        "children": [],
    },
]


def test_order_for_load_puts_dictionaries_first():
    files = [
        {"kind": "coordinates", "label": "model"},
        {"kind": "map_2fofc", "label": "2fofc"},
        {"kind": "dictionary", "label": "LIG"},
        {"kind": "coordinates", "label": "model2"},
    ]
    ordered = api_client.order_for_load(files)
    assert [f["kind"] for f in ordered] == [
        "dictionary", "coordinates", "coordinates", "map_2fofc"]
    # stable within a priority band: the two coordinates keep their order
    assert [f["label"] for f in ordered if f["kind"] == "coordinates"] == [
        "model", "model2"]


def test_filter_rows_case_insensitive_substring():
    rows = [("2: Refine model", 2), ("1: Import data", 1),
            ("3: refmac run", 3)]
    key = lambda r: r[0]
    # blank query -> everything, order preserved
    assert api_client.filter_rows(rows, "", key) == rows
    assert api_client.filter_rows(rows, "   ", key) == rows
    # case-insensitive substring
    assert [r[1] for r in api_client.filter_rows(rows, "ref", key)] == [2, 3]
    assert [r[1] for r in api_client.filter_rows(rows, "IMPORT", key)] == [1]
    assert api_client.filter_rows(rows, "zzz", key) == []


def test_display_label_prefixes_only_cross_project():
    # same project -> unchanged
    assert api_client.display_label("model", "projA", "projA") == "model"
    # other project -> prefixed, so the two stay distinct in Coot
    assert api_client.display_label("model", "projB", "projA") == \
        "projB: model"
    # unknown source project -> unchanged (no misleading prefix)
    assert api_client.display_label("model", None, "projA") == "model"


def test_browse_model_orders_dictionaries_before_coordinates():
    tree = [{
        "id": 1, "uuid": "u1", "number": "1", "title": "build", "status": 6,
        "files": [
            {"id": 1, "uuid": "f1", "name": "XYZOUT.pdb",
             "type": "chemical/x-pdb", "sub_type": 1, "annotation": ""},
            {"id": 2, "uuid": "f2", "name": "LIG.cif",
             "type": "application/refmac-dictionary", "sub_type": None,
             "annotation": "Ligand"},
        ],
        "children": [],
    }]
    files = api_client.browse_model(tree)[0]["files"]
    # dictionary must lead so coordinate parsing sees ligand geometry
    assert files[0]["kind"] == "dictionary"
    assert files[1]["kind"] == "coordinates"


def test_browse_model_flattens_and_filters():
    rows = api_client.browse_model(JOB_TREE)
    assert [row["label"] for row in rows] == ["2: Refine", "1: Import"]
    first = rows[0]
    assert [f["kind"] for f in first["files"]] == ["coordinates", "map_2fofc"]
    assert first["files"][0]["label"] == "Refined model"
    assert first["files"][1]["label"] == "FPHIOUT.mtz"  # no annotation
    # The empty sub-job is skipped by default but kept when asked for.
    assert len(api_client.browse_model(JOB_TREE, include_empty=True)) == 3
    depths = [row["depth"] for row in
              api_client.browse_model(JOB_TREE, include_empty=True)]
    assert depths == [0, 1, 0]


# ---------------------------------------------------------------------------
# Drop-dir save/harvest contract
# ---------------------------------------------------------------------------


def test_output_numbering_and_harvest_order(tmp_path):
    drop = tmp_path / "COOT_FILE_DROP"
    assert api_client.next_output_number(str(drop)) == 1
    path1 = api_client.output_path(str(drop), 1, "pdb")
    open(path1, "w").close()
    (drop / "output3.cif").write_text("x")
    (drop / "outputNaN.pdb").write_text("x")  # non-contract name, ignored
    assert api_client.next_output_number(str(drop)) == 4
    numbers = [n for n, _ in api_client.harvestable_outputs(str(drop))]
    assert numbers == [1, 3]


# ---------------------------------------------------------------------------
# HTTP client: envelope handling and auth header
# ---------------------------------------------------------------------------


class _FakeResponse(object):
    def __init__(self, payload):
        self._payload = payload

    def read(self):
        return self._payload

    def close(self):
        pass


def test_client_sends_bearer_token_and_unwraps_envelope():
    seen = {}

    def opener(request, timeout=None):
        seen["url"] = request.get_full_url()
        seen["auth"] = request.get_header("Authorization")
        return _FakeResponse(json.dumps(
            {"success": True, "data": {"xml": "<x/>"}}).encode())

    config = api_client.BridgeConfig({
        "CCP4I2_API_URL": "http://127.0.0.1:4000",
        "CCP4I2_ACCESS_TOKEN": "tok",
    })
    client = api_client.CootBridgeClient(config, opener=opener)
    assert client.params_xml(7) == "<x/>"
    assert seen["url"] == "http://127.0.0.1:4000/api/ccp4i2/jobs/7/params_xml/"
    assert seen["auth"] == "Bearer tok"


def test_client_passes_bare_json_through_and_raises_on_failure():
    def opener_ok(request, timeout=None):
        return _FakeResponse(json.dumps([{"id": 1}]).encode())

    config = api_client.BridgeConfig({})
    client = api_client.CootBridgeClient(config, opener=opener_ok)
    assert client.projects() == [{"id": 1}]

    def opener_fail(request, timeout=None):
        return _FakeResponse(json.dumps(
            {"success": False, "error": "nope"}).encode())

    failing = api_client.CootBridgeClient(config, opener=opener_fail)
    with pytest.raises(api_client.BridgeError):
        failing.get_json("projects/")


def test_job_tree_unwraps_the_job_tree_key():
    """The real endpoint answers {"job_tree": [...], "total_jobs": ...}."""
    def opener(request, timeout=None):
        return _FakeResponse(json.dumps({
            "job_tree": JOB_TREE, "total_jobs": 3, "total_files": 4,
        }).encode())

    client = api_client.CootBridgeClient(
        api_client.BridgeConfig({}), opener=opener)
    tree = client.job_tree(2)
    assert isinstance(tree, list) and tree[0]["number"] == "2"
    rows = api_client.browse_model(tree)
    assert [row["label"] for row in rows] == ["2: Refine", "1: Import"]


def test_project_by_uuid_normalises_hyphens():
    def opener(request, timeout=None):
        return _FakeResponse(json.dumps([
            {"id": 1, "uuid": "AABB-CCDD"},
            {"id": 2, "uuid": "11223344"},
        ]).encode())

    client = api_client.CootBridgeClient(
        api_client.BridgeConfig({}), opener=opener)
    assert client.project_by_uuid("aabbccdd")["id"] == 1
