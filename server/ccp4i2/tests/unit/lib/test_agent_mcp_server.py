"""The MCP facade (ccp4i2/agent/mcp_server.py): a thin client of the REST
API. Here the API is faked: what is tested is what the tools ask of it and
what they hand the agent."""
import asyncio

import pytest

pytest.importorskip("mcp", reason="needs the mcp package (pip install ccp4i2[agent])")

from mcp.server.mcpserver.exceptions import ToolError  # noqa: E402

from ccp4i2.agent import mcp_server  # noqa: E402

EXPECTED_TOOLS = {
    "list_projects", "create_project", "list_tasks", "describe_task", "project_jobs",
    "create_job", "clone_job", "job_parameters", "set_parameter", "set_file",
    "upload_file", "validate", "run_job", "job_status", "wait_for_job", "judge_job",
    "what_next",
}


@pytest.fixture
def api(monkeypatch):
    calls = []
    answers = {}

    def fake(method, path, body=None, query=None, content_type="application/json"):
        calls.append((method, path, body, query))
        answer = answers.get((method, path))
        if isinstance(answer, Exception):
            raise answer
        return answer

    monkeypatch.setattr(mcp_server, "_request", fake)
    return calls, answers


def test_the_tools_offered():
    tools = asyncio.run(mcp_server.server.list_tools())
    assert {t.name for t in tools} == EXPECTED_TOOLS
    # Nothing that deletes or rewrites a job's status.
    assert not any("delete" in t.name or "status" in t.name and t.name != "job_status" for t in tools)


def test_set_file_resolves_the_reference_and_drops_the_full_path(api):
    calls, answers = api
    answers[("GET", "jobs/7")] = {"id": 7, "project": 3}
    answers[("GET", "projects/3/resolve_fileuse")] = {
        "project": "p", "baseName": "PHASER.1.pdb", "dbFileId": "f", "relPath": "CCP4_JOBS/job_5",
        "annotation": "Positioned coordinates", "fullPath": "/somewhere/PHASER.1.pdb"}
    out = mcp_server.set_file(7, "inputData.XYZIN", "[5].XYZOUT[0]")
    assert out == {"path": "inputData.XYZIN", "file": "PHASER.1.pdb",
                   "annotation": "Positioned coordinates"}
    method, path, body, _ = calls[-1]
    assert (method, path) == ("POST", "jobs/7/set_parameter")
    assert body["object_path"] == "inputData.XYZIN"
    assert "fullPath" not in body["value"] and body["value"]["dbFileId"] == "f"
    assert calls[1][3] == {"fileuse": "[5].XYZOUT[0]"}


def test_validate_splits_errors_from_warnings(api):
    _, answers = api
    answers[("GET", "jobs/7/validation")] = {"xml": (
        "<errorReportList><errorReport><severity>4</severity><description>no data</description>"
        "</errorReport><errorReport><severity>2</severity><description>no free R</description>"
        "</errorReport></errorReportList>")}
    answers[("GET", "jobs/7/run_time_validation")] = {"xml": "<errorReportList/>"}
    out = mcp_server.validate(7)
    assert [e["description"] for e in out["errors"]] == ["no data"]
    assert [w["description"] for w in out["warnings"]] == ["no free R"]


def test_a_refusal_reaches_the_agent_as_its_reason(api):
    _, answers = api
    answers[("GET", "jobs/7")] = {"id": 7, "project": 3}
    answers[("POST", "jobs/7/set_parameter")] = mcp_server.ApiError("Use clone API first")
    # The SDK hands a ToolError's message to the model (an is_error result
    # over the protocol); any other exception would reach it as a bare
    # "Error executing tool".
    with pytest.raises(ToolError, match="Use clone API first"):
        asyncio.run(mcp_server.server.call_tool(
            "set_parameter", {"job_id": 7, "path": "inputData.NCOPIES", "value": 2}))


def test_list_tasks_leaves_out_superseded_and_interactive(api):
    _, answers = api
    answers[("GET", "task_lookup")] = {
        "phaser_simple_phil": {"TASKTITLE": "Phaser basic", "DESCRIPTION": "MR"},
        "phaser_simple": {"TASKTITLE": "Phaser old", "supersededBy": "phaser_simple_phil"},
        "moorhen_rebuild": {"TASKTITLE": "Moorhen", "interactive": True},
    }
    assert [t["task"] for t in mcp_server.list_tasks("phaser")["tasks"]] == ["phaser_simple_phil"]


def test_wait_for_job_returns_when_the_job_ends(api, monkeypatch):
    _, answers = api
    answers[("GET", "jobs/7")] = {"id": 7, "number": "3", "task_name": "x", "status": 6}
    monkeypatch.setattr("time.sleep", lambda s: None)
    out = mcp_server.wait_for_job(7, max_seconds=5)
    assert out["status"] == "Finished" and out["waited_out"] is False


def test_a_list_is_set_whole_with_its_files_resolved(api):
    calls, answers = api
    answers[("GET", "jobs/7")] = {"id": 7, "project": 3}
    answers[("GET", "projects/3/resolve_fileuse")] = {
        "project": "p", "baseName": "beta.pdb", "dbFileId": "f1", "relPath": "CCP4_JOBS/job_3",
        "fullPath": "/x/beta.pdb"}
    answers[("GET", "jobs/7/parameters")] = {"parameters": [
        {"path": "inputData.ENSEMBLES", "set": True, "value": [{"label": "beta"}]}]}
    out = mcp_server.set_parameter(7, "inputData.ENSEMBLES", [
        {"label": "beta", "number": 1, "use": True,
         "pdbItemList": [{"structure": {"file": "[3].XYZOUT"}, "identity_to_target": 1.0}]}])
    posted = next(c for c in calls if c[1] == "jobs/7/set_parameter")[2]
    structure = posted["value"][0]["pdbItemList"][0]["structure"]
    assert structure == {"project": "p", "baseName": "beta.pdb", "dbFileId": "f1",
                         "relPath": "CCP4_JOBS/job_3"}
    assert posted["value"][0]["label"] == "beta"
    assert out == {"path": "inputData.ENSEMBLES", "now": [{"label": "beta"}], "set": True}


def test_a_value_echoed_from_job_parameters_round_trips(api):
    calls, answers = api
    answers[("GET", "jobs/7")] = {"id": 7, "project": 3}
    answers[("GET", "projects/3/resolve_fileuse")] = {"project": "p", "baseName": "a.mtz",
                                                      "dbFileId": "abc"}
    answers[("GET", "jobs/7/parameters")] = {"parameters": []}
    mcp_server.set_parameter(7, "inputData.F_SIGF", {"file": "a.mtz", "annotation": "x", "fileId": "abc"})
    resolve = next(c for c in calls if c[1] == "projects/3/resolve_fileuse")
    assert resolve[3] == {"fileuse": "abc"}


def test_severity_as_the_server_writes_it(api):
    # The validation endpoint writes <severity>ERROR</severity>: read as a
    # number it was 0, and an agent was told there were no errors.
    _, answers = api
    answers[("GET", "jobs/7/validation")] = {"xml": (
        "<errorReportList><errorReport><code>113</code><description>add a sequence"
        "</description><severity>ERROR</severity></errorReport><errorReport>"
        "<severity>WARNING</severity><description>no free set</description>"
        "</errorReport></errorReportList>")}
    answers[("GET", "jobs/7/run_time_validation")] = {"xml": "<errorReportList/>"}
    out = mcp_server.validate(7)
    assert [e["code"] for e in out["errors"]] == ["113"]
    assert [w["description"] for w in out["warnings"]] == ["no free set"]
