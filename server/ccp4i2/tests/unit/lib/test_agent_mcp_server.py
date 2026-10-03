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
    "what_next", "job_errors", "file_summary",
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
    # Project 3 holds job number "5" (server id 7) and its step "5.1" (id 8)
    answers[("GET", "projects/3/job_tree")] = {"job_tree": [
        {"id": 7, "number": "5", "status": 5, "task_name": "x",
         "children": [{"id": 8, "number": "5.1", "status": 5, "task_name": "y", "children": []}]}]}
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
    out = mcp_server.set_file(3, "5", "inputData.XYZIN", "[5].XYZOUT[0]")
    assert out == {"path": "inputData.XYZIN", "file": "PHASER.1.pdb",
                   "annotation": "Positioned coordinates"}
    method, path, body, _ = calls[-1]
    assert (method, path) == ("POST", "jobs/7/set_parameter")
    assert body["object_path"] == "inputData.XYZIN"
    assert "fullPath" not in body["value"] and body["value"]["dbFileId"] == "f"
    resolve = next(c for c in calls if c[1] == "projects/3/resolve_fileuse")
    assert resolve[3] == {"fileuse": "[5].XYZOUT[0]"}


def test_validate_splits_errors_from_warnings(api):
    _, answers = api
    answers[("GET", "jobs/7/validation")] = {"xml": (
        "<errorReportList><errorReport><severity>4</severity><description>no data</description>"
        "</errorReport><errorReport><severity>2</severity><description>no free R</description>"
        "</errorReport></errorReportList>")}
    answers[("GET", "jobs/7/run_time_validation")] = {"xml": "<errorReportList/>"}
    out = mcp_server.validate(3, "5")
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
            "set_parameter", {"project_id": 3, "job": "5", "path": "inputData.NCOPIES", "value": 2}))


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
    answers[("GET", "jobs/7")] = {"id": 7, "number": "5", "task_name": "x", "status": 6, "project": 3}
    monkeypatch.setattr("time.sleep", lambda s: None)
    out = mcp_server.wait_for_job(3, "5", max_seconds=5)
    assert out["status"] == "Finished" and out["waited_out"] is False


def test_a_list_is_set_whole_with_its_files_resolved(api):
    calls, answers = api
    answers[("GET", "jobs/7")] = {"id": 7, "project": 3}
    answers[("GET", "projects/3/resolve_fileuse")] = {
        "project": "p", "baseName": "beta.pdb", "dbFileId": "f1", "relPath": "CCP4_JOBS/job_3",
        "fullPath": "/x/beta.pdb"}
    answers[("GET", "jobs/7/parameters")] = {"parameters": [
        {"path": "inputData.ENSEMBLES", "set": True, "value": [{"label": "beta"}]}]}
    out = mcp_server.set_parameter(3, "5", "inputData.ENSEMBLES", [
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
    mcp_server.set_parameter(3, "5", "inputData.F_SIGF", {"file": "a.mtz", "annotation": "x", "fileId": "abc"})
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
    out = mcp_server.validate(3, "5")
    assert [e["code"] for e in out["errors"]] == ["113"]
    assert [w["description"] for w in out["warnings"]] == ["no free set"]


def test_a_value_inside_a_list_item_is_echoed(api):
    _, answers = api
    answers[("GET", "jobs/7/parameters")] = {"parameters": [{"path": "inputData.ENSEMBLES", "value": [
        {"label": "beta", "pdbItemList": [{"structure": {"file": "beta.pdb"}}]}]}]}
    assert mcp_server._value_at(7, "inputData.ENSEMBLES[0].pdbItemList[0].structure") == (
        {"file": "beta.pdb"}, True)
    assert mcp_server._value_at(7, "inputData.ENSEMBLES[1].label") == (None, False)


def test_a_job_is_named_by_its_number_in_its_project(api):
    _, answers = api
    assert mcp_server._jid(3, "5") == 7 and mcp_server._jid(3, "5.1") == 8
    # A number that is not this project's is refused, never another project's job
    with pytest.raises(mcp_server.ApiError, match="No job 9 in project 3"):
        mcp_server._jid(3, "9")


def test_job_errors_reports_the_job_and_its_failed_steps(api, monkeypatch):
    _, answers = api
    xml = ("<errorReportList><errorReport><code>201</code><description>No output"
           "</description><severityName>WARNING</severityName></errorReport></errorReportList>")
    answers[("GET", "jobs/7/diagnostic_xml")] = {"xml": xml}
    answers[("GET", "jobs/8/diagnostic_xml")] = {"xml": "<errorReportList/>"}
    logs = {"projects/3/files_by_path/CCP4_JOBS/job_5/log.txt": "a\n\nlast line\n",
            "projects/3/files_by_path/CCP4_JOBS/job_5/job_1/log.txt": "step log"}
    monkeypatch.setattr(mcp_server, "_text", lambda path: logs.get(path, ""))
    out = mcp_server.job_errors(3, "5", log_lines=1)
    assert out["status"] == "Failed" and out["errors"][0]["code"] == "201"
    assert out["log_tail"] == "last line"
    assert out["failed_steps"][0]["job"] == "5.1" and out["failed_steps"][0]["log_tail"] == "step log"


def test_a_job_with_validation_errors_is_not_started(api):
    calls, answers = api
    answers[("GET", "jobs/7/validation")] = {"xml": (
        "<errorReportList><errorReport><severity>ERROR</severity><objectPath>x.inputData.HKLIN"
        "</objectPath><description>Required value not set</description></errorReport>"
        "</errorReportList>")}
    answers[("GET", "jobs/7/run_time_validation")] = {"xml": "<errorReportList/>"}
    with pytest.raises(mcp_server.ApiError, match="Required value not set"):
        mcp_server.run_job(3, "5")
    assert not any(c[1] == "jobs/7/run" for c in calls)


def test_a_job_number_may_be_given_as_a_number(api):
    import asyncio as _asyncio
    _, answers = api
    answers[("GET", "jobs/7")] = {"id": 7, "number": "5", "task_name": "x", "status": 6, "project": 3}
    # An agent sent "job": 1 and the argument check refused it; now it is the
    # same job as "5" (the call would raise ToolError on a rejected argument)
    _asyncio.run(mcp_server.server.call_tool("job_status", {"project_id": 3, "job": 5}))
    assert mcp_server.job_status(3, 5)["job"] == "5"


def test_file_summary_digests_the_file_a_reference_names(api):
    # What an agent could not see without it: whether the heavy atoms
    # survived into a rebuilt model, which residues were left unbuilt
    calls, answers = api
    answers[("GET", "projects/3/resolve_fileuse")] = {
        "baseName": "modelcraft.cif", "dbFileId": "abc", "relPath": "CCP4_JOBS/job_8",
        "annotation": "ModelCraft model", "fullPath": "/somewhere/modelcraft.cif"}
    answers[("GET", "files_by_uuid/abc/digest")] = {"composition": {"elements": ["S"]}}
    out = mcp_server.file_summary(3, "[8].XYZOUT")
    assert out == {"file": {"name": "modelcraft.cif", "annotation": "ModelCraft model",
                            "job_directory": "CCP4_JOBS/job_8"},
                   "summary": {"composition": {"elements": ["S"]}}}
    assert ("GET", "projects/3/resolve_fileuse", None, {"fileuse": "[8].XYZOUT"}) in calls
    assert "/somewhere" not in str(out)


def test_file_summary_of_nothing_says_so(api):
    _, answers = api
    answers[("GET", "projects/3/resolve_fileuse")] = {}
    with pytest.raises(ToolError, match="no file found"):
        mcp_server.file_summary(3, "[9].XYZOUT")
