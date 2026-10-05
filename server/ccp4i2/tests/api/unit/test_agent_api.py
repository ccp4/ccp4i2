"""The endpoints an agent works from (docs/agentic-knowledge.md): a job's
parameters one line each, a task's judgement, and the verdict on a job."""
import pytest
import yaml
from rest_framework.test import APIClient

from ccp4i2.agent import judgement
from ccp4i2.db import models

API = "/api/ccp4i2"


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def project(tmp_path):
    directory = tmp_path / "agentapi"
    (directory / "CCP4_JOBS").mkdir(parents=True)
    return models.Project.objects.create(name="agentapi", directory=str(directory))


def _job(client, project, task="ProvideAsuContents"):
    r = client.post(f"{API}/projects/{project.id}/create_task/",
                    data={"task_name": task}, format="json")
    assert r.status_code == 200, r.content
    return models.Job.objects.get(id=r.json()["data"]["new_job"]["id"])


def test_parameters_are_one_line_each(client, project):
    job = _job(client, project)
    r = client.get(f"{API}/jobs/{job.id}/parameters/")
    assert r.status_code == 200, r.content
    params = r.json()["data"]["parameters"]
    by_path = {p["path"]: p for p in params}
    asu = by_path["inputData.ASU_CONTENT"]
    assert asu["class"] == "CAsuContentSeqList" and asu["set"] is False
    assert all(p["path"].split(".")[0] in ("inputData", "controlParameters", "outputData")
               for p in params)

    only_input = client.get(f"{API}/jobs/{job.id}/parameters/?section=inputData").json()
    assert {p["path"].split(".")[0] for p in only_input["data"]["parameters"]} == {"inputData"}
    found = client.get(f"{API}/jobs/{job.id}/parameters/?query=asu_content").json()
    assert "inputData.ASU_CONTENT" in [p["path"] for p in found["data"]["parameters"]]


def test_a_task_without_judgement_says_so(client, project):
    r = client.get(f"{API}/agent/tasks/freerflag/")  # a task nobody has written up
    assert r.status_code == 200
    assert r.json()["data"]["judgement"] is None
    assert client.get(f"{API}/agent/tasks/no_such_task/").status_code == 404

    job = _job(client, project, task="freerflag")
    verdict = client.get(f"{API}/jobs/{job.id}/judgement/").json()["data"]
    assert verdict["outcome"] is None and "No judgement" in verdict["note"]


def test_a_job_that_has_not_run_is_not_judged(client, project):
    # The rules would read the absent results as a result: a pending
    # ProvideAsuContents was judged "empty", no sequence recorded
    job = _job(client, project)
    verdict = client.get(f"{API}/jobs/{job.id}/judgement/").json()["data"]
    assert verdict["outcome"] is None and "has not run" in verdict["note"]
    assert "results" not in verdict


def test_judgement_is_read_from_the_job(client, project, tmp_path, monkeypatch):
    job = _job(client, project)
    job.status = models.Job.Status.FINISHED
    job.save()
    (job.directory / "program.xml").write_text(
        "<R><Final><RFree>0.24</RFree></Final><Anom>NaN</Anom></R>")
    key, _ = models.JobValueKey.objects.get_or_create(name="RFactor", defaults={"description": "R"})
    models.JobFloatValue.objects.create(job=job, key=key, value=0.2)
    path = tmp_path / "ProvideAsuContents.agent.yaml"
    path.write_text(yaml.safe_dump({
        "task": "ProvideAsuContents", "status": "draft",
        "results": {"RFREE": {"xpath": ".//Final/RFree"},
                    "R": {"file": "kpi", "kpi": "RFactor"},
                    "ANOM": {"xpath": ".//Anom"}},
        "verdict": [{"when": "ANOM > 0", "outcome": "anomalous"},
                    {"when": "RFREE < 0.3 and R < RFREE", "outcome": "good", "basis": "test"},
                    {"when": True, "outcome": "poor"}],
        "next": [{"when": 'outcome == "good"', "task": "servalcat_pipe"}],
    }))
    monkeypatch.setattr(judgement, "judgement_path", lambda task: path)

    described = client.get(f"{API}/agent/tasks/ProvideAsuContents/").json()["data"]
    assert described["judgement"]["results"]["RFREE"]["xpath"] == ".//Final/RFree"

    verdict = client.get(f"{API}/jobs/{job.id}/judgement/").json()["data"]
    # NaN (CTRUNCATE's "no anomalous limit") compares false in the rules and
    # is served as null, which JSON can carry; it was there, so not missing
    assert verdict["results"] == {"RFREE": 0.24, "R": 0.2, "ANOM": None}
    assert verdict["missing"] == []
    assert verdict["outcome"] == "good"
    assert verdict["next"] == [{"when": 'outcome == "good"', "task": "servalcat_pipe"}]
    assert verdict["note"] == judgement.DRAFT_NOTE


def test_task_lookup_says_which_tasks_have_a_judgement(client):
    # An agent choosing between a pipeline and its stages sees which tasks
    # have been written up (Haiku twice started SAD phasing with ShelxCD)
    lookup = client.get(f"{API}/task_lookup/").json()
    assert lookup["shelx"]["hasJudgement"] is True
    assert lookup["ShelxCD"]["hasJudgement"] is True
    assert lookup["editbfac"]["hasJudgement"] is False


def test_a_value_of_the_wrong_shape_is_refused_when_set(client, project):
    # An agent's ASU_CONTENT sequence given as {"text": ...} was stored as the
    # dict's repr and failed only at validation
    import json as _json
    job = _job(client, project)
    r = client.post(f"{API}/jobs/{job.id}/set_parameter/", data=_json.dumps({
        "object_path": "ProvideAsuContents.inputData.ASU_CONTENT",
        "value": [{"name": "A", "sequence": {"text": "MKV"}, "nCopies": 1}]}),
        content_type="application/json")
    assert r.status_code >= 400, r.content
    assert "ASU_CONTENT[0].sequence takes a single value" in r.content.decode()


def test_a_next_step_makes_its_job_and_says_what_it_set(client, project, monkeypatch):
    # The judgement tab's "Rerun" and "Create job": the job is made, its
    # inputs set, and left pending; an input that could not be set is said
    import json as _json
    from pathlib import Path
    from ccp4i2.lib.utils.jobs import apply_next as module
    job = _job(client, project)
    steps = [
        {"task": "ProvideAsuContents", "rerun": True, "advice": "  Two  copies. ",
         "inputs": {"ASU_CONTENT": [{"name": "A", "sequence": "MKVLAAGIV", "nCopies": 2}]}},
        {"task": "ProvideAsuContents", "inputs": {"NOT_A_PARAMETER": "x"}},
        {"advice": "read the report"},
    ]
    monkeypatch.setattr(module, "next_steps", lambda j: steps)

    def apply(index):
        return client.post(f"{API}/jobs/{job.id}/apply_next/", data=_json.dumps({"index": index}),
                           content_type="application/json")

    r = apply(0)
    assert r.status_code == 200, r.content
    out = r.json()["data"]
    assert out["rerun"] is True and out["advice"] == "Two copies."
    assert out["inputs"] == [{"name": "ASU_CONTENT", "value": steps[0]["inputs"]["ASU_CONTENT"],
                              "ok": True}]
    clone = models.Job.objects.get(id=out["job"]["id"])
    assert clone.id != job.id and clone.task_name == "ProvideAsuContents"
    assert "MKVLAAGIV" in (Path(clone.directory) / "input_params.xml").read_text()

    out = apply(1).json()["data"]
    assert out["rerun"] is False and out["inputs"][0]["ok"] is False
    assert out["inputs"][0]["error"]

    r = apply(2)
    assert r.status_code == 400 and "advice only" in r.content.decode()
    assert apply(7).status_code == 400


# --- UniProt into a pending job, through the object_method endpoint --------

def _fake_fetch(accession, residue_range=None, opener=None):
    from ccp4i2.lib.utils.sequences import uniprot
    if accession != "P14635":
        raise uniprot.UniProtError(f"{accession!r} is not a UniProt accession or entry name")
    seq = "MALRVTRNSKINAENKAKINMAGAKRVPTAPAATSKPGLRPRTALGDIGNKVSEQLQAKMPMKKEAKPSATGKVIDKKLPKPLEKVPMLVPVPVSEPVPEPEPEPEPEPVKEEKLSPEPILVDTASPSPMETSGCAPAEEDLCQAFSDVILAVNDVDAEDGADPNLCSEYVKDIYAYLRQLEEEQAVRPKYLLGREVTGNMRAILIDWLVQVQMKFRLLQETMYMTVSIIDRFMQNNCVPKKMLQLVGVTAMFIASKYEEMYPPEIGDFAFVTDNTYTKHQIRQMEMKILRALNFGLGRPLPLHFLRRASKIGEVDVEQHTLAKYLMELTMLDYDMVHFPPSQIAAGAFCLALKILDNGEWTPTLQHYLSYTEESLLPVMQHLAKNVVMVNQGLTKHMTVKNKYATSKHAKISTLPQLNSALVQDLAKAVAKV"
    cut = uniprot.parse_range(residue_range)
    if cut:
        seq = seq[cut[0] - 1:cut[1]]
    return {"accession": "P14635", "entry_name": "CCNB1_HUMAN", "protein_name": "G2/mitotic-specific cyclin-B1",
            "gene": "CCNB1", "organism": "Homo sapiens", "reviewed": True, "taxid": 9606,
            "range": f"{cut[0]}-{cut[1]}" if cut else None, "sequence": seq,
            "fasta": f">sp|P14635|CCNB1_HUMAN cyclin-B1 OS=Homo sapiens\n{seq}\n"}


def _method(client, job, task, name, args):
    import json as _json
    r = client.post(f"{API}/jobs/{job.id}/object_method/", data=_json.dumps(
        {"object_path": task, "method_name": name, "args": args}), content_type="application/json")
    assert r.status_code == 200, r.content
    return r.json()["data"]["result"]


def test_a_sequence_is_fetched_into_a_pending_job_and_saved(client, project, monkeypatch):
    from pathlib import Path
    from ccp4i2.lib.utils.sequences import uniprot
    monkeypatch.setattr(uniprot, "fetch", _fake_fetch)
    job = _job(client, project, task="ProvideSequence")
    out = _method(client, job, "ProvideSequence", "fetchUniProt", ["P14635"])
    assert out["success"] and out["added"]["gene"] == "CCNB1"
    import xml.etree.ElementTree as ET
    saved = ET.parse(Path(job.directory) / "input_params.xml").getroot().findtext(".//SEQUENCETEXT")
    assert saved.startswith(">sp|P14635|CCNB1_HUMAN")  # the record, provenance in its header, saved

    out = _method(client, job, "ProvideSequence", "fetchUniProt", ["cyclin B1"])
    assert out == {"success": False, "error": "'cyclin B1' is not a UniProt accession or entry name"}

    job.status = models.Job.Status.FINISHED
    job.save()
    out = _method(client, job, "ProvideSequence", "fetchUniProt", ["P14635"])
    assert out["success"] is False and "only a pending job" in out["error"]


def test_an_asu_entry_is_filled_from_uniprot_with_its_construct(client, project, monkeypatch):
    from pathlib import Path
    from ccp4i2.lib.utils.sequences import uniprot
    monkeypatch.setattr(uniprot, "fetch", _fake_fetch)
    job = _job(client, project, task="ProvideAsuContents")
    out = _method(client, job, "ProvideAsuContents", "fetchUniProt", ["P14635", "175-432", None, 2])
    assert out["success"] and out["added"]["length"] == 258 and out["added"]["range"] == "175-432"
    saved = (Path(job.directory) / "input_params.xml").read_text()
    assert "YAYLRQLEEE" in saved and "CCNB1" in saved
    out = _method(client, job, "ProvideAsuContents", "fetchUniProt", ["nonsense"])
    assert out["success"] is False
