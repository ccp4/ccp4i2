"""The judgement files an agent reads: the ``when`` language, reading a job's
results, and the verdict (docs/agentic-knowledge.md)."""
from pathlib import Path

import pytest

from ccp4i2.agent import condition, judgement
from ccp4i2.core.tasks import TASKS


@pytest.mark.parametrize("text, values, expected", [
    ("TFZ >= 8", {"TFZ": 8.0}, True),
    ("TFZ >= 8", {"TFZ": 7.9}, False),
    ("TFZ >= 8 and LLG > 60", {"TFZ": 9, "LLG": 50}, False),
    ("TFZ >= 8 or LLG > 60", {"TFZ": 9, "LLG": 50}, True),
    ("not (RFREE > 0.35)", {"RFREE": 0.3}, True),
    ('outcome == "solved"', {"outcome": "solved"}, True),
    ("true", {}, True),
    ("false or RFREE < 1e-1", {"RFREE": 0.05}, True),
    ("-1 < X", {"X": 0}, True),
    ("A - B >= 0.02", {"A": 0.30, "B": 0.25}, True),
    ("A - B >= 0.02", {"A": 0.30, "B": 0.29}, False),
    ("A-B>0", {"A": 2, "B": 1}, True),
    ("A / B > 1.1 and -A < 0", {"A": 0.5, "B": 0.4}, True),
    ("1 + 2 * 3 == 7", {}, True),
    ("(1 + 2) * 3 == 9", {}, True),
    ("NCS == null", {}, True),
    ("NCS == null", {"NCS": "mr"}, False),
    ("NCS != null and NCS == \"mr\"", {"NCS": "mr"}, True),
    ("NCS != null", {"NCS": None}, False),
    ("not (NCS == null) or X > 1", {}, False),
])
def test_conditions(text, values, expected):
    assert condition.holds(text, values) is expected


def test_a_missing_result_never_passes():
    # A comparison with a number that could not be read is false either way,
    # and so is its negation: a missing TFZ is never a solution.
    assert not condition.holds("TFZ >= 8", {})
    assert not condition.holds("TFZ < 8", {"TFZ": None})
    assert not condition.holds("not (TFZ < 8)", {})
    assert not condition.holds('TFZ > "x"', {"TFZ": 3})
    assert not condition.holds("A - B > 0", {"A": 1})
    assert not condition.holds("A / B > 0", {"A": 1, "B": 0})


@pytest.mark.parametrize("text", [
    "TFZ >=", "(TFZ > 8", "TFZ > 8)", "TFZ > 8 and", "__import__('os')",
    "TFZ > 8; x", "and TFZ", "A - ", "A * * B",
])
def test_not_conditions(text):
    with pytest.raises(condition.ConditionError):
        condition.parse(text)


def test_names():
    assert condition.names(condition.parse("A > 1 and (B < 2 or not C == 3)")) == {"A", "B", "C"}


PROGRAM_XML = """<PHASER><Rama>4.01%</Rama><Solution><TFZ>12.4</TFZ><LLG>310</LLG>
<Space attr="P 21 21 21"/></Solution><Solution><TFZ>5.0</TFZ></Solution></PHASER>"""

JUDGEMENT = {
    "task": "example",
    "status": "draft",
    "results": {
        "TFZ": {"xpath": ".//Solution/TFZ"},
        "LLG": {"xpath": ".//Solution/LLG", "type": "int"},
        "SG": {"xpath": ".//Solution/Space", "attribute": "attr", "type": "str"},
        "RFREE": {"file": "kpi", "kpi": "RFree"},
        "ABSENT": {"xpath": ".//Nowhere"},
        "RAMA": {"xpath": ".//Rama"},
        "NOTE": {"xpath": ".//Note", "type": "str", "optional": True},
    },
    "verdict": [
        {"when": "TFZ >= 8 and LLG > 60", "outcome": "solved", "basis": "doc"},
        {"when": "TFZ >= 6", "outcome": "ambiguous", "basis": "doc"},
        {"when": True, "outcome": "failed"},
    ],
    "next": [
        {"when": 'outcome == "solved"', "task": "servalcat_pipe"},
        {"when": 'outcome != "solved"', "advice": "try another model"},
    ],
}


@pytest.fixture
def job_dir(tmp_path):
    (tmp_path / "program.xml").write_text(PROGRAM_XML)
    return tmp_path


def test_results_are_read_from_the_job(job_dir):
    values = judgement.read_results(JUDGEMENT, job_dir, kpis={"RFree": 0.27})
    assert values == {"TFZ": 12.4, "LLG": 310, "SG": "P 21 21 21",
                      "RFREE": 0.27, "ABSENT": None, "RAMA": 4.01, "NOTE": None}


def test_verdict_first_that_holds(job_dir):
    verdict = judgement.judge("example", job_dir, judgement=JUDGEMENT)
    assert verdict["outcome"] == "solved"
    assert verdict["missing"] == ["ABSENT", "RFREE"]
    assert [n.get("task") for n in verdict["next"]] == ["servalcat_pipe"]
    assert verdict["note"] == judgement.DRAFT_NOTE


def test_a_job_with_no_program_xml_fails_through(tmp_path):
    verdict = judgement.judge("example", tmp_path, judgement=JUDGEMENT)
    assert verdict["outcome"] == "failed"
    assert verdict["next"][0]["advice"] == "try another model"


def test_reviewed_judgement_carries_no_draft_note(job_dir):
    verdict = judgement.judge("example", job_dir, judgement=dict(JUDGEMENT, status="reviewed"))
    assert "note" not in verdict


def test_problems_are_named():
    bad = {
        "results": {"TFZ": {}, "SG": {"xpath": ".//SG", "type": "text"},
                    "BAD": {"xpath": ".//A[@b!='c'"}},
        "verdict": [{"when": "TFZ > 8", "outcome": "solved"},
                    {"when": "LLG > 1", "outcome": "x", "basis": "b"},
                    {"when": "TFZ >", "basis": "b"},
                    {"when": True, "outcome": "failed"}],
    }
    assert judgement.problems(bad) == [
        "result TFZ: no xpath",
        "result SG: type 'text' is not one of ['float', 'int', 'str', 'string']",
        "result BAD: xpath \".//A[@b!='c'\" is not one ElementTree reads "
        "('NoneType' object is not callable); keep to tags, /, //, [n], [last()], "
        "[@a='v'], [tag='v']",
        "verdict[0]: a threshold with no basis",
        "verdict[1]: unknown result ['LLG']",
        "verdict[2]: 'TFZ >' ends too soon",
    ]
    assert judgement.problems(JUDGEMENT) == []


def _shipped():
    for task in TASKS:
        path = judgement.judgement_path(task)
        if path is not None and path.is_file():
            yield pytest.param(task, path, id=task)


@pytest.mark.parametrize("task, path", list(_shipped()) or [pytest.param(None, None, marks=pytest.mark.skip("none written yet"))])
def test_shipped_judgement_files_are_well_formed(task, path):
    data = judgement.load(path=path)
    assert data["task"] == task
    assert data.get("status") in ("draft", "reviewed")
    assert judgement.problems(data) == []


def test_an_xpath_list_takes_the_first_that_reads(tmp_path):
    (tmp_path / "program.xml").write_text(
        "<R><Before><RFree>0.282</RFree></Before><After><RFree></RFree></After></R>")
    spec = {"xpath": [".//After/RFree", ".//Before/RFree"]}
    assert judgement.read_result(spec, tmp_path) == 0.282   # After is empty
    (tmp_path / "program.xml").write_text(
        "<R><Before><RFree>0.282</RFree></Before><After><RFree>0.275</RFree></After></R>")
    assert judgement.read_result(spec, tmp_path) == 0.275
    assert judgement.problems({"results": {"X": {"xpath": [".//A", ".//B[@c!='d'"]}}})


def test_a_finished_job_is_judged_once_until_its_judgement_changes(tmp_path, monkeypatch):
    # Judged once and kept as judgement.json, with the judgement's version:
    # a finished job's files do not change, but the draft judgements do
    import json
    from ccp4i2.agent import judgement as J
    rules = tmp_path / "task.agent.yaml"
    rules.write_text(
        "task: t\nstatus: draft\nresults:\n  R:\n    file: program.xml\n    xpath: .//R\n"
        "    type: float\nverdict:\n  - when: R < 0.3\n    outcome: good\n"
        "  - when: true\n    outcome: bad\n")
    monkeypatch.setattr(J, "judgement_path", lambda name: rules)
    job = tmp_path / "job"
    job.mkdir()
    (job / "program.xml").write_text("<x><R>0.2</R></x>")

    first = J.judge_finished("t", job)
    assert first["outcome"] == "good" and first["judgement_version"]
    kept = json.loads((job / J.CACHE_NAME).read_text())
    assert kept["judgement_version"] == first["judgement_version"]

    (job / "program.xml").write_text("<x><R>0.5</R></x>")  # files of a finished job do not change...
    assert J.judge_finished("t", job)["outcome"] == "good"  # ...so the kept verdict stands

    rules.write_text(rules.read_text().replace("R < 0.3", "R < 0.6"))  # the judgement is revised
    again = J.judge_finished("t", job)
    assert again["judgement_version"] != first["judgement_version"]
    assert again["outcome"] == "good" and again["results"]["R"] == 0.5  # judged afresh


# --- what the app draws: clauses, tone, meanings, reruns, pinned references --

def test_clauses_give_each_comparison_its_value_and_threshold():
    out = condition.clauses("TFZ >= 8 and RFREE < 0.55", {"TFZ": 12.7, "RFREE": 0.52})
    assert [(c["name"], c["op"], c["threshold"], c["value"], c["holds"]) for c in out] == [
        ("TFZ", ">=", 8.0, 12.7, True), ("RFREE", "<", 0.55, 0.52, True)]


def test_a_clause_on_a_missing_result_is_unknown_and_arithmetic_has_no_gauge():
    out = condition.clauses("RFREE_START - RFREE >= 0.02 or X == null", {"RFREE": 0.3})
    assert out[0]["holds"] is None and "threshold" not in out[0]
    assert out[1] == {"text": "X == null", "holds": True, "values": {"X": None}}


def test_a_number_on_the_left_is_drawn_as_a_threshold_on_the_right():
    (clause,) = condition.clauses("0.4 > RFREE", {"RFREE": 0.3})
    assert (clause["name"], clause["op"], clause["threshold"]) == ("RFREE", "<", 0.4)


def _rules(tmp_path, extra=""):
    path = tmp_path / "t.agent.yaml"
    path.write_text(
        "task: t\nstatus: draft\nresults:\n  R:\n    file: program.xml\n    xpath: .//R\n"
        "    type: float\n    meaning: >\n      R-free of the job.\n"
        "verdict:\n  - when: R < 0.3\n    outcome: refined\n    basis: b\n"
        "  - when: true\n    outcome: failed\n"
        "next:\n  - when: outcome == \"refined\"\n    rerun: true\n"
        "    inputs:\n      ADD_WATERS: \"True\"\n"
        "  - when: outcome == \"refined\"\n    task: other\n"
        "    inputs:\n      XYZIN: \"t[-1].XYZOUT\"\n      HKLIN: \"other[-1].X\"\n" + extra)
    return path


def test_a_verdict_carries_tone_meanings_clauses_and_reruns_name_the_task(tmp_path, monkeypatch):
    monkeypatch.setattr(judgement, "judgement_path", lambda name: _rules(tmp_path))
    job = tmp_path / "job"
    job.mkdir()
    (job / "program.xml").write_text("<x><R>0.2</R></x>")
    v = judgement.judge("t", job)
    assert v["tone"] == "good" and v["meanings"] == {"R": "R-free of the job."}
    assert v["clauses"][0]["holds"] is True
    assert v["next"][0]["task"] == "t" and v["next"][0]["rerun"] is True


def test_only_this_tasks_latest_is_pinned_to_this_job():
    steps = [{"task": "other", "inputs": {"XYZIN": "t[-1].XYZOUT", "HKLIN": "other[-1].X",
                                          "N": 3}}]
    pinned = judgement.pin_references(steps, "t", "5")
    assert pinned[0]["inputs"] == {"XYZIN": "[5].XYZOUT", "HKLIN": "other[-1].X", "N": 3}
    assert steps[0]["inputs"]["XYZIN"] == "t[-1].XYZOUT"  # the original is untouched


def test_a_rerun_must_be_of_this_task_and_change_something(tmp_path):
    import yaml
    rules = yaml.safe_load(_rules(tmp_path).read_text())
    assert judgement.problems(rules) == []
    rules["next"][0]["task"] = "else"
    rules["next"].append({"when": "true", "rerun": True})
    found = judgement.problems(rules)
    assert any("a rerun is of this task" in p for p in found)
    assert any("a rerun with nothing changed" in p for p in found)
