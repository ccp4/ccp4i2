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


PROGRAM_XML = """<PHASER><Solution><TFZ>12.4</TFZ><LLG>310</LLG>
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
                      "RFREE": 0.27, "ABSENT": None}


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
        "results": {"TFZ": {}, "SG": {"xpath": ".//SG", "type": "text"}},
        "verdict": [{"when": "TFZ > 8", "outcome": "solved"},
                    {"when": "LLG > 1", "outcome": "x", "basis": "b"},
                    {"when": "TFZ >", "basis": "b"},
                    {"when": True, "outcome": "failed"}],
    }
    assert judgement.problems(bad) == [
        "result TFZ: no xpath",
        "result SG: type 'text' is not one of ['float', 'int', 'str', 'string']",
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
