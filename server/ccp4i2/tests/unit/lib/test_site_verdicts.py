"""What a campaign cell says, from its reviewers' verdicts
(lib/site_verdicts.py, docs/site-verdict-reconciliation.md). Pure: no
database, no CCP4."""

from types import SimpleNamespace as Row

from ccp4i2.lib.site_verdicts import cell_verdict, cell_verdicts


def rows(*pairs):
    return [Row(evaluator=who, verdict=what, project_id=1, site_id=1) for who, what in pairs]


def test_nobody_has_looked():
    assert cell_verdict([]) is None


def test_one_reviewer_is_agreed():
    cell = cell_verdict(rows(("alice", "hit")), reviewer="alice")
    assert (cell["state"], cell["verdict"], cell["mine"]) == ("agreed", "hit", "hit")


def test_any_difference_is_to_reconcile_and_reads_as_unclear():
    for other in ("empty", "unclear"):
        cell = cell_verdict(rows(("alice", "hit"), ("bob", other)))
        assert (cell["state"], cell["verdict"]) == ("to_reconcile", "unclear")


def test_unattributed_verdicts_count_only_alone():
    assert cell_verdict(rows(("", "empty")))["verdict"] == "empty"
    cell = cell_verdict(rows(("", "empty"), ("alice", "hit")), reviewer="bob")
    assert (cell["state"], cell["verdict"], cell["mine"]) == ("agreed", "hit", None)
    assert cell["verdicts"] == [{"evaluator": "", "verdict": "empty"},
                                {"evaluator": "alice", "verdict": "hit"}]


def test_cells_are_keyed_by_project_and_site():
    evaluations = rows(("alice", "hit")) + [Row(evaluator="alice", verdict="empty", project_id=2, site_id=1)]
    assert {key: c["verdict"] for key, c in cell_verdicts(evaluations).items()} == {
        (1, 1): "hit", (2, 1): "empty"}
