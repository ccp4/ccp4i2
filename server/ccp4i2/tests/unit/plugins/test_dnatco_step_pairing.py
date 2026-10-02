"""Comparing two models, the DNATCO report pairs each step with its
counterpart. Step names begin with the entry id, which differs between the
files ("1hr2" and "custom" for 1hr2 and its PDB-REDO re-refinement), so every
step was listed twice and the per-step plot put model 2 beside model 1."""
from ccp4i2.wrappers.dnatco.script.dnatco_report import _merged_steps


def _model(entry, n):
    steps = [{"name": f"{entry}_A_U_{200 + i}_C_{201 + i}", "chain": "A", "step": i}
             for i in range(n)]
    return {"ntc": {"steps": steps}}


def test_steps_pair_across_entries():
    merged = _merged_steps([_model("1hr2", 3), _model("custom", 3)])
    assert len(merged) == 3
    assert all(None not in entry["per_model"] for entry in merged)
