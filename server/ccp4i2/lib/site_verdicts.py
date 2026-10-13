"""A campaign cell's verdict, when more than one reviewer can score it.

A cell is one dataset at one site. Each reviewer keeps their own
``SiteEvaluation`` there, and what the cell *says* is derived from them
(docs/site-verdict-reconciliation.md):

* no verdicts: nobody has looked (``None``);
* every reviewer gave the same verdict: **agreed**, and that is the verdict;
* reviewers differ, however slightly: **to reconcile**, and the effective
  verdict is ``unclear``, because nobody has settled it. Any difference
  counts, hit beside unclear as much as hit beside empty.

Verdicts with no reviewer were recorded before reviewers were tracked. They
are kept, but they count only at a cell where no named reviewer has scored:
the first attributed verdict is treated as though it were the first verdict
there, so an old unattributed call never puts a cell "to reconcile".

Nothing here writes a verdict (invariant 1 of docs/pandda-campaign-design.md,
section 9.1).
"""

from typing import Dict, Iterable, Optional, Tuple

AGREED = "agreed"
TO_RECONCILE = "to_reconcile"
UNCLEAR = "unclear"


def cell_verdict(rows: Iterable, reviewer: Optional[str] = None) -> Optional[dict]:
    """What one cell says, from its ``SiteEvaluation`` rows (or anything with
    ``evaluator`` and ``verdict``), or None when it has none.

    Returns ``state``, ``verdict`` (the effective verdict), ``verdicts`` (every
    reviewer's, unattributed ones included, for the record) and ``mine`` (the
    given reviewer's verdict, or None).
    """
    rows = list(rows)
    if not rows:
        return None
    counted = [row for row in rows if row.evaluator] or rows
    distinct = {row.verdict for row in counted}
    if len(distinct) == 1:
        state, verdict = AGREED, distinct.pop()
    else:
        state, verdict = TO_RECONCILE, UNCLEAR
    mine = None
    if reviewer is not None:
        mine = next((row.verdict for row in rows if row.evaluator == reviewer), None)
    return {
        "state": state,
        "verdict": verdict,
        "verdicts": [
            {"evaluator": row.evaluator, "verdict": row.verdict}
            for row in sorted(rows, key=lambda row: row.evaluator)
        ],
        "mine": mine,
    }


def cell_verdicts(evaluations: Iterable,
                  reviewer: Optional[str] = None) -> Dict[Tuple[int, int], dict]:
    """``cell_verdict`` for every cell the evaluations touch, keyed by
    ``(project_id, site_id)``."""
    by_cell: Dict[Tuple[int, int], list] = {}
    for evaluation in evaluations:
        by_cell.setdefault((evaluation.project_id, evaluation.site_id), []).append(evaluation)
    return {cell: cell_verdict(rows, reviewer) for cell, rows in by_cell.items()}
