"""Pluggable molecular-replacement engine for ``molrep_map``.

Place a model into a *phased* cryo-EM map. The map is P1, which has no origin for
an unphased translation function -- so the engine must search against the map's
own phased density, never a bare TF. molrep does this with ``-f <map>`` (it scrapes
the "phased translation function" block), and is the default, proven engine.

Phaser is declared as an alternative -- motivated by the PHIL-driven Phaser stack
-- but guarded off until its phased-P1-map placement mode is confirmed; see
``docs/molrep-map-design.md`` sections 5 and 6. Wiring it means feeding it the box
as phased structure factors and driving its PHIL tree; the crystallography (which
Phaser mode uses the phases, not an unphased P1 TF) is the gate, not the plumbing.
"""

from __future__ import annotations

import logging
import os
import subprocess
from dataclasses import dataclass
from typing import Optional

logger = logging.getLogger(__name__)

MOLREP = "molrep"
PHASER = "phaser"
ENGINES = (MOLREP, PHASER)


@dataclass
class PlacementResult:
    """Outcome of one placement run (one hand, one engine)."""
    model_path: Optional[str]      # molrep.pdb, or None if it produced nothing
    doc_path: Optional[str]        # molrep.doc, for report scraping
    log_path: str
    score: Optional[float]         # best-effort scalar for hand recommendation
    timed_out: bool = False

    @property
    def placed(self) -> bool:
        return self.model_path is not None and os.path.exists(self.model_path)


def place(engine: str, map_path, model_path, work_dir, *,
          nmon: int = 1, np_peaks: int = 5, b_add: float = 0.0,
          resmax: Optional[float] = None,
          time_limit: Optional[float] = None) -> PlacementResult:
    """Place ``model_path`` into ``map_path`` using ``engine``.

    ``resmax`` (Angstrom) caps the search resolution -- the primary speed lever
    for triage; ``time_limit`` (seconds) bounds wall-clock. Returns a
    :class:`PlacementResult`; a run that places nothing is reported, not raised.
    """
    if engine == MOLREP:
        return _place_molrep(map_path, model_path, work_dir,
                             nmon=nmon, np_peaks=np_peaks, b_add=b_add,
                             resmax=resmax, time_limit=time_limit)
    if engine == PHASER:
        raise NotImplementedError(
            "The phaser engine is declared but not yet wired: placing a model "
            "against a phased P1 map needs Phaser's phased-map mode confirmed "
            "(see docs/molrep-map-design.md section 6). Use engine='molrep'.")
    raise ValueError(f"Unknown MR engine {engine!r}; expected one of {ENGINES}")


def _place_molrep(map_path, model_path, work_dir, *, nmon, np_peaks, b_add,
                  resmax, time_limit) -> PlacementResult:
    os.makedirs(work_dir, exist_ok=True)
    com_lines = []
    if nmon:
        com_lines.append(f"nmon {nmon}")
    if np_peaks:
        com_lines.append(f"np {np_peaks}")
    if b_add:
        com_lines.append(f"badd {b_add}")
    if resmax:
        com_lines.append(f"resmax {resmax}")
    com = ("\n".join(com_lines) + "\n") if com_lines else "\n"

    log_path = os.path.join(work_dir, "log.txt")
    cmd = [MOLREP, "-f", str(map_path), "-m", str(model_path), "-i"]
    timed_out = False
    with open(log_path, "w") as logf:
        try:
            subprocess.run(cmd, cwd=work_dir, input=com, text=True,
                           stdout=logf, stderr=subprocess.STDOUT,
                           timeout=time_limit)
        except subprocess.TimeoutExpired:
            timed_out = True
            logger.warning("molrep exceeded the %s s time limit in %s",
                           time_limit, work_dir)

    doc = os.path.join(work_dir, "molrep.doc")
    pdb = os.path.join(work_dir, "molrep.pdb")
    doc_path = doc if os.path.exists(doc) else None
    return PlacementResult(
        model_path=pdb if os.path.exists(pdb) else None,
        doc_path=doc_path,
        log_path=log_path,
        score=_scrape_molrep_score(doc_path) if doc_path else None,
        timed_out=timed_out,
    )


def _scrape_molrep_score(doc_path: str) -> Optional[float]:
    """Best-effort scalar score from ``molrep.doc`` for the hand recommendation.

    Reads the (phased) translation-function peak table and returns the largest
    correlation-like value (the last numeric column of any peak row). Advisory
    only -- the task emits both hands regardless -- so a parse miss returns None
    rather than failing.
    """
    try:
        with open(doc_path) as fh:
            lines = fh.read().split("\n")
    except OSError:
        return None

    best = None
    titles = []
    in_tf = False
    for line in lines:
        stripped = line.strip()
        if stripped in ("--- Translation function ---",
                        "--- phased translation function ---"):
            in_tf = True
            continue
        if not in_tf:
            continue
        if stripped.startswith("RF "):
            titles = (line.replace("(", " ").replace(")", "")
                      .replace("/", "_").split())
            continue
        words = (line.replace("(", " ").replace(")", "")
                 .replace("-", " -").split())
        if titles and len(words) == len(titles):
            try:
                int(words[0]); int(words[1])
                value = float(words[-1])
            except (ValueError, IndexError):
                continue
            best = value if best is None else max(best, value)
    return best
