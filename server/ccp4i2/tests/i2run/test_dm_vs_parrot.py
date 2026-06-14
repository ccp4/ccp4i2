"""Compare per-domain NCS averaging (dm_multidomain) with whole-monomer NCS
averaging (parrot), recovering from a realistic poor-MR phase set on AHIR.

Benchmark design (chosen so there is a *known truth*):
  * data        = model |Fc| (self-consistent, so the model PHIC IS the truth)
  * start phases = a 2.0 A coordinate-shaken model -> resolution-dependent error
                   (the sigmaA/Luzzati shape, not a flat per-shell error),
                   weighted by a sigmaA FOM.
  * metric      = mean |phase error| to truth PER RESOLUTION SHELL, reported
                  SEPARATELY for acentric and centric reflections. Centric
                  phases are restricted to two values 180 deg apart, so they
                  must not be averaged in with acentrics; density modification
                  also treats them specially. Global map CC is dominated by
                  strong low-res terms and hides where DM actually acts.

Both tasks run through their real i2run wrappers on identical input. Density
modification only has headroom when phases are poor (this poor-MR regime); on a
good MR solution neither helps -- the textbook behaviour.

Split into two test functions (one job each) because the in-process Django test
DB does not support two job-creating runs in a single test; parrot writes its
per-shell errors to a JSON file the dm run reads to print the comparison.
"""
import json
import os

import numpy as np
import gemmi

from .utils import demoData, i2run

_HANDOFF = "/tmp/dm_spike/compare_results.json"

# Fixed resolution-shell edges (A), equal-ish reflection counts for this cell.
_EDGES = [99.0, 4.60, 3.63, 3.16, 2.87, 2.66, 2.40]


def _read_phi(mtz_path):
    """({hkl -> phase deg}, cell, spacegroup) from the first phase column."""
    m = gemmi.read_mtz_file(str(mtz_path))
    phi_col = next(c for c in m.columns if c.type == 'P')
    hkl = m.make_miller_array()
    phi = phi_col.array
    d = {tuple(int(x) for x in hkl[i]): phi[i]
         for i in range(len(phi)) if not np.isnan(phi[i])}
    return d, m.cell, m.spacegroup


def _truth_phi():
    return _read_phi(demoData("ahir", "ahir_calc.mtz"))[0]


def _shell_errors(phi, truth, cell, sg):
    """Per-shell mean |phase error| to truth (deg), split acentric / centric."""
    go = sg.operations()
    nb = len(_EDGES) - 1
    acen = [[] for _ in range(nb)]
    cen = [[] for _ in range(nb)]
    for h, p in phi.items():
        if h not in truth:
            continue
        d = cell.calculate_d(h)
        dphi = abs((p - truth[h] + 180) % 360 - 180)
        for b in range(nb):
            if _EDGES[b] >= d > _EDGES[b + 1]:
                (cen if go.is_reflection_centric(h) else acen)[b].append(dphi)
                break
    mean = lambda L: [round(float(np.mean(x)), 1) if x else float('nan') for x in L]
    return {"acentric": mean(acen), "centric": mean(cen)}


def _run(task_args):
    truth = _truth_phi()
    start_phi, cell, sg = _read_phi(demoData("ahir", "ahir_phases_mr.mtz"))
    start = _shell_errors(start_phi, truth, cell, sg)
    with i2run(task_args) as job:
        out_phi, _, _ = _read_phi(job / "FPHIOUT.mtz")
    out = _shell_errors(out_phi, truth, cell, sg)
    return start, out


def test_parrot_recovery():
    """Whole-monomer NCS averaging (parrot, NCS auto-detected from the model)."""
    start, parrot = _run([
        "parrot",
        "--F_SIGF", demoData("ahir", "ahir_fmodel.mtz"),
        "--ABCD", demoData("ahir", "ahir_phases_mr.mtz"),
        "--ASUIN", demoData("ahir", "ahir.asu.xml"),
        "--XYZIN_MODE", "mr",
        "--XYZIN_MR", demoData("ahir", "ahir_model.cif"),
        "--CYCLES", "10"])
    os.makedirs(os.path.dirname(_HANDOFF), exist_ok=True)
    with open(_HANDOFF, "w", encoding="utf-8") as fh:
        json.dump({"start": start, "parrot": parrot}, fh)
    print("\nparrot acentric:", parrot["acentric"], " centric:", parrot["centric"])
    assert len(parrot["acentric"]) == len(_EDGES) - 1


def test_dm_multidomain_recovery():
    """Per-domain NCS averaging; prints the centric/acentric comparison table."""
    start, dm = _run([
        "dm_multidomain",
        "--F_SIGF", demoData("ahir", "ahir_fmodel.mtz"),
        "--ABCD", demoData("ahir", "ahir_phases_mr.mtz"),
        "--XYZIN", demoData("ahir", "ahir_model.cif"),
        "--DOMAINS", "chainId=A", "firstRes=340", "lastRes=485", "mode=average",
        "--DOMAINS", "chainId=A", "firstRes=140", "lastRes=339", "mode=average",
        "--DOMAINS", "chainId=A", "firstRes=13", "lastRes=139", "mode=average",
        "--NCYCLES", "10"])

    prior = {}
    if os.path.exists(_HANDOFF):
        with open(_HANDOFF, encoding="utf-8") as fh:
            prior = json.load(fh)
    parrot = prior.get("parrot")

    shells = "  ".join(f"{_EDGES[i]:.2f}" for i in range(len(_EDGES)))
    print("\n=== mean |phase error| to truth, per resolution shell (deg) ===")
    print(f"  shell edges (A): {shells}")
    for kind in ("acentric", "centric"):
        print(f"  --- {kind} ---")
        print(f"  start  : {start[kind]}")
        if parrot:
            print(f"  parrot : {parrot[kind]}  (whole-monomer)")
        print(f"  dm     : {dm[kind]}  (per-domain)")

    # the demonstration: per-domain beats whole-monomer on acentric mid-res
    # shells (where averaging has headroom), by a clear margin.
    if parrot:
        assert dm["acentric"][1] < parrot["acentric"][1] - 2.0
        assert dm["acentric"][2] < parrot["acentric"][2] - 2.0
