# Copyright (C) 2026 Newcastle University
"""End-to-end tests for the simple single-datatype Import<Type> tasks."""
from os import environ
from pathlib import Path

import gemmi

from .utils import demoData, i2run


def _find_output(job: Path, suffix: str, exclude_prefixes=("HKLIN", "XYZIN", "DICTIN", "MAPIN", "SEQIN", "ASUIN")):
    """Return the single import output file in the job dir with the given suffix."""
    candidates = [
        f for f in job.iterdir()
        if f.suffix.lower() == suffix and not f.name.startswith(exclude_prefixes)
    ]
    assert len(candidates) == 1, f"expected one *{suffix} output, found {[f.name for f in candidates]}"
    return candidates[0]


def _check_mtz(path: Path, expected_cols):
    mtz = gemmi.read_mtz_file(str(path))
    labels = set(c.label for c in mtz.columns)
    assert labels == {"H", "K", "L"} | set(expected_cols), labels
    assert mtz.nreflections > 0


# ---------------------------------------------------------------------------
# Pass-through imports (no columns)
# ---------------------------------------------------------------------------


def test_import_coordinate():
    args = ["ImportCoordinate", "--XYZIN", demoData("gamma", "gamma_model.pdb")]
    with i2run(args) as job:
        out = _find_output(job, ".pdb")
        structure = gemmi.read_structure(str(out))
        assert len(structure) > 0


def test_import_dictionary():
    cif = str(Path(environ["CLIBD_MON"], "a", "A1LU6.cif"))
    args = ["ImportDictionary", "--DICTIN", cif]
    with i2run(args) as job:
        out = _find_output(job, ".cif")
        assert out.exists() and out.stat().st_size > 0


# ---------------------------------------------------------------------------
# Mini-MTZ imports — auto-detect (single unambiguous group)
# ---------------------------------------------------------------------------


def test_import_freer_auto():
    args = ["ImportFreeR", "--HKLIN", demoData("gamma", "freeR.mtz")]
    with i2run(args) as job:
        _check_mtz(_find_output(job, ".mtz"), ["FREER"])


def test_import_phases_auto():
    args = ["ImportPhases", "--HKLIN", demoData("gamma", "initial_phases.mtz")]
    with i2run(args) as job:
        _check_mtz(_find_output(job, ".mtz"), ["HLA", "HLB", "HLC", "HLD"])


def test_import_obs_auto():
    args = ["ImportObs", "--HKLIN", demoData("gamma", "merged_intensities_native.mtz")]
    with i2run(args) as job:
        _check_mtz(_find_output(job, ".mtz"), ["Iplus", "SIGIplus", "Iminus", "SIGIminus"])


# ---------------------------------------------------------------------------
# Mini-MTZ imports — ambiguous file, explicit column selection (hybrid)
# ---------------------------------------------------------------------------


def test_import_mapcoeffs_explicit_columns():
    # mergedForMakingCif.mtz has two FPHI groups (FPHIOUT_* and DIFFPHIOUT_*)
    args = ["ImportMapCoeffs", "--HKLIN", demoData("gamma", "mergedForMakingCif.mtz")]
    args += ["--COLUMNS", "FPHIOUT_F,FPHIOUT_PHI"]
    with i2run(args) as job:
        _check_mtz(_find_output(job, ".mtz"), ["F", "PHI"])


def test_import_mapcoeffs_ambiguous_fails():
    # No COLUMNS given and >1 candidate group -> must fail at validation
    args = ["ImportMapCoeffs", "--HKLIN", demoData("gamma", "mergedForMakingCif.mtz")]
    failed = False
    try:
        with i2run(args):
            pass
    except Exception:
        failed = True
    assert failed, "ambiguous map-coefficient import should have been rejected"
