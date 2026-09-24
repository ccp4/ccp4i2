"""The PanDDA 2 invocation contract, as code: argv, environment, progress
symbol, failure catalogue and sizing.

Design note: docs/pandda-campaign-design.md, sections 4.1-4.5 and 6.
Source of truth: ``CCP4I2_PANDDA_INVOCATION_CONTRACT.md`` in the Materia
repository (adopted verbatim, decision 4). Nothing here runs anything; the
plugin composes these into a job.
"""
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

#: What this module implements. Recorded in every run's provenance (4.4).
CONTRACT_VERSION = "CCP4I2_PANDDA_INVOCATION_CONTRACT.md (materia, 2026-06-06)"

#: The name PanDDA gives its output tree, which fan-out and Reinspect both
#: read. Part of the contract; never renamed (4.3).
OUT_DIR_NAME = "pandda2_out"

PROGRAM = "pandda2.analyse"

#: PanDDA takes the crystal number from the last run of digits in a directory
#: name and its default range is 0-99999 on the CCP4 build; anything above
#: is silently dropped. Always explicit.
DATASET_RANGE = "0-999999999"

#: Never matches anything we stage. PanDDA's default matches final.pdb and
#: takes the protein model for a ligand.
LIGAND_PDB_REGEX = "ligand.pdb"

MODEL_REGEX = "final.pdb"
REFLECTIONS_REGEX = "final.mtz"
DICT_REGEX = "dict.cif"

#: Below this PanDDA refuses to characterise a ground state at all
#: (``--min_characterisation_datasets`` default).
MIN_DATASETS = 25


def build_argv(data_dirs, out_dir, local_cpus: int) -> List[str]:
    """The contract's argv, argument for argument, both defensive literals
    included. ``data_dirs`` is the staged ``datasets/`` directory."""
    return [
        "--data_dirs", str(data_dirs),
        "--out_dir", str(out_dir),
        "--local_cpus", str(int(local_cpus)),
        "--pdb_regex", MODEL_REGEX,
        "--mtz_regex", REFLECTIONS_REGEX,
        "--ligand_cif_regex", DICT_REGEX,
        "--ligand_pdb_regex", LIGAND_PDB_REGEX,
        "--dataset_range", DATASET_RANGE,
    ]


def build_env(base_env: Dict[str, str], scratch_dir) -> Dict[str, str]:
    """The child's environment: the caller's, with ``RAY_TMPDIR`` set and
    every ``PANDDA_*`` switch removed.

    Removed, not left alone: the experimental switches of section 6.5 are
    presence-checked, so a stray ``PANDDA_...=0`` in a developer's shell would
    *enable* an unvalidated code path. The run's environment is what this
    function builds, not what the shell happened to hold.
    """
    env = {k: v for k, v in base_env.items() if not k.startswith("PANDDA_")}
    env["RAY_TMPDIR"] = str(scratch_dir)
    return env


# --- progress (4.5) ---------------------------------------------------------

PROGRESS_RE = re.compile(r"^PANDDA_PROGRESS: dataset (\d+)/(\d+)$", re.MULTILINE)


def parse_progress(text: str) -> Optional[Tuple[int, int]]:
    """The last ``(done, total)`` PanDDA reported, or None: absence means
    progress unknown (an older build), never failure."""
    matches = PROGRESS_RE.findall(text or "")
    if not matches:
        return None
    done, total = matches[-1]
    return int(done), int(total)


# --- failure catalogue (4.5), seeded from the contract ---------------------

#: (name, pattern over stderr+log, ERROR_CODES code, what to do about it)
FAILURE_CATALOGUE = [
    ("oom", re.compile(r"MemoryError|Killed: 9|SIGKILL|Out of memory", re.I), 210,
     "PanDDA ran out of memory. Fewer datasets per run is the one lever that "
     "moves this (section 6); or run it on a larger machine and fan out from the tree"),
    ("free_r_label", re.compile(r"No RFree Flag found", re.I), 211,
     "PanDDA did not recognise the free-R column; staging should have relabelled it "
     "- report this with the manifest"),
    ("ligand_block", re.compile(r"KeyError: [\"']block 'comp_", re.I), 212,
     "This PanDDA build reads the ligand block by name; a newer build reads it by content"),
    ("dataset_range_zeroed", re.compile(r"0/\d+ datasets passed range filter", re.I), 213,
     "No dataset passed PanDDA's range filter: the staged names were not clean xtal-NNNN"),
    ("ray_scratch_full", re.compile(r"No space left on device", re.I), 214,
     "The scratch disk filled. Point SCRATCH_DIR at a disk with tens of GB free"),
    ("ccp4_missing", re.compile(r"(refmac5|gemmi): command not found", re.I), 215,
     "PanDDA's environment could not see CCP4; the launcher did not source it"),
]

UNCLASSIFIED = ("unclassified_crash", None, 216,
                "PanDDA exited with an error this task does not recognise; read the log")


def classify_failure(text: str):
    """``(name, code, prompt)`` for the first catalogue entry whose pattern
    matches ``text`` (stderr followed by the log), else the unclassified
    entry. Never the exit-code value: that carries only success/failure."""
    for name, pattern, code, prompt in FAILURE_CATALOGUE:
        if pattern.search(text or ""):
            return name, code, prompt
    name, _pattern, code, prompt = UNCLASSIFIED
    return name, code, prompt


# --- sizing (6.1, 6.4) ------------------------------------------------------

#: Rough per-map cost by cell class, MB. Map size scales with cell volume:
#: ~4 MB for a small bromodomain, ~22-25 MB for a long-axis cell (6.1).
MAP_MB = {"small": 4.0, "medium": 12.0, "large": 25.0}
#: Maps PanDDA holds per comparator dataset in a shell (xmap, mean, sigma,
#: z-map ...): the multiplier that turns the 6.1 formula into the observed
#: 30-40 GB single shell.
MAPS_PER_DATASET = 4
BASELINE_GIB = 6.0
#: PanDDA characterises a shell from at most this many comparators.
MAX_COMPARATORS = 60


def cell_volume_class(cell) -> str:
    """``small``/``medium``/``large`` from a unit cell ``(a, b, c, alpha,
    beta, gamma)``. Volume moves memory more than dataset count does, and a
    long axis moves it most (6.1)."""
    a, b, c = float(cell[0]), float(cell[1]), float(cell[2])
    longest = max(a, b, c)
    volume = a * b * c
    if longest >= 150.0 or volume >= 8.0e5:
        return "large"
    if longest >= 100.0 or volume >= 5.0e5:
        return "medium"
    return "small"


_CLASS_ORDER = {"small": 0, "medium": 1, "large": 2}


def sizing_hint(n_datasets: int, cells) -> Dict[str, object]:
    """The run contract's ``sizing_hint``: computed once, used locally as a
    warning and, if the run is delegated, as the payload field (6.4)."""
    worst = "small"
    for cell in cells:
        klass = cell_volume_class(cell)
        if _CLASS_ORDER[klass] > _CLASS_ORDER[worst]:
            worst = klass
    return {"datasets": int(n_datasets), "cell_volume_class": worst}


def estimate_peak_gib(n_datasets: int, volume_class: str, local_cpus: int) -> float:
    """A coarse peak-RAM estimate from the 6.1 formula: comparators x map
    size x workers, plus a baseline. Meant to say "this will not fit", not
    to be right to the gigabyte."""
    comparators = min(int(n_datasets), MAX_COMPARATORS)
    maps_mb = MAP_MB.get(volume_class, MAP_MB["large"])
    return BASELINE_GIB + comparators * maps_mb * MAPS_PER_DATASET * max(1, int(local_cpus)) / 1024.0


# --- the executable (4.6) ---------------------------------------------------

def probe_executable(path) -> Dict[str, object]:
    """What can be said about the PanDDA behind ``path`` without running it.

    Both candidate builds report version 0.0.1 and carry no commit, so the
    resolved path plus a capability probe is the only usable provenance
    there is (4.6). The probe looks for the ``PANDDA_PROGRESS`` symbol in the
    installed source, which is the one behaviour this task adapts to.
    """
    path = Path(path)
    probe = {"path": str(path), "launcher": None, "version": None,
             "progress_signal": None, "site_packages": None}
    try:
        head = path.read_text(errors="replace")[:2000]
    except OSError:
        return probe
    probe["launcher"] = head.strip().splitlines()[-1][:200] if head.strip() else None
    env_root = None
    match = re.search(r"-r\s+(\S+)\s+-n\s+(\S+)", head)
    if match:
        root, name = match.group(1), match.group(2)
        root = os.path.expandvars(root)
        env_root = Path(root) / "envs" / name
    if env_root and env_root.is_dir():
        site = next(iter(sorted(env_root.glob("lib/python*/site-packages"))), None)
        if site is not None:
            probe["site_packages"] = str(site)
            dist = next(iter(sorted(site.glob("pandda_gemmi-*.dist-info"))), None)
            if dist is not None:
                probe["version"] = dist.name[len("pandda_gemmi-"):-len(".dist-info")]
            source = site / "pandda_gemmi" / "pandda" / "pandda.py"
            if source.is_file():
                try:
                    probe["progress_signal"] = "PANDDA_PROGRESS" in source.read_text(errors="replace")
                except OSError:
                    pass
    return probe
