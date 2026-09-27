"""Read one dataset's share of a PanDDA 2 output tree.

Design note: docs/pandda-campaign-design.md, sections 4.3, 7.1-7.4, 8.4.

The tree is the invocation contract's, as PanDDA writes it
(``pandda_gemmi/constants``)::

    <out>/
    ├── analyses/pandda_analyse_events.csv        # run-level; may be absent
    └── processed_datasets/<dtag>/
        ├── <dtag>-pandda-input.pdb               # symlink to the staged apo model
        ├── <dtag>-z_map.native.ccp4
        ├── <dtag>-ground-state-average-map.native.ccp4
        ├── events.yaml                           # {} when there are no events
        ├── <dtag>-event_<idx>_1-BDC_<token>_map.native.ccp4
        ├── <dtag>_event_<idx>_best_autobuild.pdb
        ├── autobuild/…                           # every build tried
        └── modelled_structures/<dtag>-pandda-model.pdb   # merged; may be absent

Two facts that must not be rediscovered (section 7.3): the apo
``-pandda-input.pdb`` is the model of record and every pose is a candidate
merged onto it; and ``Optimal Contour`` is in absolute map units. One more
from the filenames: the ``BDC_<token>`` in an event map's name is *not* the
BDC in ``events.yaml`` (it is ``1 - BDC``, rounded), so maps are matched by
event index and never by that token.

This module reads and reports; it moves no bytes and touches no database.
"""
import csv
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(f"ccp4i2:{__name__}")

ANALYSES_DIR = "analyses"
EVENTS_TABLE = "pandda_analyse_events.csv"
PROCESSED_DIR = "processed_datasets"
EVENTS_YAML = "events.yaml"
MODELLED_DIR = "modelled_structures"
LIGAND_FILES_DIR = "ligand_files"
#: The staged tree's own names, as pandda_campaign.pandda_staging writes them.
#: Imported rather than guessed would be better, but the receipt must import on
#: a machine with no campaign wrapper present, so they are mirrored with this
#: note: pandda_staging.MODEL_NAME, .DICT_NAME and .LIGAND_DIR_NAME.
STAGED_MODEL_NAME = "final.pdb"
STAGED_DICT_NAME = "dict.cif"
STAGED_LIGAND_DIR = "compound"
#: Where that tree sits relative to the run's output directory: both are
#: children of the campaign job's own directory (pandda_campaign writes
#: staging/ beside pandda2_out/).
STAGING_DIR_NAME = "staging"
STAGED_DATASETS_DIR = "datasets"
#: What PanDDA names every residue it builds, whatever the dictionary said.
PANDDA_RESIDUE_NAME = "LIG"


class DatasetNotFound(Exception):
    """The tree has no processed directory for this dtag: nothing to receipt."""


@dataclass
class BuildRecord:
    path: Optional[Path]           # the pose file PanDDA selected, if it exists
    build_score: Optional[float]
    rscc: Optional[float]
    optimal_contour: Optional[float]
    raw: dict = field(default_factory=dict)


@dataclass
class EventRecord:
    idx: int
    bdc: Optional[float]
    score: Optional[float]
    centroid: Optional[tuple]      # (x, y, z) orthogonal Angstroms
    build: Optional[BuildRecord]
    event_map: Optional[Path]      # None when PanDDA recorded it but the file is missing
    site_idx: Optional[int] = None
    hit_probability: Optional[float] = None
    raw: dict = field(default_factory=dict)

    @property
    def pose(self) -> Optional[Path]:
        return self.build.path if self.build else None


@dataclass
class DatasetOutputs:
    dtag: str
    directory: Path
    apo_model: Optional[Path]
    zmap: Optional[Path]
    mean_map: Optional[Path]
    pandda_model: Optional[Path]
    events: List[EventRecord]
    events_table_present: bool
    #: The component code of the dictionary PanDDA was given (``MZ0``), read
    #: from the copy it keeps under ligand_files/; None when there was none.
    ligand_id: Optional[str] = None
    #: That dictionary itself, so the receipt can carry it with the poses.
    dictionary: Optional[Path] = None
    #: What PanDDA linked but this machine cannot follow, by name, with the
    #: target it recorded. See ``dangling_target``.
    unresolved: Dict[str, str] = field(default_factory=dict)
    #: What was taken from the staged input tree instead, by name.
    from_staging: Dict[str, Path] = field(default_factory=dict)

    # -- what was declared vs what arrived -------------------------------
    @property
    def n_events(self) -> int:
        return len(self.events)

    @property
    def n_event_maps(self) -> int:
        return sum(1 for e in self.events if e.event_map is not None)

    @property
    def n_poses_expected(self) -> int:
        return sum(1 for e in self.events if e.build is not None)

    @property
    def n_poses(self) -> int:
        return sum(1 for e in self.events if e.pose is not None)

    def shortfalls(self) -> List[str]:
        """Every declared thing that did not arrive, in words. Empty means
        the receipt is complete."""
        missing = []
        if self.apo_model is None:
            missing.append(self._why("apo model (-pandda-input.pdb)", "apo_model"))
        if self.zmap is None:
            missing.append(self._why("Z-map", "zmap"))
        for event in self.events:
            if event.event_map is None:
                missing.append(f"event {event.idx}: event map")
            if event.build is not None and event.build.path is None:
                missing.append(f"event {event.idx}: candidate pose")
        return missing

    def _why(self, what: str, key: str) -> str:
        """A shortfall, said precisely: absent, or linked to somewhere this
        machine cannot reach."""
        target = self.unresolved.get(key)
        if not target:
            return what
        return (f"{what}: PanDDA linked it to {target}, which does not exist "
                "here -- the run used a different path for the share")


def _float(value) -> Optional[float]:
    try:
        return None if value is None else float(value)
    except (TypeError, ValueError):
        return None


def dangling_target(path: Path) -> Optional[str]:
    """Where a symlink points, when it points nowhere reachable from here.

    None when ``path`` is not a symlink, or is one that resolves. This is the
    difference between "PanDDA did not write it" and "PanDDA wrote a link this
    machine cannot follow", and the two want different answers: the first is a
    failed run, the second is a mount that does not match the one the run
    used. Reported as "missing", the second sends a reader looking for a file
    that is there.

    It happens whenever PanDDA runs somewhere the share sits at a different
    absolute path -- an Azure Batch node mounting it under
    $AZ_BATCH_NODE_MOUNTS_DIR while the job's own machine has it at
    /mnt/projects. PanDDA links its input model and its ligand dictionary into
    its output tree using the paths it was given, so both links carry the
    other machine's prefix.
    """
    try:
        if not path.is_symlink():
            return None
        if path.exists():           # follows the link
            return None
        return os.readlink(path)
    except OSError:
        return None


def _existing(path: Path) -> Optional[Path]:
    """``path`` if it names a real file (through any symlink), else None.

    PanDDA's ``-pandda-input.pdb`` is a symlink into the staging tree; the
    receipt must copy what it points at, and a dangling link is a missing
    file, not a present one. ``dangling_target`` above says which of the two
    happened, and ``read_dataset`` can fall back to the staged original.
    """
    try:
        resolved = path.resolve(strict=True)
    except (OSError, RuntimeError):
        return None
    return resolved if resolved.is_file() else None


def read_events_yaml(path: Path) -> Dict[int, dict]:
    """``events.yaml`` as ``{event_idx: record}``; ``{}`` for no events."""
    import yaml

    with open(path) as handle:
        data = yaml.safe_load(handle)
    if not data:
        return {}
    if not isinstance(data, dict):
        raise ValueError(f"{path}: expected a mapping of event index to record")
    return {int(k): (v or {}) for k, v in data.items()}


def read_events_table(tree_root: Path) -> Dict[tuple, dict]:
    """The run-level events table as ``{(dtag, event_idx): row}``, or ``{}``
    when the run did not get as far as writing it (a partial tree)."""
    path = Path(tree_root) / ANALYSES_DIR / EVENTS_TABLE
    if not path.is_file():
        return {}
    rows = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            try:
                key = (row["dtag"], int(row["event_idx"]))
            except (KeyError, ValueError):
                continue
            rows[key] = row
    return rows


#: Multiple of an event map's non-zero spread at which to open it.
DISPLAY_SIGMA = 1.5


def display_contour(event_map: Path, optimal: Optional[float]) -> Optional[float]:
    """A level to open an event map at, in absolute map units.

    An event map is a box that is flat zero outside the event, so its
    whole-box rmsd is meaningless (Reinspect found ~0.13 for BAZ2B and a
    sigma slider pinned at the rail). The spread over the non-zero region is
    the map's own scale; 1.5 of it shows the event without the noise. The
    optimal contour caps it (a build that scored best low is a low event)
    but never raises it: on a poorly characterised run it can sit above the
    map's peak.
    """
    try:
        import gemmi
        import numpy as np
        grid = np.array(gemmi.read_ccp4_map(str(event_map)).grid, copy=False)
        nonzero = grid[grid != 0]
        spread = float(nonzero.std()) if nonzero.size > 1 else 0.0
    except Exception:      # noqa: BLE001 - no level is better than a crash
        return optimal
    level = DISPLAY_SIGMA * spread if spread > 0 else None
    if level is None:
        return optimal
    if optimal is not None and optimal > 0:
        level = min(level, optimal)
    return level


def find_event_map(dataset_dir: Path, dtag: str, idx: int) -> Optional[Path]:
    """The event map for ``idx``, matched on the index and never on the BDC
    token in the name."""
    matches = sorted(dataset_dir.glob(f"{dtag}-event_{idx}_1-BDC_*_map.native.ccp4"))
    for candidate in matches:
        found = _existing(candidate)
        if found is not None:
            return found
    return None


def find_pose(dataset_dir: Path, dtag: str, idx: int, build: dict) -> Optional[Path]:
    """The selected pose: the ``best_autobuild`` copy PanDDA leaves at the top
    of the dataset directory, else the ``Build Path`` it recorded."""
    best = _existing(dataset_dir / f"{dtag}_event_{idx}_best_autobuild.pdb")
    if best is not None:
        return best
    recorded = build.get("Build Path")
    if recorded:
        return _existing(Path(recorded))
    return None


def find_dictionary(dataset_dir: Path) -> Optional[Path]:
    """The dictionary PanDDA copied under ligand_files/, if any."""
    ligand_dir = Path(dataset_dir) / LIGAND_FILES_DIR
    if not ligand_dir.is_dir():
        return None
    for path in sorted(ligand_dir.glob("*.cif")):
        if path.is_file():
            return path
    return None


def read_ligand_id(dataset_dir: Path) -> Optional[str]:
    """The component code of the dictionary PanDDA used for this dataset.

    PanDDA copies the dictionary it was given into ``ligand_files/``. The
    code is the name of the restraint block (``comp_MZ0`` -> ``MZ0``), taking
    the first block after the ``comp_list`` header that is not the
    ``comp_LIG`` alias staging appends for the name-based reader.
    """
    import gemmi

    ligand_dir = Path(dataset_dir) / LIGAND_FILES_DIR
    if not ligand_dir.is_dir():
        return None
    for path in sorted(ligand_dir.glob("*.cif")):
        code = _ligand_id_from(path)
        if code:
            return code
    return None


def _ligand_id_from(path) -> Optional[str]:
    """That component code, from one dictionary file.

    Split out because the dictionary is not always the copy under
    ligand_files/: when PanDDA's copy is a link this machine cannot follow,
    the receipt reads the staged original instead, and the code has to come
    from whichever file it actually used.
    """
    import gemmi

    if path is None:
        return None
    try:
        doc = gemmi.cif.read(str(path))
    except Exception:      # noqa: BLE001 - an unreadable copy is no code
        return None
    blocks = [b for b in doc if b.name != "comp_list"
              and b.find_values("_chem_comp_atom.atom_id")]
    preferred = [b for b in blocks if b.name != f"comp_{PANDDA_RESIDUE_NAME}"] or blocks
    if not preferred:
        return None
    name = preferred[0].name
    return name[len("comp_"):] if name.startswith("comp_") else name


def _centroid(record: dict) -> Optional[tuple]:
    value = record.get("Centroid")
    if isinstance(value, (list, tuple)) and len(value) == 3:
        coords = [_float(v) for v in value]
        if None not in coords:
            return tuple(coords)
    return None


def read_dataset(tree_root, dtag: str, staged_dir=None) -> DatasetOutputs:
    """Everything the tree holds for ``dtag``, and everything it says it
    should hold. Raises ``DatasetNotFound`` when there is no dataset
    directory at all; every lesser absence is a shortfall, reported.

    ``staged_dir`` is this dataset's directory in the tree the campaign
    staged -- the record of what PanDDA was given. Two of the things a
    receipt wants are not written by PanDDA but linked by it from there: the
    apo model and the ligand dictionary. When those links do not resolve (a
    run on a machine that mounts the share elsewhere), the staged originals
    are the same files, present and readable, so the receipt takes them and
    says it did. Without it a whole campaign harvests with no reference
    coordinates and no dictionary, which is what happened on DDU on
    2026-09-27.
    """
    tree_root = Path(tree_root)
    dataset_dir = tree_root / PROCESSED_DIR / dtag
    if not dataset_dir.is_dir():
        raise DatasetNotFound(f"{dataset_dir} is not a directory")

    table = read_events_table(tree_root)
    yaml_path = dataset_dir / EVENTS_YAML
    records = read_events_yaml(yaml_path) if yaml_path.is_file() else {}
    if not yaml_path.is_file():
        logger.warning("%s: no %s; treating as zero events", dataset_dir, EVENTS_YAML)

    events = []
    for idx in sorted(records):
        record = records[idx]
        build_raw = record.get("Build") or None
        build = None
        if build_raw:
            build = BuildRecord(
                path=find_pose(dataset_dir, dtag, idx, build_raw),
                build_score=_float(build_raw.get("Build Score", build_raw.get("Score"))),
                rscc=_float(build_raw.get("RSCC")),
                optimal_contour=_float(build_raw.get("Optimal Contour")),
                raw=build_raw,
            )
        row = table.get((dtag, idx), {})
        site_idx = row.get("site_idx")
        events.append(EventRecord(
            idx=idx,
            bdc=_float(record.get("BDC")),
            score=_float(record.get("Score")),
            centroid=_centroid(record),
            build=build,
            event_map=find_event_map(dataset_dir, dtag, idx),
            site_idx=int(site_idx) if site_idx not in (None, "") else None,
            hit_probability=_float(row.get("hit_in_site_probability")),
            raw=record,
        ))

    apo_link = dataset_dir / f"{dtag}-pandda-input.pdb"
    apo_model = _existing(apo_link)
    dictionary = find_dictionary(dataset_dir)
    unresolved, from_staging = {}, {}

    target = dangling_target(apo_link)
    if target:
        unresolved["apo_model"] = target
    for path in sorted((dataset_dir / LIGAND_FILES_DIR).glob("*.cif")) if (
            dataset_dir / LIGAND_FILES_DIR).is_dir() else []:
        target = dangling_target(path)
        if target:
            unresolved["dictionary"] = target
            break

    staged = Path(staged_dir) if staged_dir else None
    if staged is not None and staged.is_dir():
        if apo_model is None:
            fallback = _existing(staged / STAGED_MODEL_NAME)
            if fallback is not None:
                apo_model, from_staging["apo_model"] = fallback, fallback
        if dictionary is None:
            for candidate in (staged / STAGED_LIGAND_DIR / STAGED_DICT_NAME,
                              staged / STAGED_DICT_NAME):
                fallback = _existing(candidate)
                if fallback is not None:
                    dictionary, from_staging["dictionary"] = fallback, fallback
                    break

    return DatasetOutputs(
        dtag=dtag,
        directory=dataset_dir,
        apo_model=apo_model,
        zmap=_existing(dataset_dir / f"{dtag}-z_map.native.ccp4"),
        mean_map=_existing(dataset_dir / f"{dtag}-ground-state-average-map.native.ccp4"),
        pandda_model=_existing(dataset_dir / MODELLED_DIR / f"{dtag}-pandda-model.pdb"),
        events=events,
        events_table_present=bool(table),
        ligand_id=read_ligand_id(dataset_dir) or _ligand_id_from(dictionary),
        dictionary=dictionary,
        unresolved=unresolved,
        from_staging=from_staging,
    )
