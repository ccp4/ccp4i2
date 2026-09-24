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
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(f"ccp4i2:{__name__}")

ANALYSES_DIR = "analyses"
EVENTS_TABLE = "pandda_analyse_events.csv"
PROCESSED_DIR = "processed_datasets"
EVENTS_YAML = "events.yaml"
MODELLED_DIR = "modelled_structures"


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
            missing.append("apo model (-pandda-input.pdb)")
        if self.zmap is None:
            missing.append("Z-map")
        for event in self.events:
            if event.event_map is None:
                missing.append(f"event {event.idx}: event map")
            if event.build is not None and event.build.path is None:
                missing.append(f"event {event.idx}: candidate pose")
        return missing


def _float(value) -> Optional[float]:
    try:
        return None if value is None else float(value)
    except (TypeError, ValueError):
        return None


def _existing(path: Path) -> Optional[Path]:
    """``path`` if it names a real file (through any symlink), else None.

    PanDDA's ``-pandda-input.pdb`` is a symlink into the staging tree; the
    receipt must copy what it points at, and a dangling link is a missing
    file, not a present one.
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


def _centroid(record: dict) -> Optional[tuple]:
    value = record.get("Centroid")
    if isinstance(value, (list, tuple)) and len(value) == 3:
        coords = [_float(v) for v in value]
        if None not in coords:
            return tuple(coords)
    return None


def read_dataset(tree_root, dtag: str) -> DatasetOutputs:
    """Everything the tree holds for ``dtag``, and everything it says it
    should hold. Raises ``DatasetNotFound`` when there is no dataset
    directory at all; every lesser absence is a shortfall, reported."""
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

    return DatasetOutputs(
        dtag=dtag,
        directory=dataset_dir,
        apo_model=_existing(dataset_dir / f"{dtag}-pandda-input.pdb"),
        zmap=_existing(dataset_dir / f"{dtag}-z_map.native.ccp4"),
        mean_map=_existing(dataset_dir / f"{dtag}-ground-state-average-map.native.ccp4"),
        pandda_model=_existing(dataset_dir / MODELLED_DIR / f"{dtag}-pandda-model.pdb"),
        events=events,
        events_table_present=bool(table),
    )
