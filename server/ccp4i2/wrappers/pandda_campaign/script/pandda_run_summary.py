"""Recapitulate a PanDDA run's own analysis from the files it wrote.

PanDDA 1 rendered an HTML summary (events per site, histograms of event
fraction, resolution and R-free); the CCP4-bundled PanDDA 2 writes none
(``analyses/html_summaries/`` stays empty). Reinspect recovered the same
picture from the tables PanDDA leaves behind, and this module does the same
for the campaign report:

* ``analyses/pandda_analyse_events.csv``: one row per event (dtag, event_idx,
  site_idx, BDC and ``1-BDC``, cluster_size, z_peak, z_mean, x/y/z,
  hit_in_site_probability, r_free, r_work, analysed_resolution, interesting).
* ``analyses/pandda_analyse_sites.csv``: one row per site (site_idx,
  centroid as a ``"(x, y, z)"`` string).
* ``processed_datasets/<dtag>/processed_dataset.yaml``: per dataset the
  processing resolution, comparator set, selected model and, per model, the
  event counts at each filtering stage. The file carries numpy scalars as
  ``!!python/object/apply`` tags, so ``yaml.safe_load`` refuses it; the
  loader here accepts the tags and yields ``None`` for what they wrap.
* ``processed_datasets/<dtag>/events.yaml``: the selected model's events with
  build score, RSCC and optimal contour (read through ``pandda_tree``).

Everything returned is plain Python (dicts, lists, floats) so it can be
written to program.xml and unit-tested against a synthetic tree without
PanDDA.
"""
from __future__ import annotations

import csv
import math
import re
from pathlib import Path
from typing import Dict, List, Optional
from xml.etree import ElementTree as ET

from ccp4i2.wrappers.pandda_events.script.pandda_tree import (
    ANALYSES_DIR, EVENTS_TABLE, EVENTS_YAML, PROCESSED_DIR, read_events_yaml,
)

SITES_TABLE = "pandda_analyse_sites.csv"
DATASET_YAML = "processed_dataset.yaml"

#: Bins for the histograms; Reinspect's default, wide enough for a campaign
#: and readable for a handful of events.
HISTOGRAM_BINS = 10


# --- PanDDA's YAML -----------------------------------------------------------

def load_pandda_yaml(path):
    """``processed_dataset.yaml`` as plain Python.

    numpy scalars are serialised as ``!!python/object/apply:numpy...`` with
    the value in a binary blob; reconstructing them would need numpy's
    unpickling, and nothing here needs the numbers they wrap (event scores are
    re-read from ``events.yaml``, which writes them as floats). They load as
    ``None``.
    """
    import yaml

    class Tolerant(yaml.SafeLoader):
        pass

    def unknown(loader, tag_suffix, node):
        if isinstance(node, yaml.SequenceNode):
            loader.construct_sequence(node, deep=True)
        elif isinstance(node, yaml.MappingNode):
            loader.construct_mapping(node, deep=True)
        return None

    Tolerant.add_multi_constructor("tag:yaml.org,2002:python/", unknown)
    with open(path) as handle:
        return yaml.load(handle, Loader=Tolerant)


def _float(value) -> Optional[float]:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(result) else result


def _int(value) -> Optional[int]:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _bool(value) -> Optional[bool]:
    if isinstance(value, bool):
        return value
    text = str(value or "").strip().lower()
    if text in ("true", "1", "yes"):
        return True
    if text in ("false", "0", "no"):
        return False
    return None


# --- per-dataset -------------------------------------------------------------

def read_dataset_summary(dataset_dir) -> Dict[str, object]:
    """What ``processed_dataset.yaml`` says about one dataset, or an empty
    record when the file is absent or unreadable (a dataset PanDDA loaded
    but never characterised)."""
    dataset_dir = Path(dataset_dir)
    record = {
        "dtag": dataset_dir.name,
        "resolution": None,
        "n_comparators": None,
        "selected_model": None,
        "n_models": None,
        "n_events": None,
        "n_initial_events": None,
        "n_size_filtered_events": None,
        "n_score_filtered_events": None,
        "analysed": (dataset_dir / f"{dataset_dir.name}-z_map.native.ccp4").is_file(),
    }
    path = dataset_dir / DATASET_YAML
    if not path.is_file():
        return record
    try:
        data = load_pandda_yaml(path) or {}
    except Exception:       # noqa: BLE001 - a summary is never worth failing over
        return record
    summary = data.get("Summary") or {}
    record["resolution"] = _float(summary.get("Processing Resolution"))
    comparators = summary.get("Comparator Datasets")
    record["n_comparators"] = len(comparators) if isinstance(comparators, list) else None
    record["selected_model"] = _int(summary.get("Selected Model"))
    selected_events = summary.get("Selected Model Events")
    record["n_events"] = len(selected_events) if isinstance(selected_events, list) else None
    models = data.get("Models") or {}
    if isinstance(models, dict):
        record["n_models"] = len(models)
        selected = models.get(record["selected_model"]) or models.get(str(record["selected_model"]))
        if isinstance(selected, dict):
            record["n_initial_events"] = _int(selected.get("Number of Initial Events"))
            record["n_size_filtered_events"] = _int(selected.get("Number of Size Filtered Events"))
            record["n_score_filtered_events"] = _int(selected.get("Number of Score Filtered Events"))
    return record


# --- run-level tables --------------------------------------------------------

EVENT_FLOATS = {
    "bdc": "bdc", "cluster_size": "cluster_size", "z_peak": "z_peak", "z_mean": "z_mean",
    "x": "x", "y": "y", "z": "z", "hit_probability": "hit_in_site_probability",
    "r_free": "r_free", "r_work": "r_work", "resolution": "analysed_resolution",
    "high_resolution": "high_resolution",
}


def read_events(tree) -> List[Dict[str, object]]:
    """The run's events table, typed, joined with each dataset's
    ``events.yaml`` for the build columns (score, build score, RSCC,
    optimal contour). ``[]`` when the table was never written."""
    tree = Path(tree)
    path = tree / ANALYSES_DIR / EVENTS_TABLE
    if not path.is_file():
        return []
    builds: Dict[str, Dict[int, dict]] = {}
    events = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            dtag = row.get("dtag")
            idx = _int(row.get("event_idx"))
            if not dtag or idx is None:
                continue
            record = {"dtag": dtag, "event_idx": idx, "site_idx": _int(row.get("site_idx")),
                      "interesting": _bool(row.get("interesting"))}
            for key, column in EVENT_FLOATS.items():
                record[key] = _float(row.get(column))
            if record["bdc"] is not None:
                record["event_fraction"] = round(1.0 - record["bdc"], 4)
            else:
                record["event_fraction"] = _float(row.get("1-BDC"))
            if dtag not in builds:
                builds[dtag] = _dataset_builds(tree / PROCESSED_DIR / dtag)
            build = builds[dtag].get(idx) or {}
            record["score"] = _float(build.get("Score", build.get("score")))
            inner = build.get("Build") or {}
            record["build_score"] = _float(inner.get("Build Score"))
            record["rscc"] = _float(inner.get("RSCC"))
            record["optimal_contour"] = _float(inner.get("Optimal Contour"))
            events.append(record)
    events.sort(key=lambda e: (e["dtag"], e["event_idx"]))
    return events


def _dataset_builds(dataset_dir: Path) -> Dict[int, dict]:
    path = dataset_dir / EVENTS_YAML
    if not path.is_file():
        return {}
    try:
        return read_events_yaml(path)
    except Exception:       # noqa: BLE001
        return {}


CENTROID_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def read_sites(tree) -> List[Dict[str, object]]:
    """The sites table as ``[{site_idx, centroid}]``; ``[]`` when absent."""
    path = Path(tree) / ANALYSES_DIR / SITES_TABLE
    if not path.is_file():
        return []
    sites = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            idx = _int(row.get("site_idx"))
            if idx is None:
                continue
            numbers = [float(m) for m in CENTROID_RE.findall(row.get("centroid") or "")]
            sites.append({"site_idx": idx, "centroid": tuple(numbers[:3]) if len(numbers) >= 3 else None})
    sites.sort(key=lambda s: s["site_idx"])
    return sites


# --- aggregation -------------------------------------------------------------

def bin_values(values, nbins: int = HISTOGRAM_BINS) -> List[Dict[str, float]]:
    """Equal-width bins over ``values`` as ``[{centre, count}]``, the shape a
    bar chart takes. A single distinct value gives one bin."""
    values = [v for v in values if v is not None]
    if not values:
        return []
    lo, hi = min(values), max(values)
    if hi == lo:
        return [{"centre": lo, "count": len(values)}]
    width = (hi - lo) / nbins
    counts = [0] * nbins
    for v in values:
        counts[min(nbins - 1, int((v - lo) / width))] += 1
    return [{"centre": round(lo + width * (i + 0.5), 4), "count": counts[i]} for i in range(nbins)]


def summarise_run(tree) -> Dict[str, object]:
    """Everything the report shows about a run, from the tree alone."""
    tree = Path(tree)
    processed = tree / PROCESSED_DIR
    dirs = sorted(d for d in processed.iterdir() if d.is_dir()) if processed.is_dir() else []
    datasets = [read_dataset_summary(d) for d in dirs]
    events = read_events(tree)
    by_dataset: Dict[str, List[dict]] = {}
    for event in events:
        by_dataset.setdefault(event["dtag"], []).append(event)
    for dataset in datasets:
        own = by_dataset.get(dataset["dtag"], [])
        if own and dataset["n_events"] is None:
            dataset["n_events"] = len(own)
        # r_free / r_work / high resolution are only ever written per event
        first = own[0] if own else {}
        dataset["r_free"] = first.get("r_free")
        dataset["r_work"] = first.get("r_work")
        dataset["best_score"] = max((e["score"] for e in own if e["score"] is not None), default=None)
        dataset["best_hit_probability"] = max(
            (e["hit_probability"] for e in own if e["hit_probability"] is not None), default=None)

    sites = {s["site_idx"]: dict(s, n_events=0, n_datasets=0, best_score=None,
                                 best_hit_probability=None, n_interesting=0)
             for s in read_sites(tree)}
    site_datasets: Dict[int, set] = {}
    for event in events:
        idx = event["site_idx"]
        if idx is None:
            continue
        site = sites.setdefault(idx, {"site_idx": idx, "centroid": None, "n_events": 0, "n_datasets": 0,
                                      "best_score": None, "best_hit_probability": None, "n_interesting": 0})
        site["n_events"] += 1
        site_datasets.setdefault(idx, set()).add(event["dtag"])
        if event["interesting"]:
            site["n_interesting"] += 1
        for key in ("score", "hit_probability"):
            value = event[key]
            best = f"best_{key}"
            if value is not None and (site[best] is None or value > site[best]):
                site[best] = value
    for idx, site in sites.items():
        site["n_datasets"] = len(site_datasets.get(idx, ()))

    analysed = [d for d in datasets if d["analysed"]]
    stats = {
        "n_datasets": len(datasets),
        "n_analysed": len(analysed),
        "n_events": len(events),
        "n_sites": len(sites),
        "n_datasets_with_events": len(by_dataset),
        "n_interesting": sum(1 for e in events if e["interesting"]),
        "best_hit_probability": max((e["hit_probability"] for e in events
                                     if e["hit_probability"] is not None), default=None),
        "best_score": max((e["score"] for e in events if e["score"] is not None), default=None),
        "median_resolution": _median([d["resolution"] for d in analysed]),
    }
    histograms = {
        "event_fraction": bin_values([e["event_fraction"] for e in events]),
        "hit_probability": bin_values([e["hit_probability"] for e in events]),
        "resolution": bin_values([d["resolution"] for d in analysed]),
        "r_free": bin_values([d["r_free"] for d in datasets]),
    }
    return {"datasets": datasets, "events": events,
            "sites": [sites[k] for k in sorted(sites)], "stats": stats, "histograms": histograms}


def _median(values) -> Optional[float]:
    values = sorted(v for v in values if v is not None)
    if not values:
        return None
    mid = len(values) // 2
    return values[mid] if len(values) % 2 else round((values[mid - 1] + values[mid]) / 2, 4)


# --- program.xml -------------------------------------------------------------

def _text(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.4g}" if abs(value) < 1e-3 or abs(value) >= 1e4 else f"{value:.4f}".rstrip("0").rstrip(".")
    return str(value)


DATASET_ATTRS = ("dtag", "analysed", "resolution", "n_comparators", "n_models", "selected_model",
                 "n_initial_events", "n_size_filtered_events", "n_score_filtered_events", "n_events",
                 "r_free", "r_work", "best_score", "best_hit_probability")
EVENT_ATTRS = ("dtag", "event_idx", "site_idx", "bdc", "event_fraction", "cluster_size", "z_peak", "z_mean",
               "score", "hit_probability", "build_score", "rscc", "optimal_contour", "resolution",
               "r_free", "r_work", "interesting", "x", "y", "z")
SITE_ATTRS = ("site_idx", "n_events", "n_datasets", "n_interesting", "best_score", "best_hit_probability")


def analysis_to_xml(parent: ET.Element, analysis: Dict[str, object]) -> ET.Element:
    """Append ``<analysis>`` to ``parent``: stats as children, then one
    element per dataset, event, site and histogram bin, values as attributes
    (absent ones written as empty strings so the columns line up)."""
    root = ET.SubElement(parent, "analysis")
    stats = ET.SubElement(root, "stats")
    for key, value in analysis["stats"].items():
        ET.SubElement(stats, key).text = _text(value)
    datasets = ET.SubElement(root, "datasets")
    for record in analysis["datasets"]:
        ET.SubElement(datasets, "dataset", {k: _text(record.get(k)) for k in DATASET_ATTRS})
    events = ET.SubElement(root, "events")
    for record in analysis["events"]:
        ET.SubElement(events, "event", {k: _text(record.get(k)) for k in EVENT_ATTRS})
    sites = ET.SubElement(root, "sites")
    for record in analysis["sites"]:
        attrs = {k: _text(record.get(k)) for k in SITE_ATTRS}
        centroid = record.get("centroid")
        attrs["centroid"] = " ".join(f"{c:.2f}" for c in centroid) if centroid else ""
        ET.SubElement(sites, "site", attrs)
    histograms = ET.SubElement(root, "histograms")
    for name, bins in analysis["histograms"].items():
        element = ET.SubElement(histograms, "histogram", name=name)
        for entry in bins:
            ET.SubElement(element, "bin", centre=_text(entry["centre"]), count=str(entry["count"]))
    return root
