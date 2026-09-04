"""
    dnatco_report.py: CCP4 GUI Project
    Copyright (C) 2025 MRC-LMB
    Author: Martin Maly

    This library is free software: you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    version 3, modified in accordance with the provisions of the
    license to address the requirements of UK law.

    You should have received a copy of the modified GNU Lesser General
    Public License along with this library.  If not, copies may be
    downloaded from http://www.ccp4.ac.uk/ccp4license.php

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

The report is drawn by :func:`draw_dnatco_report` from the files DNATCO
wrote (extended mmCIF + NAVAL JSON), for one model or for two side by side.
``dnatco_pipe_report`` reuses it for the pipeline's one-or-two models.
"""

import os

from ccp4i2.core import CCP4Modules
from ccp4i2.report import Report
from ccp4i2.wrappers.dnatco.script import dnatco_data

DNATCO_SERVER_NOTE = (
    "A more detailed analysis can be performed at the"
    " <a href='https://dnatco.datmos.org'>DNATCO web server</a>."
)
MAX_CONCERN_ROWS = 100
GRAPH_STYLE = "height:330px; width:585px; border:0px; padding:10px; padding-left:15px; margin-right:15px;"


def _fmt(value, digits=2):
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def _load_models(cif_paths, json_paths):
    """One entry per model: its parsed NtC data, NAVAL data, paths and a label."""
    count = max(len(cif_paths), len(json_paths))
    models = []
    for index in range(count):
        cif_path = cif_paths[index] if index < len(cif_paths) else None
        json_path = json_paths[index] if index < len(json_paths) else None
        models.append({
            "label": f"model {index + 1}" if count > 1 else "",
            "cif_path": cif_path,
            "json_path": json_path,
            "ntc": dnatco_data.read_ntc_data(cif_path),
            "naval": dnatco_data.read_naval_json(json_path),
        })
    return models


def _suffix(model):
    return f" ({model['label']})" if model["label"] else ""


def _expected(path):
    return f" (expected path: {path})" if path else ""


def _merged_steps(models):
    """Steps of all models aligned by CIF step name, in first-model order.

    Each entry: {"name", "chain", "step", "per_model": [step-or-None, ...]}.
    """
    merged = []
    by_name = {}
    for model_index, model in enumerate(models):
        ntc = model["ntc"]
        if ntc is None:
            continue
        for step in ntc["steps"]:
            key = step["name"] or f"{step['chain']} {step['step']}"
            entry = by_name.get(key)
            if entry is None:
                entry = {"name": key, "chain": step["chain"], "step": step["step"],
                         "per_model": [None] * len(models)}
                by_name[key] = entry
                merged.append(entry)
            entry["per_model"][model_index] = step
    return merged


def draw_dnatco_report(parent, cif_paths, json_paths):
    """Draw the DNATCO report for one model, or two models compared.

    ``cif_paths`` and ``json_paths`` are parallel lists (one entry per model)
    of the extended mmCIF and NAVAL JSON files DNATCO wrote.
    """
    models = _load_models(cif_paths, json_paths)
    compare = len(models) > 1

    for model in models:
        if model["ntc"] is None:
            note = parent.addDiv(style="font-size:110%;color:red;")
            note.append(
                f"DNATCO did not produce an extended mmCIF with NtC assignments{_suffix(model)}"
                f"{_expected(model['cif_path'])}. Please check the log file.")
        if model["naval"] is None:
            note = parent.addDiv(style="font-size:110%;color:red;")
            note.append(
                f"DNATCO did not produce the NAVAL validation file{_suffix(model)}"
                f"{_expected(model['json_path'])}. Please check the log file.")

    if any(model["naval"] is not None for model in models):
        _draw_naval(parent, models, compare)
    if any(model["ntc"] is not None for model in models):
        _draw_ntc(parent, models, compare)
    parent.addDiv(style="clear:both;")


def _draw_naval(parent, models, compare):
    fold = parent.addFold(label="NAVAL bond lengths and angles validation", initiallyOpen=True)
    body = fold.addDiv(style="margin-left:1.5em;")
    note = body.addDiv(style="font-size:110%;")
    note.append(
        "Every bond length and bond angle of the nucleotides is compared with the NAVAL"
        " reference distributions and classed as Preferred, Allowed or Of Concern.")

    table = body.addTable(title="NAVAL overall validation", transpose=True)
    header = []
    for model in models:
        header += [f"Bond lengths{_suffix(model)}", f"Bond angles{_suffix(model)}"]
    table.addData(title="", data=header)
    for tier in dnatco_data.NAVAL_TIERS:
        row = []
        for model in models:
            for kind in ("lengths", "angles"):
                count, percentage = dnatco_data.naval_tier_counts(model["naval"], kind).get(tier, (None, None))
                row.append("-" if count is None else f"{count} ({_fmt(percentage, 1)}%)")
        table.addData(title=tier, data=row)

    # Geometry terms worth a look are listed for the last model only: in a
    # before/after comparison that is the model being judged.
    model = models[-1]
    if model["naval"] is None:
        return
    for title, key, name_label, value_label in [
        ("Bond lengths", "lengths", "Bond", "Value (&#197;)"),
        ("Bond angles", "angles", "Angle", "Value (&#176;)"),
    ]:
        concerned = dnatco_data.concerned_items(model["naval"].get(key, []))
        concern_fold = body.addFold(label=f"{title} of concern{_suffix(model)}", initiallyOpen=True)
        if not concerned:
            concern_fold.append(f"No {title.lower()} of concern found.")
            continue
        concern_fold.append(
            f"{len(concerned)} questionable {title.lower()} found."
            f" Up to {MAX_CONCERN_ROWS} are listed below, sorted by NAVAL tier and then by"
            " ProSco (probability score, lower is less likely).<br />"
            "<i>pGroup is the probability percentile group;"
            " a superscript (-1) marks an atom of the preceding nucleotide.</i>")
        concerned = concerned[:MAX_CONCERN_ROWS]
        table = concern_fold.addTable(title=f"{title} with concerns")
        table.addData(title="Residue", data=[item["residue"] for item in concerned])
        table.addData(title=name_label, data=[item["name"] for item in concerned])
        table.addData(title=value_label, data=[_fmt(item["value"]) for item in concerned])
        table.addData(title="NAVAL", data=[item["naval_tier"] for item in concerned])
        table.addData(title="ProSco", data=[_fmt(item["prosco"]) for item in concerned])
        table.addData(title="pGroup", data=[item["pGroup"] for item in concerned])


def _draw_ntc(parent, models, compare):
    fold = parent.addFold(label="DNATCO dinucleotide conformer class validation", initiallyOpen=True)
    body = fold.addDiv(style="margin-left:1.5em;")
    merged = _merged_steps(models)

    overall_fold = body.addFold(label="Overall structure quality", initiallyOpen=True)
    note = overall_fold.addDiv(style="font-size:110%;")
    note.append(
        "Assignment of a conformer class (NtC) to each dinucleotide step, and the quality"
        " of fit between the step and the NtC reference structure measured by RMSD and"
        " confal score.")
    table = overall_fold.addTable(title="Overall structure quality", transpose=True)
    if compare:
        table.addData(title="", data=[model["label"].capitalize() for model in models])
    rows = [
        ("No. NtC assigned", lambda o, s: o["num_classified"]),
        ("No. NtC close", lambda o, s: o["num_unclassified_rmsd_close"]),
        ("No. NtC unassigned", lambda o, s: o["num_unclassified"]),
        ("No. steps with RMSD &#8804; 0.5 &#197;", lambda o, s: dnatco_data.rmsd_bins(s)[0]),
        ("No. steps with RMSD 0.5-1 &#197;", lambda o, s: dnatco_data.rmsd_bins(s)[1]),
        ("No. steps with RMSD &gt; 1 &#197;", lambda o, s: dnatco_data.rmsd_bins(s)[2]),
        ("Confal score", lambda o, s: o["confal_score"]),
        ("Confal score percentile", lambda o, s: o["confal_percentile"]),
    ]
    for title, getter in rows:
        table.addData(title=title, data=[
            "-" if model["ntc"] is None else _fmt(getter(model["ntc"]["overall"], model["ntc"]["steps"]))
            for model in models])

    if merged:
        graph = overall_fold.addFlotGraph(title="RMSD per step", style=GRAPH_STYLE)
        graph.addData(title="step", data=list(range(1, len(merged) + 1)))
        for index, model in enumerate(models):
            graph.addData(
                title=f"RMSD{'_' + model['label'].replace(' ', '') if model['label'] else ''}(A)",
                data=[None if entry["per_model"][index] is None else entry["per_model"][index]["rmsd"]
                      for entry in merged])
        plot = graph.addPlotObject()
        plot.append("title", "RMSD to the closest NtC representative, per step")
        plot.append("plottype", "xy")
        plot.append("xlabel", "step")
        plot.append("ylabel", "RMSD (A)")
        plot.append("legendposition", x=1, y=0)
        for index, colour in zip(range(len(models)), ("gray", "blue")):
            line = plot.append("plotline", xcol=1, ycol=index + 2)
            line.append("symbolsize", "1")
            line.append("linestyle", ".")
            line.append("colour", colour)

    _draw_step_table(
        body, models, merged, compare,
        label="Dinucleotide outliers", flag="is_outlier",
        description=(
            "All unassigned dinucleotide steps. The conformer (NtC) listed is the closest"
            " NtC that would be assigned to the step if all assignment criteria were met."),
        none_text="No dinucleotide outliers found.")
    _draw_step_table(
        body, models, merged, compare,
        label="Improvable dinucleotide outliers", flag="is_improvable",
        description=(
            "Unassigned dinucleotide steps that are sufficiently close to a representative"
            f" from the Golden Set: RMSD at most {dnatco_data.IMPROVABLE_RMSD} &#197;."),
        none_text="No improvable dinucleotide outliers found.")

    all_fold = body.addFold(label="All dinucleotides", initiallyOpen=False)
    table = all_fold.addTable(title="All dinucleotides")
    table.addData(title="Chain", data=[entry["chain"] for entry in merged])
    table.addData(title="Step", data=[entry["step"] for entry in merged])
    for index, model in enumerate(models):
        steps = [entry["per_model"][index] for entry in merged]
        table.addData(title=f"Assigned NtC{_suffix(model)}",
                      data=[s["assigned_NtC"] or "-" if s else "-" for s in steps])
        table.addData(title=f"Confal score{_suffix(model)}",
                      data=[_fmt(s["confal_score"]) if s else "-" for s in steps])
        table.addData(title=f"RMSD to closest NtC (&#197;){_suffix(model)}",
                      data=[_closest_text(s) if s else "-" for s in steps])


def _closest_text(step):
    """'0.42 (BB00)': RMSD and, for an unassigned step, which NtC it is closest to."""
    if step["is_assigned"]:
        return _fmt(step["rmsd"])
    return f"{_fmt(step['rmsd'])} ({step['closest_NtC'] or '-'})"


def _draw_step_table(body, models, merged, compare, label, flag, description, none_text):
    fold = body.addFold(label=label, initiallyOpen=True)
    note = fold.addDiv(style="font-size:110%;")
    note.append(description)
    selected = [entry for entry in merged
                if any(step is not None and step[flag] for step in entry["per_model"])]
    if not selected:
        fold.addDiv(style="font-size:110%;").append(none_text)
        return
    if compare:
        counts = [sum(1 for entry in selected if entry["per_model"][i] is not None and entry["per_model"][i][flag])
                  for i in range(len(models))]
        note.append(" ".join(f"{count} in {model['label']}." for count, model in zip(counts, models)))
    else:
        note.append(f"{len(selected)} reported.")
    table = fold.addTable(title=label)
    table.addData(title="Chain", data=[entry["chain"] for entry in selected])
    table.addData(title="Step", data=[entry["step"] for entry in selected])
    for index, model in enumerate(models):
        steps = [entry["per_model"][index] for entry in selected]
        table.addData(title=f"Closest NtC{_suffix(model)}",
                      data=[(s["closest_NtC"] or "-") if s else "-" for s in steps])
        table.addData(title=f"RMSD to closest NtC (&#197;){_suffix(model)}",
                      data=[_fmt(s["rmsd"]) if s else "-" for s in steps])
    if compare:
        def status(entry):
            flags = [step is not None and step[flag] for step in entry["per_model"]]
            if all(flags):
                return "both"
            return ", ".join(model["label"] for model, on in zip(models, flags) if on) + " only"
        table.addData(title="Where", data=[status(entry) for entry in selected])
    fold.addDiv(style="font-size:110%;").append(DNATCO_SERVER_NOTE)


class dnatco_report(Report):
    TASKNAME = 'dnatco'
    RUNNING = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        if jobStatus is None or jobStatus.lower() == 'nooutput':
            return
        if jobStatus.lower() == 'running':
            self.runningReport(parent=self)
        else:
            self.defaultReport(parent=self)

    def runningReport(self, parent=None):
        if parent is None:
            parent = self
        fold = parent.addFold(label="DNATCO log", initiallyOpen=True)
        fold.addPre("DNATCO is running...")

    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        parent.addDiv(style="clear:both;")  # gives space for the title
        filenames = self.jobInfo.get('filenames', {}) if self.jobInfo else {}
        draw_dnatco_report(parent, [filenames.get('CIFOUT')], [filenames.get('JSONOUT')])
