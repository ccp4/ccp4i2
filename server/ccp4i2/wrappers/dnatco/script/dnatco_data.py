"""
Readers for what DNATCO writes: the NtC assignment categories it adds to an
extended mmCIF, and the NAVAL bond-length / bond-angle JSON.

Shared by the dnatco wrapper (which summarises a run into program.xml), the
dnatco and dnatco_pipe reports, and the pre-run nucleic-acid check. Pure
gemmi + stdlib so the reports can import it on the server.
"""

import json
from pathlib import Path


# Threshold below which an unassigned step is "close" to an NtC representative
# and hence considered improvable (DNATCO's own RMSD closeness criterion).
IMPROVABLE_RMSD = 0.5

NAVAL_TIERS = ("Preferred", "Allowed", "Of Concern")

# NAVAL lists geometry as either Preferred, Allowed, or Of Concern. The report
# lists everything that is not Preferred, and everything whose probability
# group is not "Common" (rare values can still sit inside the allowed range).
CONCERN_TIERS = ("Of Concern", "Allowed")
CONCERN_PGROUPS = ("Rare", "Unique", "Ambiguous", "Outlier", None)


def _num(value, cast=float):
    """A CIF value as a number, or None for '.', '?' and junk."""
    if value is None:
        return None
    text = str(value).strip()
    if text in ("", ".", "?"):
        return None
    try:
        return cast(text)
    except (TypeError, ValueError):
        return None


def _text(value):
    """A CIF value as unquoted text, '' for '.'/'?'."""
    import gemmi
    if value is None:
        return ""
    text = str(value)
    if text in (".", "?"):
        return ""
    try:
        return gemmi.cif.as_string(text)
    except Exception:
        return text.strip("'\"")


def model_has_nucleic_acid(path):
    """True if the first model has at least one nucleotide residue.

    Returns True (rather than guessing wrong) when gemmi is not importable or
    the file cannot be read: the job itself will then say what is wrong.
    """
    try:
        import gemmi
        structure = gemmi.read_structure(str(path))
    except Exception:
        return True
    if len(structure) == 0:
        return True
    for chain in structure[0]:
        for residue in chain:
            info = gemmi.find_tabulated_residue(residue.name)
            if info is not None and info.is_nucleic_acid():
                return True
    return False


def read_naval_json(path):
    """The NAVAL per-residue JSON as a dict, or None if absent/unreadable."""
    if path is None or not Path(str(path)).is_file():
        return None
    try:
        with open(str(path), "r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def read_ntc_data(path):
    """The NtC assignment categories of a DNATCO extended mmCIF.

    Returns None if the file is missing or carries no DNATCO categories.
    Otherwise a dict with the overall block and one entry per dinucleotide
    step (in file order), each step a dict keyed on the CIF step name.
    """
    if path is None or not Path(str(path)).is_file():
        return None
    import gemmi
    try:
        document = gemmi.cif.read_file(str(path))
    except Exception:
        return None
    if len(document) == 0:
        return None
    block = document[0]

    overall_prefix = "_ndb_struct_ntc_overall."
    overall = {
        "confal_score": _num(block.find_value(overall_prefix + "confal_score")),
        "confal_percentile": _num(block.find_value(overall_prefix + "confal_percentile"), int),
        "num_steps": _num(block.find_value(overall_prefix + "num_steps"), int),
        "num_classified": _num(block.find_value(overall_prefix + "num_classified"), int),
        "num_unclassified": _num(block.find_value(overall_prefix + "num_unclassified"), int),
        "num_unclassified_rmsd_close": _num(
            block.find_value(overall_prefix + "num_unclassified_rmsd_close"), int),
    }

    step_table = block.find("_ndb_struct_ntc_step.", [
        "id", "name",
        "label_asym_id_1", "auth_asym_id_1",
        "label_comp_id_1", "auth_seq_id_1",
        "label_comp_id_2", "auth_seq_id_2",
    ])
    summary_table = block.find("_ndb_struct_ntc_step_summary.", [
        "step_id", "assigned_NtC", "confal_score",
        "cartesian_rmsd_closest_NtC_representative", "closest_NtC",
    ])
    # gemmi's find() yields an empty (falsy) table when any tag is absent
    if not step_table or not summary_table:
        if all(v is None for v in overall.values()):
            return None
        return {"overall": overall, "steps": []}

    summaries = {}
    for row in summary_table:
        summaries[_text(row[0])] = {
            "assigned_NtC": _text(row[1]),
            "confal_score": _num(row[2]),
            "rmsd": _num(row[3]),
            "closest_NtC": _text(row[4]),
        }

    steps = []
    for row in step_table:
        step_id = _text(row[0])
        label_chain = _text(row[2])
        auth_chain = _text(row[3])
        summary = summaries.get(step_id, {})
        assigned = summary.get("assigned_NtC", "")
        closest = summary.get("closest_NtC", "")
        is_assigned = bool(assigned) and assigned != "NANT"
        rmsd = summary.get("rmsd")
        steps.append({
            "step_id": step_id,
            "name": _text(row[1]),
            "chain": label_chain if label_chain == auth_chain or not auth_chain
                     else f"{label_chain} (auth: {auth_chain})",
            "step": f"{_text(row[4])}{_text(row[5])} {_text(row[6])}{_text(row[7])}",
            "assigned_NtC": assigned if is_assigned else "",
            "closest_NtC": closest or (assigned if is_assigned else ""),
            "confal_score": summary.get("confal_score"),
            "rmsd": rmsd,
            "is_assigned": is_assigned,
            "is_outlier": not is_assigned,
            "is_improvable": (not is_assigned) and rmsd is not None and rmsd <= IMPROVABLE_RMSD,
        })

    return {"overall": overall, "steps": steps}


def rmsd_bins(steps):
    """How many steps have RMSD <= 0.5, in (0.5, 1], and > 1 angstrom."""
    below, between, above = 0, 0, 0
    for step in steps:
        rmsd = step.get("rmsd")
        if rmsd is None:
            continue
        if rmsd <= 0.5:
            below += 1
        elif rmsd <= 1.0:
            between += 1
        else:
            above += 1
    return below, between, above


def format_geometry_name(name):
    """C2'-C3' style names; 'r2' marks an atom of the previous residue."""
    if name is None:
        return ""
    parts = [part.strip() for part in str(name).split("-") if part and part.strip()]
    return "-".join(part.replace("r2", "<sup>(-1)</sup>") for part in parts)


def format_residue(item):
    """'A/G 12.B' from a NAVAL per-residue entry: chain/compound seqId[.ins][ alt]."""
    residue = f"{item.get('authChain', '')}/{item.get('compound', '')} {item.get('authSeqId', '')}"
    ins_code = item.get("insCode")
    if ins_code and str(ins_code).strip() not in ("", "."):
        residue += f".{ins_code}"
    alt_id = item.get("altId")
    if alt_id and str(alt_id).strip() not in ("", "."):
        residue += f" {alt_id}"
    return residue


def concerned_items(entries):
    """Flatten NAVAL per-residue entries to the geometry terms worth a look.

    Keeps terms whose NAVAL tier is not Preferred or whose probability group
    is not Common, and orders them most worrying first: Of Concern before
    Allowed, then ascending ProSco (a low probability score is worse).
    """
    tier_order = {"Of Concern": 0, "Allowed": 1}
    concerned = []
    for item in entries or []:
        for detail in item.get("details", []) or []:
            tier = detail.get("naval_tier")
            pgroup = detail.get("pGroup")
            if tier not in CONCERN_TIERS and pgroup not in CONCERN_PGROUPS:
                continue
            prosco = detail.get("prosco")
            concerned.append({
                "residue": format_residue(item),
                "name": format_geometry_name(detail.get("name")),
                "value": detail.get("value"),
                "pGroup": pgroup if pgroup is not None else "Unique",
                "prosco": prosco if isinstance(prosco, (int, float)) else 0.0,
                "naval_tier": tier or "",
            })
    concerned.sort(key=lambda x: (tier_order.get(x["naval_tier"], 2), x["prosco"]))
    return concerned


def naval_tier_counts(naval, kind):
    """{tier: (count, percentage)} for 'lengths' or 'angles'; {} if absent."""
    if not naval:
        return {}
    key = "navalLengthsStats" if kind == "lengths" else "navalAnglesStats"
    counts = {}
    for row in naval.get(key, []) or []:
        tier = row.get("navalTier")
        if tier:
            counts[tier] = (row.get("exclusiveCount"), row.get("exclusivePercentage"))
    return counts
