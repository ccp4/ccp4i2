"""Data reduction statistics for the deposited model.

Aimless (through aimless_pipe) writes ``DataStatistics.cif`` with the
``_reflns`` and ``_reflns_shell`` categories OneDep wants in the model file;
those are copied into the model block category by category, each written
once. Older scaling jobs have the statistics only in the Aimless XML: that is
shown in the report but not merged into the model (OneDep's form takes them
by hand), which the task says in its log.
"""
from lxml import etree

from .mmcif_edit import category_rows, replace_loop, replace_pairs

# Categories taken from DataStatistics.cif. Everything else there (_entry,
# _cell, _symmetry) the model already has, and a second copy is what used to
# crash MolProbity and EDS on the validation server.
STATISTICS_CATEGORIES = ("_reflns", "_reflns_shell")
_SHELL_ORDINAL = "pdbx_ordinal"


# Aimless writes the "_all" statistics items (by policy: no weak-reflection
# cut was made) but OneDep's form reads the "_obs" ones -- unique counts,
# Rmerge and I/sigma per shell came up blank on the Ahir3252 deposition, and
# number_unique_all was not read at all. Both are written: the _obs item is
# copied from its _all twin when absent.
OBS_FROM_ALL = {
    "_reflns_shell": (("number_unique_obs", "number_unique_all"),
                      ("Rmerge_I_obs", "Rmerge_I_all"),
                      ("meanI_over_sigI_obs", "pdbx_netI_over_sigmaI_all"),
                      ("percent_possible_obs", "percent_possible_all"),
                      ("number_measured_obs", "number_measured_all")),
    "_reflns": (("number_obs", "number_all"),
                ("pdbx_Rmerge_I_obs", "pdbx_Rmerge_I_all"),
                ("pdbx_netI_over_sigmaI", "pdbx_netI_over_sigmaI_all"),
                ("percent_possible_obs", "percent_possible_all")),
}


def with_obs_items(tags, rows, category):
    """``(tags, rows)`` with the _obs items OneDep reads added from their
    _all twins where the file lacks them."""
    tags = list(tags)
    rows = [list(r) for r in rows]
    for obs, all_ in OBS_FROM_ALL.get(category, ()):
        if obs in tags or all_ not in tags:
            continue
        src = tags.index(all_)
        tags.append(obs)
        for row in rows:
            row.append(row[src])
    return tags, rows


def outer_shell(tags, rows):
    """The one row with the highest resolution (smallest d_res_high),
    renumbered 1: the shell OneDep's form expects, rather than every Aimless
    bin from 80 A outwards (each low-resolution bin raised a warning)."""
    if "d_res_high" not in tags or len(rows) <= 1:
        return rows
    col = tags.index("d_res_high")

    def d_high(row):
        try:
            return float(row[col])
        except ValueError:
            return float("inf")

    best = list(min(rows, key=d_high))
    if "pdbx_ordinal" in tags:
        best[tags.index("pdbx_ordinal")] = "1"
    return [best]


def merge_statistics(model_block, stats_path, all_shells=False):
    """Copy the statistics categories from ``stats_path`` into ``model_block``,
    replacing any the model already has, with the _obs items OneDep reads and
    (unless ``all_shells``) only the outermost resolution shell. Returns the
    categories written."""
    import gemmi

    stats_block = gemmi.cif.read(str(stats_path))[0]
    written = []
    for category in STATISTICS_CATEGORIES:
        tags, rows = category_rows(stats_block, category)
        if not rows:
            continue
        tags, rows = with_obs_items(tags, rows, category)
        if category == "_reflns_shell" and not all_shells:
            rows = outer_shell(tags, rows)
        # Values come across as CIF tokens already (quoted as needed); write
        # them back verbatim rather than re-quoting.
        if category == "_reflns" and len(rows) == 1:
            _replace_pairs_raw(model_block, category, list(zip(tags, rows[0])))
        else:
            _replace_loop_raw(model_block, category, tags, rows)
        written.append(category)
    return written


def _replace_pairs_raw(block, category, items):
    replace_pairs(block, category, [(name, None) for name, _ in items])
    for name, token in items:
        block.set_pair("%s.%s" % (category, name), token)


def _replace_loop_raw(block, category, tags, rows):
    replace_loop(block, category, tags, [["?"] * len(tags)])
    loop = block.find_loop("%s.%s" % (category, tags[0])).get_loop()
    loop.set_all_values([[row[j] for row in rows] for j in range(len(tags))])


STATS_TABLE = [
    ("Low resolution limit", "_reflns.d_resolution_low", "_reflns_shell.d_res_low"),
    ("High resolution limit", "_reflns.d_resolution_high", "_reflns_shell.d_res_high"),
    ("Rmerge (all I+ and I-)", "_reflns.pdbx_Rmerge_I_all", "_reflns_shell.Rmerge_I_all"),
    ("Rmeas (all I+ & I-)", "_reflns.pdbx_Rrim_I_all", "_reflns_shell.pdbx_Rrim_I_all"),
    ("Rpim (all I+ & I-)", "_reflns.pdbx_Rpim_I_all", "_reflns_shell.pdbx_Rpim_I_all"),
    ("Number of observations", "_reflns.pdbx_number_measured_all", "_reflns_shell.number_measured_all"),
    ("Number unique", "_reflns.number_obs", "_reflns_shell.number_unique_all"),
    ("Mean((I)/sd(I))", "_reflns.pdbx_netI_over_sigmaI", "_reflns_shell.pdbx_netI_over_sigmaI_all"),
    ("Half-set correlation CC(1/2)", "_reflns.pdbx_CC_half", "_reflns_shell.pdbx_CC_half"),
    ("Completeness %", "_reflns.percent_possible_obs", "_reflns_shell.percent_possible_all"),
    ("Multiplicity", "_reflns.pdbx_redundancy", "_reflns_shell.pdbx_redundancy"),
    ("Anomalous completeness %", "_reflns.pdbx_percent_possible_anomalous",
     "_reflns_shell.pdbx_percent_possible_anomalous"),
    ("Anomalous multiplicity", "_reflns.pdbx_redundancy_anomalous", "_reflns_shell.pdbx_redundancy_anomalous"),
    ("Anomalous CC(1/2)", "_reflns.pdbx_CC_half_anomalous", "_reflns_shell.pdbx_CC_half_anomalous"),
]


def cifstats_to_xml(cifstats_path):
    """A ``<CifStats>`` element (overall / inner / outer shell per statistic)
    for the report, from ``DataStatistics.cif``."""
    import gemmi

    block = gemmi.cif.read(str(cifstats_path))[0]
    shell_loop = None
    for item in block:
        if item.loop is not None and any("_reflns_shell" in t for t in item.loop.tags):
            shell_loop = item.loop
            break
    root = etree.Element("CifStats")
    for label, reflns_tag, shell_tag in STATS_TABLE:
        overall = block.find_value(reflns_tag)
        if not overall or overall == "?":
            continue
        stat = etree.SubElement(root, "Stat")
        stat.set("label", label)
        stat.set("overall", overall)
        if shell_loop is not None:
            shell_tags = list(shell_loop.tags)
            if shell_tag in shell_tags:
                ci = shell_tags.index(shell_tag)
                stat.set("inner", shell_loop[0, ci])
                stat.set("outer", shell_loop[shell_loop.length() - 1, ci])
    wilson = block.find_value("_reflns.B_iso_Wilson_estimate")
    if wilson and wilson != "?":
        stat = etree.SubElement(root, "Stat")
        stat.set("label", "Wilson B estimate")
        stat.set("overall", wilson)
    return root
