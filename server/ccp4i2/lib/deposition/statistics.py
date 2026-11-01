"""Data reduction statistics for the deposited model.

Aimless (through aimless_pipe) writes ``DataStatistics.cif`` with the
``_reflns`` and ``_reflns_shell`` categories OneDep wants in the model file;
older jobs have them only in the Aimless XML. Either way the statistics are
copied into the model block category by category, each written once.
"""
from lxml import etree

from .mmcif_edit import category_rows, replace_loop, replace_pairs

# Categories taken from DataStatistics.cif. Everything else there (_entry,
# _cell, _symmetry) the model already has, and a second copy is what used to
# crash MolProbity and EDS on the validation server.
STATISTICS_CATEGORIES = ("_reflns", "_reflns_shell")
_SHELL_ORDINAL = "pdbx_ordinal"


def merge_statistics(model_block, stats_path):
    """Copy the statistics categories from ``stats_path`` into ``model_block``,
    replacing any the model already has. Returns the categories written."""
    import gemmi

    stats_block = gemmi.cif.read(str(stats_path))[0]
    written = []
    for category in STATISTICS_CATEGORIES:
        tags, rows = category_rows(stats_block, category)
        if not rows:
            continue
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
