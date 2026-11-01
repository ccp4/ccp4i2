"""Sanitising the mmCIF a refinement program wrote, before it is deposited.

BUSTER, Refmac/refmacat and servalcat each leave things in their mmCIF that
OneDep either rejects, misreads or silently carries into the archive. The
fixers here are selected from the evidence in the file (``_software``,
``_refine``, the program's private categories) and every change is recorded
as a :class:`Fix` ("found X, did Y") for the report, so the depositor can see
what was done to their file. Checks that find something the task cannot
repair come back as *problems*, and the task fails on them rather than
deposit a file it knows to be wrong.

What is verified against real output and what is defensive is said per
fixer. Refmac and servalcat cases come from the Ahir3252 deposition and DDU
jobs (2026-10/11); no BUSTER mmCIF was reachable when this was written, so
the BUSTER fixer is written from the program's documented conventions and
says so.
"""
import json
import os
from collections import namedtuple

import gemmi

from .mmcif_edit import category_rows, replace_pairs, rewrite_category, set_pairs_in_place, value_of

Fix = namedtuple("Fix", "program found did")

PLACEHOLDER_ENTRY = "XXXX"

# Refmac/refmacat's private categories: not in the PDBx dictionary, dropped.
REFMAC_PRIVATE = ("_ccp4_refine_ls", "_ccp4_refine_tls", "_ccp4_struct_ncs",
                  "_ccp4_form_factor", "_ccp4_reflns_twin")

# Header categories inherited from whatever PDB file the model started life
# as (a search model, an earlier deposition): never ours to keep.
INHERITED_HEADER = ("_pdbx_database_status", "_database_2", "_pdbx_database_PDB_obs_spr",
                    "_database_PDB_rev", "_database_PDB_rev_record", "_pdbx_audit_revision_history",
                    "_pdbx_audit_revision_details", "_pdbx_audit_revision_group",
                    "_pdbx_audit_revision_category", "_pdbx_audit_revision_item",
                    "_audit_author", "_citation", "_citation_author", "_struct_keywords",
                    "_pdbx_contact_author", "_struct", "_pdbx_database_related",
                    "_pdbx_SG_project", "_audit_conform", "_pdbx_deposit_group")


def detect_program(block):
    """``refmac``, ``servalcat``, ``buster`` or ``unknown``, from ``_software``
    and the program's own categories."""
    names = [row.str(0).lower() for row in block.find("_software.", ["name"])]
    if any(n.startswith("refmac") for n in names) or \
            any(block.find_mmcif_category(c + ".") for c in REFMAC_PRIVATE):
        return "refmac"
    if any("servalcat" in n for n in names):
        return "servalcat"
    if any("buster" in n for n in names) or \
            block.find_mmcif_category("_pdbx_refine_tls_group.") and \
            any("BUSTER" in (block.find_value("_refine.details") or "") for _ in [0]):
        return "buster"
    details = (block.find_value("_refine.details") or "").lower()
    if "buster" in details:
        return "buster"
    return "unknown"


def _erase_category(block, category):
    table = block.find_mmcif_category(category + ".")
    if table:
        table.erase()
        return True
    return False


def keep_model_block(doc):
    """Drop every data block but the first (Refmac appends the monomer
    library's dictionary blocks, 27 on Ahir3252). Returns a Fix or None."""
    names = [b.name for b in doc]
    if len(names) <= 1:
        return None
    for i in range(len(names) - 1, 0, -1):  # from the end, so indices stay valid
        del doc[i]
    return Fix("any", "%d data blocks (model + dictionary blocks)" % len(names),
               "kept only the model block %s" % names[0])


def strip_inherited_header(block):
    """Remove header categories inherited from an earlier file and reset the
    entry id. Returns the Fixes made."""
    fixes = []
    for category in INHERITED_HEADER:
        table = block.find_mmcif_category(category + ".")
        if not table:
            continue
        rows = len(table)
        sample = ""
        if category in ("_audit_author", "_citation_author"):
            names = [r.str(0) for r in block.find(category + ".", ["name"])][:3]
            sample = " (%s%s)" % (", ".join(names), "..." if rows > 3 else "")
        elif category == "_citation":
            sample = " (%s)" % gemmi.cif.as_string(block.find_value("_citation.title") or "?")[:60]
        table.erase()
        fixes.append(Fix("any", "inherited %s with %d row(s)%s" % (category, rows, sample),
                         "removed it; the task writes this entry's own"))
    entry = block.find_value("_entry.id")
    if entry is None or gemmi.cif.as_string(entry) != PLACEHOLDER_ENTRY:
        replace_pairs(block, "_entry", [("id", PLACEHOLDER_ENTRY)])
        fixes.append(Fix("any", "_entry.id %s" % (entry or "absent"),
                         "set the placeholder %s" % PLACEHOLDER_ENTRY))
    return fixes


def fix_refmac(block):
    """Refmac/refmacat: drop the private ``_ccp4_*`` categories and give the
    wavelength the id ``1`` (refmac names it after the dataset, ``FROMDIALS``).
    Verified on the Ahir3252 refmacat 5.8.0431 output."""
    fixes = []
    dropped = [c for c in REFMAC_PRIVATE if _erase_category(block, c)]
    if dropped:
        fixes.append(Fix("refmac", "private categories %s" % ", ".join(dropped),
                         "dropped them (not in the PDBx dictionary)"))
    # Refmac writes the wavelength as a one-row loop named after the dataset.
    # Rewriting the whole category keeps the wavelength beside the new id
    # (set_pair on a looped tag would drop the loop and the value with it).
    tags, rows = category_rows(block, "_diffrn_radiation_wavelength")
    if rows and "id" in tags:
        idx = tags.index("id")
        old = gemmi.cif.as_string(rows[0][idx])
        if old != "1":
            rows[0][idx] = "1"
            rewrite_category(block, "_diffrn_radiation_wavelength", tags, rows)
            rtags, rrows = category_rows(block, "_diffrn_radiation")
            if rrows and "wavelength_id" in rtags:
                rrows[0][rtags.index("wavelength_id")] = "1"
                rewrite_category(block, "_diffrn_radiation", rtags, rrows)
            fixes.append(Fix("refmac", "_diffrn_radiation_wavelength.id %s" % old,
                             "renamed it 1 (and _diffrn_radiation.wavelength_id with it)"))
    return fixes


def _has_value(block, tag):
    return value_of(block, tag) is not None


def servalcat_overall_statistics(job_dir):
    """Overall R-work and R-free from the last cycle of the ``refined_stats.json``
    servalcat writes in its job directory (searched recursively; the pipeline
    keeps the refinement in a sub-job).

    Looks for whichever of Rwork/Rfree (conventional) or R1work/R1free (the
    intensity-based R1 older servalcat reports) the job provides, preferring
    the conventional pair. Returns ``{"R_work", "R_free", "keys": (work, free),
    "source": path}`` or None.
    """
    # The job's own file first (a bare servalcat job); else the newest among
    # the pipeline's direct sub-jobs, which is the refinement that wrote the
    # final XYZOUT (servalcat_pipe runs its cycles in job_N_1, job_N_2, ...).
    candidates = []
    own = os.path.join(str(job_dir), "refined_stats.json")
    if os.path.exists(own):
        candidates.append((float("inf"), own))
    try:
        subdirs = [os.path.join(str(job_dir), d) for d in os.listdir(str(job_dir))]
    except OSError:
        subdirs = []
    for sub in subdirs:
        path = os.path.join(sub, "refined_stats.json")
        if os.path.isdir(sub) and os.path.exists(path):
            candidates.append((os.path.getmtime(path), path))
    for _mtime, path in sorted(candidates, reverse=True):
        try:
            with open(path) as handle:
                stats = json.load(handle)
            summary = stats[-1]["data"]["summary"]
        except (OSError, ValueError, IndexError, KeyError, TypeError):
            continue
        for work, free in (("Rwork", "Rfree"), ("R1work", "R1free"), ("R", None), ("R1", None)):
            if work in summary:
                return {"R_work": float(summary[work]),
                        "R_free": float(summary[free]) if free and free in summary else None,
                        "keys": (work, free), "source": path}
    return None


def fix_servalcat(block, refinement_job_dir=None):
    """servalcat up to 0.4.1xx writes ``_refine`` with no overall R-work/R-free
    (only ``_refine_ls_shell``), so wwPDB reports PDB-R NotAvailable and skips
    the EDS/RSRZ analysis. Verified on DDU job 60 (Ahir3252, servalcat 0.4.123).

    The fixer is conditional: when the file already carries the overall values
    (a forthcoming servalcat writes them itself) ``_refine`` is left as written
    and only verified. Otherwise :func:`_harvest_servalcat_r` fills them from
    the refinement job -- that branch is the part to delete once every
    servalcat in use writes them.

    Returns ``(fixes, problems)``: a problem when no R can be found.
    """
    if _has_value(block, "_refine.ls_R_factor_R_work"):
        fixes = [Fix("servalcat", "_refine carries overall R-work %s / R-free %s" % (
            gemmi.cif.as_string(block.find_value("_refine.ls_R_factor_R_work")),
            gemmi.cif.as_string(block.find_value("_refine.ls_R_factor_R_free") or "?")),
            "left as written")]
        return fixes, []
    return _harvest_servalcat_r(block, refinement_job_dir)


def _harvest_servalcat_r(block, refinement_job_dir):
    """The harvest branch of :func:`fix_servalcat` (removable once servalcat
    writes the overall R itself)."""
    fixes, problems = [], []
    stats = servalcat_overall_statistics(refinement_job_dir) if refinement_job_dir else None
    if not stats:
        problems.append(
            "servalcat wrote no overall R-work/R-free into _refine and none could be "
            "harvested from the refinement job%s; the file must not be deposited without "
            "them (wwPDB would report PDB-R NotAvailable and skip the EDS analysis)."
            % (" directory %s" % refinement_job_dir if refinement_job_dir else ""))
        return fixes, problems
    work_key, free_key = stats["keys"]
    counts = _sum_shell_counts(block)
    if work_key.startswith("R1"):
        # The intensity-based R1 (I/sigma > 2 cut) is not an R-work: it goes
        # where the dictionary keeps an R on intensities, and the file still
        # has no overall R-work, which the report says.
        items = [("ls_R_I_factor_obs", "%.4f" % stats["R_work"])]
        items += list(counts.items())
        set_pairs_in_place(block, "_refine", items)
        fixes.append(Fix("servalcat", "_refine without overall R-work/R-free; the job reports only "
                         "the intensity-based R1 (%s)" % "/".join(k for k in (work_key, free_key) if k),
                         "wrote R1 %.4f%s into _refine.ls_R_I_factor_obs from %s; no R-work/R-free "
                         "could be written -- expect wwPDB to report PDB-R NotAvailable" % (
                             stats["R_work"],
                             " (R1-free %.4f)" % stats["R_free"] if stats["R_free"] is not None else "",
                             os.path.basename(stats["source"]))))
    else:
        items = [("ls_R_factor_R_work", "%.4f" % stats["R_work"])]
        if stats["R_free"] is not None:
            items.append(("ls_R_factor_R_free", "%.4f" % stats["R_free"]))
        items += list(counts.items())
        set_pairs_in_place(block, "_refine", items)
        fixes.append(Fix("servalcat", "_refine without overall R-work/R-free",
                         "wrote R-work %.4f%s from %s of %s" % (
                             stats["R_work"],
                             " / R-free %.4f" % stats["R_free"] if stats["R_free"] is not None else "",
                             "/".join(k for k in (work_key, free_key) if k),
                             os.path.basename(stats["source"]))))
    if counts:
        fixes.append(Fix("servalcat", "_refine without reflection counts",
                         "summed _refine_ls_shell: " + ", ".join("%s=%s" % kv for kv in counts.items())))
    return fixes, problems


def _sum_shell_counts(block):
    totals = {}
    for tag, key in (("number_reflns_R_work", "ls_number_reflns_R_work"),
                     ("number_reflns_R_free", "ls_number_reflns_R_free"),
                     ("number_reflns_obs", "ls_number_reflns_obs")):
        column = block.find_loop("_refine_ls_shell." + tag)
        if not column:
            continue
        try:
            totals[key] = sum(int(float(v)) for v in column if v not in ("?", "."))
        except ValueError:
            continue
    if "ls_number_reflns_obs" not in totals and "ls_number_reflns_R_work" in totals:
        totals["ls_number_reflns_obs"] = totals["ls_number_reflns_R_work"] + \
            totals.get("ls_number_reflns_R_free", 0)
    return totals


def fix_buster(block):
    """BUSTER: written defensively, from the program's documented output,
    because no BUSTER mmCIF was reachable to verify against (the 2016 runs in
    the Ahir3252 archive are PDB format). ``_pdbx_refine_tls_group`` keeps its
    BUSTER selection syntax (OneDep accepts free text there);
    ``_refine.pdbx_refine_id`` and ``_refine_ls_restr.pdbx_refine_id`` are
    added when absent, which BUSTER has been seen to omit.
    """
    fixes = []
    if block.find_mmcif_category("_refine.") and \
            block.find_pair_item("_refine.pdbx_refine_id") is None:
        set_pairs_in_place(block, "_refine", [("pdbx_refine_id", "X-RAY DIFFRACTION")])
        fixes.append(Fix("buster", "_refine without pdbx_refine_id", "set X-RAY DIFFRACTION"))
    return fixes


def general_checks(block, keep_hydrogens=True):
    """Fixes and problems every deposition file gets: ``_exptl.method``,
    ``_refine.pdbx_refine_id``, a wavelength, occupancies, altlocs, hydrogens,
    anisotropic displacement parameters consistent with the model.

    Returns ``(fixes, problems, structure)``; the structure is re-read from
    the block after any atom edit.
    """
    fixes, problems = [], []
    if block.find_pair_item("_exptl.method") is None:
        set_pairs_in_place(block, "_exptl", [("entry_id", PLACEHOLDER_ENTRY),
                                             ("method", "X-RAY DIFFRACTION")])
        fixes.append(Fix("any", "no _exptl.method", "set X-RAY DIFFRACTION"))
    if block.find_mmcif_category("_refine.") and block.find_pair_item("_refine.pdbx_refine_id") is None:
        set_pairs_in_place(block, "_refine", [("pdbx_refine_id", "X-RAY DIFFRACTION")])
        fixes.append(Fix("any", "_refine without pdbx_refine_id", "set X-RAY DIFFRACTION"))
    # No wavelength here is not a problem yet: the task writes it from the
    # form or the reflection data (a PDB-format model never carries one), and
    # checks for it after that.

    st = gemmi.make_structure_from_block(block)
    st.setup_entities()
    bad_occ = 0
    altloc_singletons = 0
    n_h = 0
    anisou = 0
    for chain in st[0]:
        for res in chain:
            alts = {}
            for atom in res:
                if atom.occ < 0.0 or atom.occ > 1.0001:
                    bad_occ += 1
                if atom.is_hydrogen():
                    n_h += 1
                if atom.aniso.nonzero():
                    anisou += 1
                if atom.altloc != "\0":
                    alts.setdefault(atom.name, set()).add(atom.altloc)
            altloc_singletons += sum(1 for s in alts.values() if len(s) == 1)
    if bad_occ:
        problems.append("%d atom(s) have an occupancy outside [0, 1]" % bad_occ)
    if altloc_singletons:
        fixes.append(Fix("any", "%d atom(s) carry an altloc with no partner" % altloc_singletons,
                         "left as they are; check them in Coot before depositing"))
    if n_h and not keep_hydrogens:
        st.remove_hydrogens()
        groups = gemmi.MmcifOutputGroups(False)
        groups.atoms = True
        st.update_mmcif_block(block, groups)
        fixes.append(Fix("any", "%d hydrogen atoms (riding)" % n_h,
                         "removed them as asked; note that the R-factors and restraint counts in "
                         "_refine were computed with them present"))
        st = gemmi.make_structure_from_block(block)
        st.setup_entities()
    elif n_h:
        fixes.append(Fix("any", "%d hydrogen atoms (riding)" % n_h,
                         "kept: the deposited model is the refined and validated one"))
    anisotrop = block.find_mmcif_category("_atom_site_anisotrop.")
    has_tls = bool(block.find_mmcif_category("_pdbx_refine_tls."))
    if anisotrop and anisou == 0:
        anisotrop.erase()
        fixes.append(Fix("any", "_atom_site_anisotrop with no anisotropic atoms",
                         "dropped the empty table"))
    elif anisou and not anisotrop:
        problems.append("%d atoms are anisotropic but the file has no _atom_site_anisotrop" % anisou)
    elif anisou and not has_tls:
        fixes.append(Fix("any", "%d anisotropic atoms without TLS groups" % anisou,
                         "kept (individual anisotropic refinement, or TLS contributions "
                         "added without the groups being recorded)"))
    return fixes, problems, st


def sanitise(doc, refinement_job_dir=None, keep_hydrogens=True):
    """Run every fixer on ``doc`` in place.

    Returns ``(program, fixes, problems, structure)``.
    """
    fixes, problems = [], []
    kept = keep_model_block(doc)
    if kept:
        fixes.append(kept)
    block = doc[0]
    program = detect_program(block)
    fixes.extend(strip_inherited_header(block))
    if program == "refmac":
        fixes.extend(fix_refmac(block))
    elif program == "servalcat":
        more, bad = fix_servalcat(block, refinement_job_dir)
        fixes.extend(more)
        problems.extend(bad)
    elif program == "buster":
        fixes.extend(fix_buster(block))
    more, bad, st = general_checks(block, keep_hydrogens)
    fixes.extend(more)
    problems.extend(bad)
    return program, fixes, problems, st


def cell_matches_reflections(block, mtz_path, tolerance=0.002):
    """A problem string when the model's cell or space group disagrees with
    the reflections MTZ, else None."""
    try:
        mtz = gemmi.read_mtz_file(str(mtz_path))
    except Exception as err:
        return "could not read the reflections to compare the cell: %s" % err
    st = gemmi.make_structure_from_block(block)
    model, data = st.cell, mtz.cell
    for a, b, name in ((model.a, data.a, "a"), (model.b, data.b, "b"), (model.c, data.c, "c"),
                       (model.alpha, data.alpha, "alpha"), (model.beta, data.beta, "beta"),
                       (model.gamma, data.gamma, "gamma")):
        if b and abs(a - b) / b > tolerance:
            return "cell %s differs between model (%.3f) and reflections (%.3f)" % (name, a, b)
    model_sg = st.find_spacegroup()
    data_sg = mtz.spacegroup
    if model_sg and data_sg and model_sg.number != data_sg.number:
        return "space group differs between model (%s) and reflections (%s)" % (
            model_sg.hm, data_sg.hm)
    return None
