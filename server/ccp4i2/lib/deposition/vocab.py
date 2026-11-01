"""Program names and the other vocabularies the deposition task maps onto
OneDep's controlled values. The values themselves come from
:mod:`enumerations` (pdb_extract's list, the PDBx dictionary); this module
holds what those sources lack: how people spell program names, which names
beyond the list OneDep is known to accept, and the consistency checks.
"""
import re

from . import enumerations

# Lower-case aliases -> the spelling in pdb_extract's _software.name list
# (the guaranteed floor of what OneDep accepts).
_LISTED_ALIASES = {
    "aimless": "Aimless", "dials": "DIALS", "xia2": "xia2", "xds": "XDS", "xscale": "XSCALE",
    "autoproc": "AutoPROC", "mosflm": "MOSFLM", "imosflm": "iMOSFLM", "scala": "SCALA",
    "hkl-2000": "HKL-2000", "hkl2000": "HKL-2000", "hkl-3000": "HKL-3000", "hkl3000": "HKL-3000",
    "molrep": "MOLREP", "phaser": "PHASER", "phasertng": "PHASER", "phaser-tng": "PHASER",
    "amore": "AMoRE", "shelxd": "SHELXD", "shelxe": "SHELXE", "shelxl": "SHELXL",
    "shelxcd": "SHELXCD", "parrot": "PARROT", "buccaneer": "BUCCANEER", "coot": "Coot",
    "refmac": "REFMAC", "refmac5": "REFMAC", "refmacat": "REFMAC", "phenix": "PHENIX",
    "phenix.refine": "PHENIX", "buster": "BUSTER", "buster-tnt": "BUSTER-TNT", "cns": "CNS",
    "pdb_extract": "PDB_EXTRACT", "arp/warp": "ARP", "arp_warp": "ARP", "arpwarp": "ARP",
    "mrbump": "MrBUMP", "ample": "AMPLE", "arcimboldo": "Arcimboldo", "morda": "MoRDa",
    "rosetta": "Rosetta", "mr-rosetta": "MR-Rosetta", "balbes": "BALBES", "sharp": "SHARP",
    "autosharp": "autoSHARP", "crystfel": "CrystFEL", "dm": "DM", "acorn": "ACORN",
}

# Names the live DepUI accepts though pdb_extract's list lacks them, each
# with how it was verified.
LIVE_DEPUI_SEEN = {
    # DepUI software-name menu, deposition D_1292162556, 2026-11-01 (offered as "pointless").
    "pointless": "pointless",
}

# Programs people use that neither source above lists. They are written
# (the archive records what was used), but validity() and the report warn
# that the DepUI may not offer them: pick the nearest in the form, or say so
# through OneDep's Communication tab.
UNVERIFIED = {
    "ctruncate": "CTRUNCATE", "truncate": "TRUNCATE", "servalcat": "Servalcat",
    "staraniso": "STARANISO", "modelcraft": "ModelCraft", "crank2": "Crank2", "crank": "Crank2",
    "shelxc": "SHELXC", "molprobity": "MolProbity", "pisa": "PISA", "ccp4": "CCP4",
    "alphafold": "AlphaFold", "alphafold2": "AlphaFold", "alphafold3": "AlphaFold",
    "rosettafold": "RoseTTAFold", "esmfold": "ESMFold", "swiss-model": "SWISS-MODEL",
    "modeller": "MODELLER", "rfdiffusion": "RFdiffusion", "proteinmpnn": "ProteinMPNN",
    "chai-1": "Chai-1", "boltz": "Boltz",
}

# Never written: OneDep has no use for the GUI that ran the programs.
SOFTWARE_NOT_WRITTEN = {"ccp4i2"}

LISTED, LIVE, UNKNOWN = "listed", "live", "unverified"


def software_name(name):
    """``(spelling, status)`` for a program name: status ``listed`` (in
    pdb_extract's list), ``live`` (accepted by the live DepUI), ``unverified``
    (written, with a warning); ``(None, None)`` for a name never written."""
    key = str(name or "").strip().lower()
    if not key or key in SOFTWARE_NOT_WRITTEN:
        return None, None
    listed = enumerations.canonical("_software.name", name)
    if listed:
        return listed, LISTED
    if key in _LISTED_ALIASES and enumerations.is_allowed("_software.name", _LISTED_ALIASES[key]):
        return _LISTED_ALIASES[key], LISTED
    if key in LIVE_DEPUI_SEEN:
        return LIVE_DEPUI_SEEN[key], LIVE
    if key in UNVERIFIED:
        return UNVERIFIED[key], UNKNOWN
    return str(name).strip(), UNKNOWN


def canonical_software_name(name):
    """The spelling to write, or None for a name that is never written."""
    return software_name(name)[0]


def unverified_software_warning(spelling):
    return ("program %s is in neither pdb_extract's software list nor the names the live DepUI "
            "is known to offer: it is written, but if the DepUI rejects it choose the nearest "
            "name there or tell the annotators through the Communication tab" % spelling)


def _values(item):
    return enumerations.allowed(item) or ()


def classifications():
    return _values("_struct_keywords.pdbx_keywords")


def assembly_evidence():
    return _values("_pdbx_struct_assembly_auth_evidence.experimental_support")


def starting_model_types():
    return _values("_pdbx_initial_refinement_model.type")


def starting_model_sources():
    return _values("_pdbx_initial_refinement_model.source_name")


_PHASING_FROM_SHORT = {"MR": "MOLECULAR REPLACEMENT", "AB INITIO": "AB INITIO PHASING"}


def phasing_method(value):
    """The _refine.pdbx_method_to_determine_struct value for a phasing
    method given as the form or an older job had it ("MR"), or None."""
    if not value:
        return None
    text = str(value).strip().upper()
    text = _PHASING_FROM_SHORT.get(text, text)
    return enumerations.canonical("_refine.pdbx_method_to_determine_struct", text)


# Home-source anode lines (A): Cu, Mo, Ga (liquid metal jet), Cr, Ag K-alpha.
HOME_SOURCE_LINES = (1.5418, 0.7107, 1.3414, 2.2909, 0.5609)
HOME_SOURCES = ("ROTATING ANODE", "SEALED TUBE")


def wavelength_consistency(source, wavelength):
    """A warning when the wavelength does not look like the source: an anode
    line at a synchrotron, or a non-anode wavelength from a home source."""
    try:
        wl = float(wavelength)
    except (TypeError, ValueError):
        return None
    near_line = any(abs(wl - line) < 0.003 for line in HOME_SOURCE_LINES)
    src = str(source or "").upper()
    if src == "SYNCHROTRON" and near_line:
        return ("wavelength %.5f A is a home-source anode line (Cu/Mo/Ga K-alpha) but the "
                "source is SYNCHROTRON: check the source or the wavelength" % wl)
    if src in HOME_SOURCES and not near_line:
        return ("wavelength %.5f A is not a Cu, Mo or Ga K-alpha line but the source is %s: "
                "check the source or the wavelength" % (wl, src))
    return None


AUTHOR_FORMAT = re.compile(
    r"^[^\s,][^,]*, (?:[A-Z][a-z]?\.(?:-?[A-Z][a-z]?\.)*)(?: (?:Jr|Sr|II|III|IV)\.?)?$")


def author_problem(name):
    """Why ``name`` is not in the PDB's "Surname, I.I." form, or None."""
    text = str(name or "").strip()
    if not text or AUTHOR_FORMAT.match(text):
        return None
    if "," not in text:
        return ("'%s' has no initials: write it as Surname, I.I. (for example Noble, M.E.M.)"
                % text)
    return ("'%s' is not in the form Surname, I.I.: initials with full stops and no spaces "
            "between them (for example Noble, M.E.M.)" % text)
