"""The reflections mmCIF for deposition, written with gemmi's MtzToCif."""
import os

# CCP4i2's mini-MTZ column labels, as makeHklin0 prefixes them, to the names
# gemmi's deposition spec recognises.
COLUMN_RENAME = {
    "FREERFLAG_FREER": "FREE",
    "F_SIGF_Iplus": "I(+)", "F_SIGF_SIGIplus": "SIGI(+)",
    "F_SIGF_Iminus": "I(-)", "F_SIGF_SIGIminus": "SIGI(-)",
    "F_SIGF_Fplus": "F(+)", "F_SIGF_SIGFplus": "SIGF(+)",
    "F_SIGF_Fminus": "F(-)", "F_SIGF_SIGFminus": "SIGF(-)",
    "F_SIGF_F": "FP", "F_SIGF_SIGF": "SIGFP",
    "F_SIGF_Fmean": "FP", "F_SIGF_SIGFmean": "SIGFP",
    "F_SIGF_I": "IMEAN", "F_SIGF_SIGI": "SIGIMEAN",
    "F_SIGF_Imean": "IMEAN", "F_SIGF_SIGImean": "SIGIMEAN",
    "FPHIOUT_F": "FWT", "FPHIOUT_PHI": "PHWT",
    "DIFFPHIOUT_F": "DELFWT", "DIFFPHIOUT_PHI": "PHDELWT",
}


def write_reflections_cif(mtz_path, cif_path, unmerged_mtz=None, log=print):
    """Write ``cif_path`` from a merged MTZ (and the scaled unmerged MTZ when
    given, which is what OneDep asks for). Returns the list of columns renamed."""
    import gemmi

    mtz = gemmi.read_mtz_file(os.path.normpath(str(mtz_path)))
    renamed = []
    for col in mtz.columns:
        if col.label in COLUMN_RENAME:
            old = col.label
            col.label = COLUMN_RENAME[old]
            renamed.append("%s -> %s" % (old, col.label))
    if renamed:
        log("Renamed columns: " + ", ".join(renamed))
    conv = gemmi.MtzToCif()
    conv.with_history = True
    if unmerged_mtz:
        unmerged = gemmi.read_mtz_file(os.path.normpath(str(unmerged_mtz)))
        text = conv.write_cif_to_string(mtz, unmerged)
    else:
        text = conv.write_cif_to_string(mtz)
    with open(str(cif_path), "w") as handle:
        handle.write(text)
    return renamed
