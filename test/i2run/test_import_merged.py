import gemmi
from .urls import pdbe_sfcif
from .utils import demoData, download, i2run


def test_2ceu_cif():
    with download(pdbe_sfcif("2ceu")) as cif:
        args = ["import_merged"]
        args += ["--HKLIN", cif]
        args += ["--SPACEGROUP", "I 2 2 2"]
        with i2run(args) as job:
            check_output(job, cif)


def test_gamma_mtz():
    ianom = demoData("gamma", "merged_intensities_Xe.mtz")
    freer = demoData("gamma", "freeR.mtz")
    args = ["import_merged"]
    args += ["--HKLIN", ianom]
    args += ["--FREERFLAG", freer]
    with i2run(args) as job:
        check_output(job, freer)


def check_output(job, freerin):
    gemmi.read_mtz_file(str(job / "OBSOUT.mtz"))
    input_flags = freer_flag_dict(freerin)
    output_flags = freer_flag_dict(str(job / "FREEOUT.mtz"))
    for hkl, free in input_flags.items():
        assert free == output_flags[hkl], f"Free-R flag mismatch at {hkl}"


def freer_flag_dict(hklin):
    if hklin.endswith(".mtz"):
        mtz = gemmi.read_mtz_file(hklin)
    else:
        doc = gemmi.cif.read(hklin)
        rblock = gemmi.as_refln_blocks(doc)[0]
        mtz = gemmi.CifToMtz().convert_block_to_mtz(rblock)
    mtz.ensure_asu()
    hcol = mtz.column_with_label("H")
    kcol = mtz.column_with_label("K")
    lcol = mtz.column_with_label("L")
    freecol = mtz.rfree_column()
    return {
        (h, k, l): min(free, 1)
        for h, k, l, free in zip(hcol, kcol, lcol, freecol)
    }
