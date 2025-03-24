from .urls import pdbe_fasta, redo_mtz
from .utils import download, i2run

import gemmi


def test_arpwarp():
    args = ["arp_warp_classic"]
    with download(redo_mtz("1o6a")) as mtz:
        with download(pdbe_fasta("1o6a")) as fasta:
            args += ["--AWA_FOBS", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]"]
            args += ["--AWA_FREE", f"fullPath={mtz}", "columnLabels=/*/*/[FREE]"]
            args += ["--AWA_PHINI", f"fullPath={mtz}", "columnLabels=/*/*/[PHIC_ALL,FOM]"]
            args += ["--AWA_PHREF", f"fullPath={mtz}", "columnLabels=/*/*/[PHIC_ALL,FOM]"]
            args += ["--AWA_SEQIN", f"seqFile={fasta}"]
            args += ["--AWA_BIG_CYCLES", "1"]
            args += ["--AWA_SMALL_CYCLES", "1"]
            with i2run(args) as directory:
                for name in ["XYZDUM", "XYZOUT"]:
                    gemmi.read_pdb(str(directory / f"{name}.pdb"))
                for name in ["DIFFPHIOUT", "FPHIOUT"]:
                    gemmi.read_mtz_file(str(directory / f"{name}.mtz"))
