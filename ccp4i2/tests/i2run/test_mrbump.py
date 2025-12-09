import gemmi
from .utils import demoData, i2run


def test_mrbump():
    args = ["mrbump_basic"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--FREERFLAG", demoData("gamma", "freeR.mtz")]
    args += ["--ASUIN", demoData("gamma", "gamma.asu.xml")]
    args += ["--MRMAX", "5"]
    args += ["--REDUNDANCYLEVEL", "95"]
    args += ["--PJOBS", "2"]
    args += ["--NCYC", "10"]
    with i2run(args) as job:
        gemmi.read_pdb(str(job / "output_mrbump_1.pdb"))
        gemmi.read_mtz_file(str(job / "output_mrbump_1.mtz"))
        path = job / "search_mrbump_1" / "results" / "results.txt"
        lines = path.read_text().splitlines()
        best_rfree = 1.0
        for i, line in enumerate(lines):
            if line.endswith("Rfree"):
                try:
                    rfree = float(lines[i + 1].split()[-1])
                    best_rfree = min(best_rfree, (rfree))
                except ValueError:
                    pass
        assert best_rfree < 0.35
