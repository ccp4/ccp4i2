from pathlib import Path
import gemmi
from .urls import rcsb_mmcif
from .utils import download, i2run


def test_rRNA():
    with download(rcsb_mmcif("1q93")) as mmcif_1q93:
        args = ["dnatco_pipe"]
        args += ["--XYZIN1", mmcif_1q93]
        args += ["--GENERATE_RESTRAINTS", "True"]
        with i2run(args) as job:
            check_output(job)


def test_2rRNA():
    with download(rcsb_mmcif("1q93")) as mmcif_1q93, download(rcsb_mmcif("1q96")) as mmcif_1q96:
        args = ["dnatco_pipe"]
        args += ["--XYZIN1", mmcif_1q93]
        args += ["--TOGGLE_XYZIN2", "True"]
        args += ["--XYZIN2", mmcif_1q96]
        args += ["--GENERATE_RESTRAINTS", "True"]
        args += ["--MAX_RMSD", "0.49"]
        args += ["--RESTRAINTS_SIGMA", "0.9"]
        with i2run(args) as job:
            check_output(job, nCiffiles=1)


def check_output(job: Path, nCiffiles=1):
    assert (job / "RESTRAINTS.txt").is_file(), "Restraints output file not found"
    for i in range(nCiffiles):
        ciffilePath = job / f"CIFOUT{i + 1}.pdb"
        assert ciffilePath.is_file(), "Extended CIF output file not found"
        gemmi.read_structure(str(ciffilePath), format=gemmi.CoorFormat.Mmcif)
        ciffile = gemmi.cif.read_file(str(ciffilePath))
        values = ['_ndb_struct_ntc_overall.confal_score',
                   '_ndb_struct_ntc_overall.confal_percentile',
                   '_ndb_struct_ntc_overall.num_classified',
                   '_ndb_struct_ntc_overall.num_unclassified',
                   '_ndb_struct_ntc_overall.num_unclassified_rmsd_close']
        for value in values:
            assert ciffile[0].find_value(value), f"Missing entry {value} in CIF file"
    assert (job / "RESTRAINTS.txt").is_file(), "Restraints output file not found"
