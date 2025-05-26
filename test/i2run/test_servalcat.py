import xml.etree.ElementTree as ET
from .utils import hasLongLigandName, i2run, testData


def test_8xfm(cif8xfm, mtz8xfm):
    args = ["servalcat_pipe"]
    args += ["--XYZIN", cif8xfm]
    args += ["--DATA_METHOD", "xtal"]
    args += ["--HKLIN", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FREE]"]
    args += ["--NCYCLES", "2"]
    args += ["--F_SIGF_OR_I_SIGI", "F_SIGF"]
    args += ["--VALIDATE_IRIS", "False"]
    args += ["--VALIDATE_BAVERAGE", "False"]
    args += ["--VALIDATE_RAMACHANDRAN", "False"]
    args += ["--VALIDATE_MOLPROBITY", "False"]
    with i2run(args) as job:
        assert hasLongLigandName(job / "CIFFILE.pdb")
        xml = ET.parse(job / "program.xml")
        rworks = [float(e.text) for e in xml.findall(".//summary/Rwork")]
        rfrees = [float(e.text) for e in xml.findall(".//summary/Rfree")]
        assert len(rworks) == len(rfrees) == 3
        assert rworks[-1] < rworks[0]
        assert rworks[-1] < 0.18
        assert rfrees[-1] < 0.25


def test_8ola(cif8ola, mtz8ola):
    args = ["servalcat_pipe"]
    args += ["--XYZIN", cif8ola]
    args += ["--DATA_METHOD", "xtal"]
    args += ["--HKLIN", f"fullPath={mtz8ola}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz8ola}", "columnLabels=/*/*/[FREE]"]
    args += ["--DICT_LIST", f"fullPath={testData('8OLA_A1H69.cif')}"]
    args += ["--NCYCLES", "2"]
    args += ["--F_SIGF_OR_I_SIGI", "F_SIGF"]
    args += ["--B_REFINEMENT_MODE", "aniso"]
    args += ["--H_OUT", "True"]
    args += ["--RANDOMIZEUSE", "True"]
    args += ["--RANDOMIZE", "0.05"]
    args += ["--VALIDATE_IRIS", "True"]
    args += ["--VALIDATE_BAVERAGE", "False"]
    args += ["--VALIDATE_RAMACHANDRAN", "True"]
    args += ["--VALIDATE_MOLPROBITY", "True"]
    args += ["--RUN_ADP_ANALYSIS", "True"]
    with i2run(args) as job:
        assert hasLongLigandName(job / "job_1" / "refined.mmcif")  # or job / 1_servalcat_test_0_ciffile_servalcat_pipe.cif ?
        xml = ET.parse(job / "program.xml")
        rworks = [float(e.text) for e in xml.findall(".//summary/Rwork")]
        rfrees = [float(e.text) for e in xml.findall(".//summary/Rfree")]
        assert len(rworks) == len(rfrees) == 3
        assert rworks[-1] < rworks[0]
        assert rworks[-1] < 0.15
        assert rfrees[-1] < 0.20
