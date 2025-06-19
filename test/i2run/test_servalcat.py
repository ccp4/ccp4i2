import os
import xml.etree.ElementTree as ET
from .urls import redo_cif, redo_mtz, pdbe_mmcif
from .utils import download, hasLongLigandName, i2run, demoData


# x-ray diffraction data
# monomer with 5-letter code
# add waters
def test_8xfm(cif8xfm, mtz8xfm):
    args = ["servalcat_pipe"]
    args += ["--XYZIN", cif8xfm]
    args += ["--DATA_METHOD", "xtal"]
    args += ["--HKLIN", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FREE]"]
    args += ["--NCYCLES", "2"]
    args += ["--ADD_WATERS", "True"]
    args += ["--NCYCLES_AFTER_ADD_WATERS", "2"]
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
        assert len(rworks) == len(rfrees) == 6
        assert rworks[-1] < rworks[0]
        assert rworks[-1] < 0.18
        assert rfrees[-1] < 0.25


# x-ray diffraction data
# input as unmerged data
# refinement against intensities
def test_1gyu_unmerged():
    with download(pdbe_mmcif("1gyu")) as cif1gyu:
        ncycle = 2
        args = ["servalcat_pipe"]
        args += ["--XYZIN", cif1gyu]
        args += ["--DATA_METHOD", "xtal"]
        args += ["--MERGED_OR_UNMERGED", "unmerged"]
        args += ["--HKLIN_UNMERGED", f"fullPath={demoData(os.path.join('gamma', 'HKLOUT_unmerged.mtz'))}"]
        args += ["--FREERFLAG", f"fullPath={demoData(os.path.join('gamma', 'freeR.mtz'))}", "columnLabels=/*/*/[FREER]"]
        args += ["--NCYCLES", str(ncycle)]
        args += ["--RES_CUSTOM", "True"]
        args += ["--RES_CUSTOM", "True"]
        args += ["--RES_MIN", "1.9"]
        args += ["--B_REFINEMENT_MODE", "aniso"]
        args += ["--H_OUT", "True"]
        args += ["--RANDOMIZEUSE", "True"]
        args += ["--RANDOMIZE", "0.05"]
        args += ["--BFACSETUSE", "True"]
        args += ["--VALIDATE_IRIS", "True"]
        args += ["--VALIDATE_BAVERAGE", "False"]
        args += ["--VALIDATE_RAMACHANDRAN", "True"]
        args += ["--VALIDATE_MOLPROBITY", "True"]
        args += ["--RUN_ADP_ANALYSIS", "True"]
        with i2run(args) as job:
            xml = ET.parse(job / "program.xml")
            r1works = [float(e.text) for e in xml.findall(".//data/summary/R1work")]
            r1frees = [float(e.text) for e in xml.findall(".//data/summary/R1free")]
            cciworkavgs = [float(e.text) for e in xml.findall(".//data/summary/CCIworkavg")]
            ccifreeavgs = [float(e.text) for e in xml.findall(".//data/summary/CCIfreeavg")]
            assert len(r1works) == len(r1frees) == len(cciworkavgs) == len(ccifreeavgs) == ncycle + 1
            assert r1works[-1] < 0.20
            assert r1frees[-1] < 0.25
            assert r1works[0] > r1works[-1]
            assert r1frees[0] > r1frees[-1]
            assert cciworkavgs[-1] > 0.90
            assert ccifreeavgs[-1] > 0.80
            assert cciworkavgs[0] < cciworkavgs[-1]
            assert ccifreeavgs[0] < ccifreeavgs[-1]


# electron diffraction data
# prosmart reference protein
def test_7beq_electron():
    with download(redo_cif("7beq")) as cif7beq, download(redo_mtz("7beq")) as mtz7beq, download(redo_cif("7ber")) as cif7ber:
        ncycle = 3
        args = ["servalcat_pipe"]
        args += ["--XYZIN", cif7beq]
        args += ["--DATA_METHOD", "xtal"]
        args += ["--HKLIN", f"fullPath={mtz7beq}", "columnLabels=/*/*/[FP,SIGFP]"]
        args += ["--FREERFLAG", f"fullPath={mtz7beq}", "columnLabels=/*/*/[FREE]"]
        args += ["--NCYCLES", str(ncycle)]
        args += ["--SCATTERING_FACTORS", "electron"]
        args += ["--USE_JELLY", "True"]
        args += ["--prosmartProtein.TOGGLE", "True"]
        args += ["--prosmartProtein.MODE", "SELECTED"]
        args += ["--prosmartProtein.REFERENCE_MODELS", f"fullPath={cif7ber}"]
        args += ["--prosmartProtein.CHAINLIST_1", "A"]
        args += ["--prosmartProtein.ALL_BEST", "ALL"]
        args += ["--prosmartProtein.SEQID", "75"]
        args += ["--prosmartProtein.SIDE_MAIN", "SIDE"]
        args += ["--VALIDATE_IRIS", "True"]
        args += ["--VALIDATE_BAVERAGE", "True"]
        args += ["--VALIDATE_RAMACHANDRAN", "True"]
        args += ["--VALIDATE_MOLPROBITY", "True"]
        args += ["--RUN_ADP_ANALYSIS", "True"]
        with i2run(args) as job:
            assert os.path.isfile(os.path.join(job, "job_1", "RESTRAINTS.txt"))
            assert os.path.isfile(os.path.join(job, "job_1", "ProSMART_Results.html"))
            xml = ET.parse(job / "program.xml")
            rworks = [float(e.text) for e in xml.findall(".//data/summary/Rwork")]
            rfrees = [float(e.text) for e in xml.findall(".//data/summary/Rfree")]
            ccfworkavgs = [float(e.text) for e in xml.findall(".//data/summary/CCFworkavg")]
            ccffreeavgs = [float(e.text) for e in xml.findall(".//data/summary/CCFfreeavg")]
            assert len(rworks) == len(rfrees) == len(ccfworkavgs) == len(ccffreeavgs) == ncycle + 1
            assert rworks[-1] < 0.25
            assert rworks[0] > rworks[-1]
            assert rfrees[-1] < 0.31
            assert rfrees[0] > rfrees[-1]
            assert ccfworkavgs[-1] > 0.85
            assert ccfworkavgs[0] < ccfworkavgs[-1]
            assert ccffreeavgs[-1] > 0.75
            assert ccffreeavgs[0] < ccffreeavgs[-1]
            xml_validation = ET.parse(job / "job_3" / "program.xml")
            check_program_xml_validate(xml_validation)


# neutron diffraction data
# metalCoord - Ca
# TODO: deuterium fraction
def test_7prg_neutron():
    with download(redo_cif("7prg")) as cif7prg, download(redo_mtz("7prg")) as mtz7prg:
        ncycle = 3
        args = ["servalcat_pipe"]
        args += ["--XYZIN", cif7prg]
        args += ["--DATA_METHOD", "xtal"]
        args += ["--HKLIN", f"fullPath={mtz7prg}", "columnLabels=/*/*/[FP,SIGFP]"]
        args += ["--FREERFLAG", f"fullPath={mtz7prg}", "columnLabels=/*/*/[FREE]"]
        args += ["--NCYCLES", str(ncycle)]
        args += ["--SCATTERING_FACTORS", "neutron"]
        args += ["--H_OUT", "True"]
        args += ["--USE_NCS", "True"]
        args += ["--RUN_METALCOORD", "True"]
        args += ["--GENERATE_OR_USE", "GENERATE"]
        args += ["--LIGAND_CODES_SELECTED", "CA"]
        args += ["--VALIDATE_IRIS", "True"]
        args += ["--VALIDATE_BAVERAGE", "True"]
        args += ["--VALIDATE_RAMACHANDRAN", "True"]
        args += ["--VALIDATE_MOLPROBITY", "True"]
        args += ["--RUN_ADP_ANALYSIS", "True"]
        with i2run(args) as job:
            assert os.path.isfile(os.path.join(job, "job_1", "CA.json"))
            assert os.path.isfile(os.path.join(job, "job_1", "CA_restraints.txt"))
            assert os.path.isfile(os.path.join(job, "metal_restraints.txt"))
            assert os.path.isfile(os.path.join(job, "metal_restraints.mmcif"))
            xml = ET.parse(job / "program.xml")
            rworks = [float(e.text) for e in xml.findall(".//data/summary/Rwork")]
            rfrees = [float(e.text) for e in xml.findall(".//data/summary/Rfree")]
            ccfworkavgs = [float(e.text) for e in xml.findall(".//data/summary/CCFworkavg")]
            ccffreeavgs = [float(e.text) for e in xml.findall(".//data/summary/CCFfreeavg")]
            assert len(rworks) == len(rfrees) == len(ccfworkavgs) == len(ccffreeavgs) == ncycle + 1
            assert rworks[-1] < 0.35
            assert rworks[0] > rworks[-1]
            assert rfrees[-1] < 0.40
            assert rfrees[0] > rfrees[-1]
            assert ccfworkavgs[-1] > 0.65
            assert ccfworkavgs[0] < ccfworkavgs[-1]
            assert ccffreeavgs[-1] > 0.55
            # assert ccffreeavgs[0] < ccffreeavgs[-1]
            xml_validation = ET.parse(job / "job_3" / "program.xml")
            check_program_xml_validate(xml_validation, expected_chain_count=4)


def check_program_xml_validate(xml, expected_chain_count=1):
    for xml_tag in ["Chain_count", "Iris", "Panel_svg", "Molprobity", "B_factors", "Ramachandran"]:
        xml_tag_findall = xml.findall(f".//{xml_tag}")
        assert len(xml_tag_findall) > 0, f"XML tag missing in program.xml: {xml_tag}"
        if xml_tag == "Chain_count":
            assert xml_tag_findall[0].text == str(expected_chain_count), f"Expected {expected_chain_count} chain, got {xml_tag_findall[0].text}"
