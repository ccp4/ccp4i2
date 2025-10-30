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
    args += ["--RUN_ADP_ANALYSIS", "False"]
    args += ["--RUN_COORDADPDEV_ANALYSIS", "False"]
    with i2run(args) as job:
        assert hasLongLigandName(job / "CIFFILE.pdb")
        xml = ET.parse(job / "program.xml")
        check_program_xml_pipeline(xml, args)
        check_r_and_cc(xml, expectedLen=6, maxRwork=0.18, maxRfree=0.25)


# x-ray diffraction data
# input as unmerged data
# refinement against intensities
# add hydrogens
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
        args += ["--H_OUT", "True"]
        args += ["--VALIDATE_IRIS", "True"]
        args += ["--VALIDATE_BAVERAGE", "False"]
        args += ["--VALIDATE_RAMACHANDRAN", "True"]
        args += ["--VALIDATE_MOLPROBITY", "True"]
        args += ["--RUN_ADP_ANALYSIS", "True"]
        args += ["--RUN_COORDADPDEV_ANALYSIS", "True"]
        with i2run(args) as job:
            xml = ET.parse(job / "program.xml")
            check_program_xml_pipeline(xml, args)
            check_r_and_cc(
                xml,
                expectedLen=ncycle + 1,
                intensities=True,
                maxRwork=0.20,
                maxRfree=0.25,
                minCCwork=0.90,
                minCCfree=0.75,
            )


# electron diffraction data
# prosmart reference protein
def test_7beq_electron(cif7beq, mtz7beq):
    with download(redo_cif("7ber")) as cif7ber:
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
        args += ["--RUN_COORDADPDEV_ANALYSIS", "True"]
        with i2run(args) as job:
            assert os.path.isfile(os.path.join(job, "job_1", "RESTRAINTS.txt"))
            assert os.path.isfile(os.path.join(job, "job_1", "ProSMART_Results.html"))
            xml = ET.parse(job / "program.xml")
            check_program_xml_pipeline(xml, args)
            check_program_xml_validate(job)
            check_r_and_cc(
                xml,
                expectedLen=ncycle + 1,
                maxRwork=0.25,
                maxRfree=0.31,
                minCCwork=0.85,
                minCCfree=0.75,
            )


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
        args += ["--RUN_COORDADPDEV_ANALYSIS", "True"]
        with i2run(args) as job:
            assert os.path.isfile(os.path.join(job, "job_1", "CA.json"))
            assert os.path.isfile(os.path.join(job, "job_1", "CA_restraints.txt"))
            assert os.path.isfile(os.path.join(job, "metal_restraints.txt"))
            assert os.path.isfile(os.path.join(job, "metal_restraints.mmcif"))
            xml = ET.parse(job / "program.xml")
            check_program_xml_pipeline(xml, args)
            check_program_xml_validate(job, expected_chain_count=4)
            check_r_and_cc(
                xml,
                expectedLen=ncycle + 1,
                maxRwork=0.35,
                maxRfree=0.40,
                minCCwork=0.65,
                minCCfree=0.55,
            )


def check_program_xml_pipeline(xml, args):
    def has_two_following_args(arg1, arg2):
        return any(args[i] == arg1 and args[i + 1] == arg2 for i in range(len(args) - 1))

    root = xml.getroot()
    assert root.tag == "SERVALCAT"
    tags = ["SERVALCAT_FIRST"]
    if not has_two_following_args("--RUN_ADP_ANALYSIS", "False"):
        tags.append("ADP_ANALYSIS")
    if not has_two_following_args("--RUN_COORDADPDEV_ANALYSIS", "False"):
        tags.append("COORD_ADP_DEV")
    if not (
        has_two_following_args("--VALIDATE_IRIS", "False")
        and has_two_following_args("--VALIDATE_BAVERAGE", "False")
        and has_two_following_args("--VALIDATE_RAMACHANDRAN", "False")
        and has_two_following_args("--VALIDATE_MOLPROBITY", "False")
    ):
        tags.append("Validation")
    if has_two_following_args("--ADD_WATERS", "True"):
        tags.append("CootAddWaters")
        tags.append("SERVALCAT_WATERS")
    for tag in tags:
        element = root.find(tag)
        assert element is not None, f"Missing tag in program.xml: {tag}"
        assert len(element) or element.text is not None, f"Tag {tag} in program.xml has no content"


def check_program_xml_validate(job, expected_chain_count=1):
    xml = ET.parse(job / "job_3" / "program.xml")
    for xml_tag in ["Chain_count", "Iris", "Panel_svg", "Molprobity", "B_factors", "Ramachandran"]:
        xml_tag_findall = xml.findall(f".//{xml_tag}")
        assert len(xml_tag_findall) > 0, f"XML tag missing in program.xml: {xml_tag}"
        if xml_tag == "Chain_count":
            assert xml_tag_findall[0].text == str(expected_chain_count), f"Expected {expected_chain_count} chain, got {xml_tag_findall[0].text}"


def check_r_and_cc(
    xml: ET.ElementTree,
    expectedLen: int,
    intensities: bool = False,
    maxRwork: float = None,
    maxRfree: float = None,
    minCCwork: float = None,
    minCCfree: float = None,
):
    rLabel = "R1" if intensities else "R"
    ccLabel = "CCI" if intensities else "CCF"
    rworks = [float(e.text) for e in xml.findall(f".//data/summary/{rLabel}work")]
    rfrees = [float(e.text) for e in xml.findall(f".//data/summary/{rLabel}free")]
    ccworks = [float(e.text) for e in xml.findall(f".//data/summary/{ccLabel}workavg")]
    ccfrees = [float(e.text) for e in xml.findall(f".//data/summary/{ccLabel}freeavg")]
    assert len(rworks) == len(rfrees) == len(ccworks) == len(ccfrees) == expectedLen
    assert rworks[0] >= rworks[-1]
    assert ccworks[0] <= ccworks[-1]
    assert maxRwork is None or rworks[-1] <= maxRwork
    assert maxRfree is None or rfrees[-1] <= maxRfree
    assert minCCwork is None or ccworks[-1] >= minCCwork
    assert minCCfree is None or ccfrees[-1] >= minCCfree
