import xml.etree.ElementTree as ET
from .urls import redo_cif, redo_mtz, rcsb_mmcif, rcsb_pdb
from .utils import download, i2run


def test_7beq(cif7beq, mtz7beq):  # electron diffr.
    args = ["validate_protein"]
    args += ["--XYZIN_1", cif7beq]
    args += ["--F_SIGF_1", f"fullPath={mtz7beq}", "columnLabels=/*/*/[FP,SIGFP]"]
    with i2run(args) as job:
        xml = ET.parse(job / "program.xml")
        check_program_xml(xml)


def test_8ola_8olf():
    with download(redo_cif("8ola")) as cif8ola, download(redo_mtz("8ola")) as mtz8ola:
        with download(redo_cif("8olf")) as cif8olf, download(redo_mtz("8olf")) as mtz8olf:
            args = ["validate_protein"]
            args += ["--XYZIN_1", cif8ola]
            args += ["--XYZIN_2", cif8olf]
            args += ["--F_SIGF_1", f"fullPath={mtz8ola}", "columnLabels=/*/*/[FP,SIGFP]"]
            args += ["--F_SIGF_2", f"fullPath={mtz8olf}", "columnLabels=/*/*/[FP,SIGFP]"]
            args += ["--TWO_DATASETS", "True"]
            with i2run(args) as job:
                xml = ET.parse(job / "program.xml")
                check_program_xml(xml)


def test_8rk1_cif():  # neutron diffr., deuterium fraction
    with download(rcsb_mmcif("8rk1")) as cif8rk1:
        args = ["validate_protein"]
        args += ["--XYZIN_1", cif8rk1]
        with i2run(args) as job:
            xml = ET.parse(job / "program.xml")
            check_program_xml(xml)


def test_8rk1_pdb():  # neutron diffr., deuterium fraction
    with download(rcsb_pdb("8rk1")) as pdb8rk1:
        args = ["validate_protein"]
        args += ["--XYZIN_1", pdb8rk1]
        with i2run(args) as job:
            xml = ET.parse(job / "program.xml")
            check_program_xml(xml)


def check_program_xml(xml, expected_chain_count=1):
    for tag in ["Chain_count", "Iris", "Panel_svg", "Molprobity", "B_factors", "Ramachandran"]:
        elements = xml.findall(f".//{tag}")
        assert len(elements) > 0, f"XML tag missing in program.xml: {tag}"
        if tag == "Chain_count":
            assert elements[0].text == str(expected_chain_count), f"Expected {expected_chain_count} chain, got {elements[0].text}"


# The following test raises Segmentation fault and hangs forever.

# def test_8rk1_cif_sffrom_8c4y(cif8rk1):  # neutron diffr., deuterium fraction
#     with download(redo_mtz("8c4y")) as mtz8c4y:
#         args = ["validate_protein"]
#         args += ["--XYZIN_1", cif8rk1]
#         args += ["--F_SIGF_1", f"fullPath={mtz8c4y}", "columnLabels=/*/*/[FP,SIGFP]"]
#         with i2run(args) as job:
#             xml = ET.parse(job / "program.xml")
#             check_program_xml(xml)
