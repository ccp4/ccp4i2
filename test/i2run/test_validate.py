import xml.etree.ElementTree as ET
from .utils import i2run


def check_program_xml(xml, expected_chain_count=1):
    for xml_tag in ["Chain_count", "Iris", "Panel_svg", "Molprobity", "B_factors", "Ramachandran"]:
        xml_tag_findall = xml.findall(f".//{xml_tag}")
        assert len(xml_tag_findall) > 0, f"XML tag missing in program.xml: {xml_tag}"
        if xml_tag == "Chain_count":
            assert xml_tag_findall[0].text == str(expected_chain_count), f"Expected {expected_chain_count} chain, got {xml_tag_findall[0].text}"


def test_8xfm_nosf(cif8xfm):
    args = ["validate_protein"]
    args += ["--XYZIN_1", cif8xfm]
    with i2run(args) as job:
        xml = ET.parse(job / "program.xml")
        check_program_xml(xml)


def test_8xfm_sf(cif8xfm, mtz8xfm):
    args = ["validate_protein"]
    args += ["--XYZIN_1", cif8xfm]
    args += ["--F_SIGF_1", mtz8xfm]
    with i2run(args) as job:
        xml = ET.parse(job / "program.xml")
        check_program_xml(xml)


def test_8ola_sf(cif8ola, mtz8ola):
    args = ["validate_protein"]
    args += ["--XYZIN_1", cif8ola]
    args += ["--F_SIGF_1", mtz8ola]
    with i2run(args) as job:
        xml = ET.parse(job / "program.xml")
        check_program_xml(xml)


def test_8ola_8olf_nosf(cif8ola, cif8olf):
    args = ["validate_protein"]
    args += ["--XYZIN_1", cif8ola]
    args += ["--XYZIN_2", cif8olf]
    args += ["--TWO_DATASETS", "True"]
    with i2run(args) as job:
        xml = ET.parse(job / "program.xml")
        check_program_xml(xml)


def test_8ola_8olf_sf(cif8ola, cif8olf, mtz8ola, mtz8olf):
    args = ["validate_protein"]
    args += ["--XYZIN_1", cif8ola]
    args += ["--F_SIGF_1", mtz8ola]
    args += ["--XYZIN_2", cif8olf]
    args += ["--F_SIGF_2", mtz8olf]
    args += ["--TWO_DATASETS", "True"]
    with i2run(args) as job:
        xml = ET.parse(job / "program.xml")
        check_program_xml(xml)


def test_8rk1_cif_nosf(cif8rk1):  # neutron diffr., deuterium fraction
    args = ["validate_protein"]
    args += ["--XYZIN_1", cif8rk1]
    with i2run(args) as job:
        xml = ET.parse(job / "program.xml")
        check_program_xml(xml)


def test_8rk1_pdb_nosf(pdb8rk1):  # neutron diffr., deuterium fraction
    args = ["validate_protein"]
    args += ["--XYZIN_1", pdb8rk1]
    with i2run(args) as job:
        xml = ET.parse(job / "program.xml")
        check_program_xml(xml)


def test_8rk1_cif_sffrom_8c4y(cif8rk1, mtz8c4y):  # neutron diffr., deuterium fraction
    args = ["validate_protein"]
    args += ["--XYZIN_1", cif8rk1]
    args += ["--F_SIGF_1", mtz8c4y]
    with i2run(args) as job:
        xml = ET.parse(job / "program.xml")
        check_program_xml(xml)


def test_7beq_cif_nosf(cif7beq):  # electron diffr.
    args = ["validate_protein"]
    args += ["--XYZIN_1", cif7beq]
    with i2run(args) as job:
        xml = ET.parse(job / "program.xml")
        check_program_xml(xml)


def test_7beq_cif_nosf(cif7beq, mtz7beq):  # electron diffr.
    args = ["validate_protein"]
    args += ["--XYZIN_1", cif7beq]
    args += ["--F_SIGF_1", mtz7beq]
    with i2run(args) as job:
        xml = ET.parse(job / "program.xml")
        check_program_xml(xml)