import xml.etree.ElementTree as ET
import pytest
from .urls import redo_cif, redo_mtz, rcsb_mmcif, rcsb_pdb
from .utils import download, i2run


# ===== Tests WITHOUT MolProbity requirement =====
# These tests verify core validation functionality (Iris, Ramachandran, B-factors)
# without requiring the optional rotarama_data from Top8000 database

def test_7beq_basic(cif7beq, mtz7beq):
    """Test basic validation without MolProbity (electron diffr.)"""
    args = ["validate_protein"]
    args += ["--XYZIN_1", cif7beq]
    args += ["--F_SIGF_1", f"fullPath={mtz7beq}", "columnLabels=/*/*/[FP,SIGFP]"]
    with i2run(args) as job:
        xml = ET.parse(job / "program.xml")
        check_program_xml_basic(xml)


def test_8ola_8olf_basic():
    """Test basic validation without MolProbity (two datasets)"""
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
                check_program_xml_basic(xml)


def test_8rk1_cif_basic():
    """Test basic validation without MolProbity (neutron diffr., mmcif)"""
    with download(rcsb_mmcif("8rk1")) as cif8rk1:
        args = ["validate_protein"]
        args += ["--XYZIN_1", cif8rk1]
        with i2run(args) as job:
            xml = ET.parse(job / "program.xml")
            check_program_xml_basic(xml)


def test_8rk1_pdb_basic():
    """Test basic validation without MolProbity (neutron diffr., pdb)"""
    with download(rcsb_pdb("8rk1")) as pdb8rk1:
        args = ["validate_protein"]
        args += ["--XYZIN_1", pdb8rk1]
        with i2run(args) as job:
            xml = ET.parse(job / "program.xml")
            check_program_xml_basic(xml)


# ===== Tests WITH MolProbity requirement =====
# These tests explicitly require MolProbity functionality and will fail if
# rotarama_data is not available (requires chem_data symlink - see CLAUDE.md)

def test_8ola_8olf_with_molprobity():
    """Test validation WITH MolProbity requirement (two datasets)"""
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
                check_program_xml_with_molprobity(xml)


def test_8rk1_cif_with_molprobity():
    """Test validation WITH MolProbity requirement (neutron diffr., mmcif)"""
    with download(rcsb_mmcif("8rk1")) as cif8rk1:
        args = ["validate_protein"]
        args += ["--XYZIN_1", cif8rk1]
        with i2run(args) as job:
            xml = ET.parse(job / "program.xml")
            check_program_xml_with_molprobity(xml)


def test_8rk1_pdb_with_molprobity():
    """Test validation WITH MolProbity requirement (neutron diffr., pdb)"""
    with download(rcsb_pdb("8rk1")) as pdb8rk1:
        args = ["validate_protein"]
        args += ["--XYZIN_1", pdb8rk1]
        with i2run(args) as job:
            xml = ET.parse(job / "program.xml")
            check_program_xml_with_molprobity(xml)


# ===== Helper functions =====

def check_program_xml_basic(xml, expected_chain_count=1):
    """Check validation XML output for required tags only (no MolProbity)"""
    required_tags = ["Chain_count", "Iris", "Panel_svg", "B_factors", "Ramachandran"]

    for tag in required_tags:
        elements = xml.findall(f".//{tag}")
        assert len(elements) > 0, f"XML tag missing in program.xml: {tag}"
        if tag == "Chain_count":
            assert elements[0].text == str(expected_chain_count), f"Expected {expected_chain_count} chain, got {elements[0].text}"


def check_program_xml_with_molprobity(xml, expected_chain_count=1):
    """Check validation XML output including MolProbity (requires rotarama_data)"""
    # Check basic tags first
    check_program_xml_basic(xml, expected_chain_count)

    # Now explicitly require MolProbity
    molprobity_elements = xml.findall(".//Molprobity")
    assert len(molprobity_elements) > 0, (
        "MolProbity tag missing in program.xml. "
        "This test requires rotarama_data from Top8000 database. "
        "Ensure chem_data symlink is created (see CLAUDE.md setup instructions)."
    )


# The following test raises Segmentation fault and hangs forever.

# def test_8rk1_cif_sffrom_8c4y(cif8rk1):  # neutron diffr., deuterium fraction
#     with download(redo_mtz("8c4y")) as mtz8c4y:
#         args = ["validate_protein"]
#         args += ["--XYZIN_1", cif8rk1]
#         args += ["--F_SIGF_1", f"fullPath={mtz8c4y}", "columnLabels=/*/*/[FP,SIGFP]"]
#         with i2run(args) as job:
#             xml = ET.parse(job / "program.xml")
#             check_program_xml(xml)
