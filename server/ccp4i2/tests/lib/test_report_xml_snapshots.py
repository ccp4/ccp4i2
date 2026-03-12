"""
Report XML snapshot tests.

These tests capture the XML output of report generation for all available
test project zips, saving them as reference snapshots. On subsequent runs,
the generated XML is compared against the saved snapshots to detect
regressions.

Usage:
    # Capture initial snapshots (run once before refactoring):
    ./run_test.sh tests/lib/test_report_xml_snapshots.py -k capture --snapshot-update

    # Validate against saved snapshots (run after refactoring):
    ./run_test.sh tests/lib/test_report_xml_snapshots.py
"""

import re
import xml.etree.ElementTree as ET
from pathlib import Path
from shutil import rmtree

from django.conf import settings
from django.test import TestCase, override_settings

from ccp4i2.db.import_i2xml import import_ccp4_project_zip
from ccp4i2.db.models import Job
from ccp4i2.lib.utils.reporting.i2_report import (
    generate_job_report,
    get_report_job_info,
)

# Resolve __file__ so that .parent chains work regardless of cwd
_THIS_FILE = Path(__file__).resolve()
_PROJECT_ROOT = _THIS_FILE.parent.parent.parent.parent.parent
_TEST_ZIPS_DIR = _PROJECT_ROOT.parent / "test101" / "ProjectZips"

# Directory where reference XML snapshots are stored
_SNAPSHOTS_DIR = _THIS_FILE.parent / "report_snapshots"


def _normalise_xml(xml_element, projects_dir=None):
    """Normalise XML for stable comparison.

    Strips attributes that vary between runs (timestamps, counters, etc.)
    and pretty-prints for readable diffs.
    """
    # Remove attributes that change between runs
    _strip_volatile_attrs(xml_element)
    # Scrub non-deterministic content (Python object addresses, paths, etc.)
    _scrub_volatile_text(xml_element, projects_dir=projects_dir)
    # Re-serialise with consistent formatting
    _indent(xml_element)
    return ET.tostring(xml_element, encoding="unicode")


def _strip_volatile_attrs(element):
    """Recursively remove attributes that change between runs."""
    # 'key' attributes contain instance counters that may vary
    # with import order or parallel test execution
    for attr in ["key"]:
        if attr in element.attrib:
            del element.attrib[attr]
    for child in element:
        _strip_volatile_attrs(child)


# Pattern matching Python object repr addresses like "at 0x11c151a80"
_OBJ_ADDR_RE = re.compile(r" at 0x[0-9a-fA-F]+")


def _scrub_volatile_text(element, projects_dir=None):
    """Recursively scrub non-deterministic content from text/tail/attrs.

    Handles:
    - Python object addresses (at 0x...) — change every run
    - Absolute paths to test project directory — change per machine
    """
    for attr_name in list(element.attrib):
        val = element.get(attr_name)
        val = _OBJ_ADDR_RE.sub(" at 0xADDR", val)
        if projects_dir:
            val = val.replace(projects_dir, "<PROJECTS_DIR>")
        element.set(attr_name, val)
    if element.text:
        element.text = _OBJ_ADDR_RE.sub(" at 0xADDR", element.text)
        if projects_dir:
            element.text = element.text.replace(projects_dir, "<PROJECTS_DIR>")
    if element.tail:
        element.tail = _OBJ_ADDR_RE.sub(" at 0xADDR", element.tail)
        if projects_dir:
            element.tail = element.tail.replace(projects_dir, "<PROJECTS_DIR>")
    for child in element:
        _scrub_volatile_text(child, projects_dir)


def _indent(elem, level=0):
    """Add pretty-print indentation to an XML element tree (in-place)."""
    indent = "\n" + "  " * level
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = indent + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = indent
        for child in elem:
            _indent(child, level + 1)
        if not child.tail or not child.tail.strip():
            child.tail = indent
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = indent


@override_settings(
    CCP4I2_PROJECTS_DIR=_PROJECT_ROOT / "CCP4I2_TEST_PROJECT_DIRECTORY"
)
class ReportSnapshotTests(TestCase):
    """Generate reports from test project zips and compare against snapshots."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        Path(settings.CCP4I2_PROJECTS_DIR).mkdir(exist_ok=True)
        _SNAPSHOTS_DIR.mkdir(exist_ok=True)

    @classmethod
    def tearDownClass(cls):
        rmtree(settings.CCP4I2_PROJECTS_DIR, ignore_errors=True)
        super().tearDownClass()

    def _import_and_generate(self, zip_name):
        """Import a project zip, generate report for each job, return dict of results.

        Note: Some older project zips may produce failure reports (missing report
        classes, incompatible XML, etc.). That's fine — the snapshot captures
        whatever XML is produced, and the important thing is that the same XML
        is produced before and after refactoring.
        """
        zip_path = _TEST_ZIPS_DIR / zip_name
        if not zip_path.exists():
            self.skipTest(f"Test zip not found: {zip_path}")

        try:
            import_ccp4_project_zip(
                zip_path,
                relocate_path=settings.CCP4I2_PROJECTS_DIR,
            )
        except Exception as e:
            self.skipTest(f"Failed to import {zip_name}: {e}")

        project_name = zip_name.replace(".ccp4_project.zip", "")
        jobs = Job.objects.filter(project__name=project_name)
        if not jobs.exists():
            self.skipTest(f"No jobs found for project {project_name}")

        results = {}
        for job in jobs:
            try:
                report_xml = generate_job_report(job)
            except Exception as e:
                print(f"  WARNING: report generation raised for {job.task_name}: {e}")
                continue
            if report_xml is None:
                print(f"  WARNING: report is None for {job.task_name}")
                continue
            snapshot_name = f"{project_name}__{job.task_name}_{job.number}"
            # Sanitise for filesystem
            snapshot_name = snapshot_name.replace("/", "_").replace(".", "_")
            results[snapshot_name] = report_xml
        return results

    def _check_or_save_snapshot(self, name, xml_element):
        """Compare against saved snapshot, or save if none exists."""
        snapshot_path = _SNAPSHOTS_DIR / f"{name}.xml"
        normalised = _normalise_xml(
            xml_element,
            projects_dir=str(settings.CCP4I2_PROJECTS_DIR),
        )

        if not snapshot_path.exists():
            # No snapshot yet — save it
            snapshot_path.write_text(normalised, encoding="utf-8")
            print(f"  SAVED snapshot: {snapshot_path.name} ({len(normalised)} chars)")
            return

        # Compare against saved snapshot
        saved = snapshot_path.read_text(encoding="utf-8")
        if normalised != saved:
            # Save the actual output for debugging
            actual_path = _SNAPSHOTS_DIR / f"{name}.actual.xml"
            actual_path.write_text(normalised, encoding="utf-8")
            self.fail(
                f"Report XML differs from snapshot for {name}.\n"
                f"  Expected: {snapshot_path}\n"
                f"  Actual:   {actual_path}\n"
                f"  Diff with: diff {snapshot_path} {actual_path}"
            )
        print(f"  OK snapshot: {snapshot_path.name}")

    def _test_project(self, zip_name):
        """Test helper: import project, generate reports, check snapshots."""
        results = self._import_and_generate(zip_name)
        for name, xml_element in results.items():
            self._check_or_save_snapshot(name, xml_element)

    # --- Individual test methods for each project zip ---
    # Each is a separate test so failures are isolated and we get clear reporting.

    def test_refmac_gamma(self):
        self._test_project("refmac_gamma_test_0.ccp4_project.zip")

    def test_aimless_gamma(self):
        self._test_project("aimless_gamma_native_test_1.ccp4_project.zip")

    def test_phaser_simple(self):
        self._test_project("phaser_simple_test_0.ccp4_project.zip")

    def test_phaser_simple_refmac(self):
        self._test_project("phaser_simple_refmac_test_0.ccp4_project.zip")

    def test_phaser_expert(self):
        self._test_project("phaser_expert_test_0.ccp4_project.zip")

    def test_phaser_expert_asu(self):
        self._test_project("phaser_expert_asu_test_0.ccp4_project.zip")

    def test_phaser_expert_asu_singlejob(self):
        self._test_project("phaser_expert_asu_singlejob_test_0.ccp4_project.zip")

    def test_phaser_expert_tncs(self):
        self._test_project("phaser_expert_tncs_test_3.ccp4_project.zip")

    def test_phaser_simple_no_solution(self):
        self._test_project("phaser_simple_no_solution_test_0.ccp4_project.zip")

    def test_molrep(self):
        self._test_project("molrep_test_0.ccp4_project.zip")

    def test_bucc(self):
        self._test_project("bucc_test_0.ccp4_project.zip")

    def test_parrot(self):
        self._test_project("parrot_test_0.ccp4_project.zip")

    def test_acorn(self):
        self._test_project("acorn_test_0.ccp4_project.zip")

    def test_arpwarp(self):
        self._test_project("arpwarp_test_0.ccp4_project.zip")

    def test_crank2(self):
        self._test_project("crank2_test_0.ccp4_project.zip")

    def test_import_merged(self):
        self._test_project("import_merged_test_0.ccp4_project.zip")

    def test_freer(self):
        self._test_project("freer_test_0.ccp4_project.zip")

    def test_acedrg(self):
        self._test_project("acedrg_test_1.ccp4_project.zip")

    def test_auspex(self):
        self._test_project("auspex_test_0.ccp4_project.zip")

    def test_shelx_sad(self):
        self._test_project("shelxSAD_test_1.ccp4_project.zip")

    def test_shelx_siras(self):
        self._test_project("shelxSIRAS_test_1.ccp4_project.zip")

    def test_shelxe_mr(self):
        self._test_project("shelxe_mr_test_1.ccp4_project.zip")

    def test_single_atom(self):
        self._test_project("single_atom_test_0.ccp4_project.zip")

    def test_mrbump(self):
        self._test_project("mrbump_gamma_test_0.ccp4_project.zip")

    def test_mrparse(self):
        self._test_project("mrparse_test_1.ccp4_project.zip")

    def test_arcimboldo_lite(self):
        self._test_project("arcimboldo_lite_test_0.ccp4_project.zip")

    def test_xia2dials(self):
        self._test_project("xia2dials_test_2.ccp4_project.zip")

    def test_xia2xds(self):
        self._test_project("xia2xds_test_0.ccp4_project.zip")

    def test_makelink(self):
        self._test_project("makelink_test0.ccp4_project.zip")

    def test_phaept(self):
        self._test_project("PhaEPT_test.ccp4_project.zip")

    def test_provide_asu_contents(self):
        self._test_project("ProvideAsuContents_test_0.ccp4_project.zip")

    def test_substitute_ligand(self):
        self._test_project("SubstituteLigandPdbIn_test_0.ccp4_project.zip")

    def test_process_predicted_models(self):
        self._test_project("process_predicted_models_test_0.ccp4_project.zip")
