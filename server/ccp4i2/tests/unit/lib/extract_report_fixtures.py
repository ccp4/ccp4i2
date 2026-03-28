#!/usr/bin/env python
# Copyright (C) 2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
"""
Extract report test fixtures from project zips.

This script imports each test project zip, then for every job extracts:
  1. program.xml  — the raw program output XML
  2. job_info.json — the jobInfo dict used by report classes
  3. report.xml    — the normalised expected report output

These fixtures allow report snapshot tests to run without the 661MB of
project zips or any database — just plain files.

Usage:
    cd server
    source ../../ccp4-20251105/bin/ccp4.setup-sh
    ccp4-python -m pytest ccp4i2/tests/lib/extract_report_fixtures.py -s -k extract

Or via Django test runner:
    ccp4-python manage.py test ccp4i2.tests.lib.extract_report_fixtures -v2
"""

import json
import os
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
from ccp4i2.core.tasks import get_report_class

# Resolve paths
_THIS_FILE = Path(__file__).resolve()
_PROJECT_ROOT = _THIS_FILE.parent.parent.parent.parent.parent
_TEST_ZIPS_DIR = _PROJECT_ROOT.parent / "test101" / "ProjectZips"

# Output directory for extracted fixtures
_FIXTURES_DIR = _THIS_FILE.parent / "report_fixtures"


# --- Normalisation (same as test_report_xml_snapshots.py) ---

_OBJ_ADDR_RE = re.compile(r" at 0x[0-9a-fA-F]+")


def _normalise_xml(xml_element, projects_dir=None):
    _strip_volatile_attrs(xml_element)
    _scrub_volatile_text(xml_element, projects_dir=projects_dir)
    _indent(xml_element)
    return ET.tostring(xml_element, encoding="unicode")


def _strip_volatile_attrs(element):
    for attr in ["key"]:
        if attr in element.attrib:
            del element.attrib[attr]
    for child in element:
        _strip_volatile_attrs(child)


def _scrub_volatile_text(element, projects_dir=None):
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


def _sanitise_name(name):
    """Sanitise a string for use as a filesystem directory name."""
    return name.replace("/", "_").replace(".", "_")


def _scrub_job_info(job_info, projects_dir):
    """Replace machine-specific paths in jobInfo with placeholders."""
    scrubbed = {}
    for key, value in job_info.items():
        if isinstance(value, str) and projects_dir and projects_dir in value:
            scrubbed[key] = value.replace(projects_dir, "<PROJECTS_DIR>")
        elif isinstance(value, list):
            scrubbed[key] = _scrub_list(value, projects_dir)
        elif isinstance(value, dict):
            scrubbed[key] = _scrub_job_info(value, projects_dir)
        else:
            scrubbed[key] = value
    return scrubbed


def _scrub_list(lst, projects_dir):
    """Recursively scrub paths in a list."""
    result = []
    for item in lst:
        if isinstance(item, str) and projects_dir and projects_dir in item:
            result.append(item.replace(projects_dir, "<PROJECTS_DIR>"))
        elif isinstance(item, dict):
            result.append(_scrub_job_info(item, projects_dir))
        elif isinstance(item, list):
            result.append(_scrub_list(item, projects_dir))
        elif isinstance(item, tuple):
            result.append(_scrub_list(list(item), projects_dir))
        else:
            result.append(item)
    return result


@override_settings(
    CCP4I2_PROJECTS_DIR=_PROJECT_ROOT / "CCP4I2_TEST_PROJECT_DIRECTORY"
)
class ExtractFixtures(TestCase):
    """Extract report fixtures from all test project zips."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        Path(settings.CCP4I2_PROJECTS_DIR).mkdir(exist_ok=True)
        _FIXTURES_DIR.mkdir(exist_ok=True)

    @classmethod
    def tearDownClass(cls):
        rmtree(settings.CCP4I2_PROJECTS_DIR, ignore_errors=True)
        super().tearDownClass()

    def test_extract_all_fixtures(self):
        """Extract fixtures from all available test project zips."""
        if not _TEST_ZIPS_DIR.exists():
            self.skipTest(f"Test zips directory not found: {_TEST_ZIPS_DIR}")

        zip_files = sorted(_TEST_ZIPS_DIR.glob("*.ccp4_project.zip"))
        if not zip_files:
            self.skipTest("No test project zips found")

        total_jobs = 0
        total_extracted = 0
        tasks_seen = set()

        for zip_path in zip_files:
            project_name = zip_path.name.replace(".ccp4_project.zip", "")
            print(f"\n=== {project_name} ===")

            try:
                import_ccp4_project_zip(
                    zip_path,
                    relocate_path=settings.CCP4I2_PROJECTS_DIR,
                )
            except Exception as e:
                print(f"  SKIP (import failed): {e}")
                continue

            jobs = Job.objects.filter(project__name=project_name)
            if not jobs.exists():
                print("  SKIP (no jobs)")
                continue

            for job in jobs:
                total_jobs += 1
                task_name = job.task_name

                # Create fixture name
                fixture_name = f"{project_name}__{_sanitise_name(task_name)}_{job.number}"
                fixture_name = _sanitise_name(fixture_name)
                fixture_dir = _FIXTURES_DIR / fixture_name

                # Check if report class exists
                report_class = get_report_class(task_name)
                if report_class is None:
                    print(f"  SKIP {task_name} #{job.number} (no report class)")
                    continue

                # Generate report
                try:
                    report_xml = generate_job_report(job)
                except Exception as e:
                    print(f"  SKIP {task_name} #{job.number} (report failed): {e}")
                    continue

                if report_xml is None:
                    print(f"  SKIP {task_name} #{job.number} (report is None)")
                    continue

                # Get jobInfo
                try:
                    job_info = get_report_job_info(job.uuid)
                except Exception as e:
                    print(f"  SKIP {task_name} #{job.number} (jobInfo failed): {e}")
                    continue

                # Find program.xml
                job_dir = Path(job.directory)
                program_xml_path = None
                for xml_name in ["program.xml", "XMLOUT.xml", "i2.xml"]:
                    candidate = job_dir / xml_name
                    if candidate.exists():
                        program_xml_path = candidate
                        break

                # Create fixture directory
                fixture_dir.mkdir(exist_ok=True)

                # 1. Save program.xml (if it exists)
                if program_xml_path is not None:
                    program_xml_content = program_xml_path.read_text(encoding="utf-8", errors="replace")
                    (fixture_dir / "program.xml").write_text(program_xml_content, encoding="utf-8")
                else:
                    # Some jobs have no program XML (e.g. wrapper-only tasks)
                    (fixture_dir / "program.xml").write_text("<empty/>", encoding="utf-8")

                # 2. Save job_info.json (scrubbed of machine-specific paths)
                scrubbed_info = _scrub_job_info(job_info, str(settings.CCP4I2_PROJECTS_DIR))
                # Add metadata for the test harness
                scrubbed_info["_fixture_meta"] = {
                    "task_name": task_name,
                    "report_class": report_class.__name__,
                    "report_module": report_class.__module__,
                    "job_number": job.number,
                    "job_status": Job.Status(job.status).label,
                    "project_name": project_name,
                    "has_program_xml": program_xml_path is not None,
                }
                (fixture_dir / "job_info.json").write_text(
                    json.dumps(scrubbed_info, indent=2, default=str),
                    encoding="utf-8",
                )

                # 3. Save normalised expected report XML
                normalised = _normalise_xml(
                    report_xml,
                    projects_dir=str(settings.CCP4I2_PROJECTS_DIR),
                )
                (fixture_dir / "expected_report.xml").write_text(normalised, encoding="utf-8")

                tasks_seen.add(task_name)
                total_extracted += 1
                print(f"  OK {task_name} #{job.number}")

        print(f"\n{'='*60}")
        print(f"Extracted {total_extracted}/{total_jobs} job fixtures")
        print(f"Unique task types: {len(tasks_seen)}")
        print(f"Tasks: {', '.join(sorted(tasks_seen))}")
        print(f"Output: {_FIXTURES_DIR}")
