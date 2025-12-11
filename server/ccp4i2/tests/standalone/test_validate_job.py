#!/usr/bin/env python3
"""
Test script for validate_job function using CPluginScript architecture.
"""
import os
import sys
from pathlib import Path

# Set up environment using relative path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'ccp4i2.settings')
os.environ.setdefault('CCP4I2_ROOT', str(PROJECT_ROOT))
os.environ.setdefault('CCP4_LOG_LEVEL', 'INFO')

# Add server directory to path
sys.path.insert(0, str(PROJECT_ROOT / 'server'))

import django
django.setup()

from ccp4i2.db.models import Job, Project
from ccp4i2.lib.utils.jobs.validate import validate_job

def main():
    # Find the test project
    proj = Project.objects.filter(name='set_param_test_2760').first()
    if not proj:
        print("ERROR: Project 'set_param_test_2760' not found")
        return 1

    # Get the job
    job = Job.objects.filter(project=proj, number=1).first()
    if not job:
        print(f"ERROR: Job 1 not found in project {proj.name}")
        return 1

    print(f"Testing validate_job on:")
    print(f"  Project: {proj.name}")
    print(f"  Job #: {job.number}")
    print(f"  Job UUID: {job.uuid}")
    print(f"  Task: {job.task_name}")
    print(f"  Status: {job.status}")
    print()

    # Test validate_job
    print("Calling validate_job()...")
    result = validate_job(job)

    if result.success:
        print("✅ validate_job SUCCEEDED!")
        print()

        error_tree = result.data
        error_reports = error_tree.findall('.//errorReport')

        print(f"Found {len(error_reports)} validation reports:")
        print()

        for idx, error_report in enumerate(error_reports, 1):
            class_name = error_report.findtext('className', 'Unknown')
            code = error_report.findtext('code', 'Unknown')
            description = error_report.findtext('description', 'No description')
            severity = error_report.findtext('severity', 'Unknown')
            object_name = error_report.findtext('objectName', '')

            print(f"  Report {idx}:")
            print(f"    Class: {class_name}")
            print(f"    Code: {code}")
            print(f"    Severity: {severity}")
            if object_name:
                print(f"    Object: {object_name}")
            print(f"    Description: {description}")
            print()

        # Show XML (formatted)
        from xml.etree import ElementTree as ET
        xml_str = ET.tostring(error_tree, encoding='unicode')
        print("Full XML output:")
        print("-" * 80)
        print(xml_str)
        print("-" * 80)

        return 0
    else:
        print(f"❌ validate_job FAILED: {result.error}")
        if result.error_details:
            print(f"   Details: {result.error_details}")
        return 1

if __name__ == '__main__':
    sys.exit(main())
