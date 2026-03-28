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
"""Conftest for report fixture tests."""


def pytest_addoption(parser):
    """Add --rebaseline option for updating expected fixture output."""
    parser.addoption(
        "--rebaseline",
        action="store_true",
        default=False,
        help="Update expected_report.xml with actual output for failing fixtures.",
    )
