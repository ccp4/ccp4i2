# Copyright (C) 2025 University of York
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
i2run - CCP4i2 command-line job runner

This module provides the i2run command-line tool for executing CCP4i2 jobs.
It supports both Django and legacy Qt backends.

Usage:
    python -m cli.i2run <task> [options]

Or via shell script:
    i2run <task> [options]
"""

from .main import main

__all__ = ['main']
