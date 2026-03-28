#!/usr/bin/env python3
# Copyright (C) 2025-2026 University of York
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

Main entry point for the i2run command-line tool.

Usage:
    python -m cli.i2run <task> [options]

Examples:
    python -m cli.i2run csymmatch --help
    python -m cli.i2run refmac --input-mtz data.mtz --input-pdb model.pdb
"""

import sys
from .CCP4i2RunnerDjango import CCP4i2RunnerDjango as Runner


def main():
    """
    Main entry point for i2run.

    Detects the backend environment and delegates to the appropriate runner.
    """
    command_line = ' '.join(sys.argv[1:])
    runner = Runner(command_line=command_line)
    job_id, exit_code = runner.execute()
    sys.exit(exit_code)


if __name__ == '__main__':
    main()
