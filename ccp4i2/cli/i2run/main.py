#!/usr/bin/env python3
"""
i2run - CCP4i2 command-line job runner

Main entry point for the i2run command-line tool.
Supports both Django and legacy Qt backends via the baselayer module.

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
