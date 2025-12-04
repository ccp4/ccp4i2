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
import os


def main():
    """
    Main entry point for i2run.

    Detects the backend environment and delegates to the appropriate runner.
    """
    # Ensure ccp4i2 root is in path
    script_dir = os.path.dirname(os.path.abspath(__file__))
    ccp4i2_root = os.path.dirname(os.path.dirname(script_dir))
    if ccp4i2_root not in sys.path:
        sys.path.insert(0, ccp4i2_root)

    # Detect backend
    try:
        from baselayer import DJANGO
        use_django = DJANGO()
    except ImportError:
        use_django = False

    if use_django:
        # Django mode: Use CCP4i2RunnerDjango
        try:
            from .CCP4i2RunnerDjango import CCP4i2RunnerDjango as Runner
        except ImportError:
            # Fall back to server path
            from server.ccp4x.i2run.CCP4i2RunnerDjango import CCP4i2RunnerDjango as Runner

        command_line = ' '.join(sys.argv[1:])
        runner = Runner(command_line=command_line)
        job_id, exit_code = runner.execute()
        sys.exit(exit_code)
    else:
        # Qt mode: Use legacy CCP4I2Runner
        try:
            from core.CCP4I2Runner import main as legacy_main
            legacy_main()
        except ImportError as e:
            print(f"Error: Could not import legacy runner: {e}")
            print("Make sure you're running in the correct environment.")
            sys.exit(1)


if __name__ == '__main__':
    main()
