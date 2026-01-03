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
