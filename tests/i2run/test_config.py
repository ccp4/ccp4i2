"""
Test configuration and utilities for i2run tests.

Provides:
- Centralized test project directory management
- Ordered, unique project naming with timestamps
- Cleanup utilities and session-end messaging
"""

import os
import secrets
from datetime import datetime
from pathlib import Path


def get_test_projects_dir() -> Path:
    """Get the directory for test projects.

    Priority:
    1. CCP4I2_TEST_PROJECTS env var (explicit override)
    2. ~/.cache/ccp4i2-tests/ (default - persistent, easy to find)

    Returns:
        Path to the test projects directory (created if it doesn't exist)
    """
    explicit = os.environ.get('CCP4I2_TEST_PROJECTS')
    if explicit:
        test_dir = Path(explicit)
    else:
        test_dir = Path.home() / '.cache' / 'ccp4i2-tests'

    # Ensure directory exists
    test_dir.mkdir(parents=True, exist_ok=True)
    return test_dir


def make_test_project_name(test_name: str) -> str:
    """Create an ordered, unique test project name.

    Format: YYYYMMDD_HHMMSS_XXXX_<test_name>

    - Timestamp ensures chronological ordering via simple `ls` sort
    - 4-char random suffix prevents collisions in parallel runs
    - Test name provides human-readable identification

    Args:
        test_name: Base name for the test (e.g., "servalcat_test_8xfm")

    Returns:
        Unique project name like "20251204_143052_a1b2_servalcat_test_8xfm"

    Examples:
        >>> make_test_project_name("phaser_test_gamma")
        '20251204_143052_a1b2_phaser_test_gamma'

        # Easy to identify runs from today:
        # ls ~/.cache/ccp4i2-tests/20251204_*

        # Easy to see chronological order:
        # ls ~/.cache/ccp4i2-tests/ | sort
    """
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    suffix = secrets.token_hex(2)  # 4 hex chars

    # Clean up test name (remove special chars that could cause issues)
    clean_name = test_name.replace('[', '_').replace(']', '_').replace('/', '_').replace('::', '_')

    return f"{timestamp}_{suffix}_{clean_name}"


def get_cleanup_command() -> str:
    """Get the shell command to clean up test cache.

    Returns:
        Shell command string for cleaning up test projects
    """
    test_dir = get_test_projects_dir()
    return f"rm -rf {test_dir}/*"


def get_cleanup_message() -> str:
    """Get a formatted message about cleaning up test cache.

    Returns:
        Multi-line string with cleanup instructions
    """
    test_dir = get_test_projects_dir()

    # Count projects and calculate size
    try:
        projects = list(test_dir.iterdir())
        num_projects = len([p for p in projects if p.is_dir()])

        # Calculate total size
        total_size = 0
        for item in test_dir.rglob('*'):
            if item.is_file():
                try:
                    total_size += item.stat().st_size
                except (OSError, IOError):
                    pass

        # Format size
        if total_size < 1024:
            size_str = f"{total_size} B"
        elif total_size < 1024 * 1024:
            size_str = f"{total_size / 1024:.1f} KB"
        elif total_size < 1024 * 1024 * 1024:
            size_str = f"{total_size / (1024 * 1024):.1f} MB"
        else:
            size_str = f"{total_size / (1024 * 1024 * 1024):.1f} GB"

        stats = f"({num_projects} project dirs, {size_str})"
    except Exception:
        stats = ""

    return f"""
================================================================================
Test Cache Cleanup
================================================================================
Test projects are stored in: {test_dir}
{stats}

To clean up all test projects:
  rm -rf {test_dir}/*

To clean up projects from a specific date (e.g., today 2024-12-04):
  rm -rf {test_dir}/20241204_*

To see failed tests (preserved for debugging):
  ls -la {test_dir}/
================================================================================
"""
