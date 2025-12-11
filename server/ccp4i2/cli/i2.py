#!/usr/bin/env python3
"""
i2 - Modern CLI for CCP4i2

A clean, resource-oriented command-line interface for CCP4i2 that wraps
Django management commands with a more intuitive syntax.

Usage:
    i2 projects [list]              List all projects
    i2 projects create <name>       Create a new project
    i2 projects show <id>           Show project details
    i2 projects tree <id>           Show project job tree

    i2 jobs <project> [list]        List jobs in a project
    i2 jobs create <project> <task> Create a new job
    i2 jobs run <job_id>            Run an existing job
    i2 jobs tree <job_id>           Show job file tree
    i2 jobs clone <job_id>          Clone a job

    i2 run <task> [options]         Create and run a task (i2run shortcut)

    i2 files <job_id> [list]        List files in a job
    i2 files cat <job_id> <name>    Display file contents

    i2 report <job_id>              Generate job report

Examples:
    i2 projects
    i2 projects create toxd
    i2 jobs toxd
    i2 run refmac --hklin data.mtz --xyzin model.pdb
    i2 report 42
"""

import os
import sys
import subprocess
from pathlib import Path


def get_manage_py():
    """Find the manage.py script."""
    # Try CCP4I2_ROOT first
    root = os.environ.get('CCP4I2_ROOT')
    if root:
        manage_py = Path(root) / 'server' / 'manage.py'
        if manage_py.exists():
            return str(manage_py)

    # Try relative to this script
    script_dir = Path(__file__).resolve().parent
    manage_py = script_dir.parent / 'server' / 'manage.py'
    if manage_py.exists():
        return str(manage_py)

    # Try current directory
    manage_py = Path.cwd() / 'server' / 'manage.py'
    if manage_py.exists():
        return str(manage_py)

    return None


def get_python():
    """Get the Python interpreter to use."""
    # Prefer ccp4-python if available
    import shutil
    ccp4_python = shutil.which('ccp4-python')
    if ccp4_python:
        return ccp4_python
    return sys.executable


def run_management_command(command, *args):
    """Run a Django management command."""
    manage_py = get_manage_py()
    if not manage_py:
        print("Error: Cannot find server/manage.py", file=sys.stderr)
        print("Set CCP4I2_ROOT or run from the ccp4i2 directory", file=sys.stderr)
        sys.exit(1)

    python = get_python()
    cmd = [python, manage_py, command] + list(args)

    # Pass through to subprocess
    result = subprocess.run(cmd)
    sys.exit(result.returncode)


def main():
    args = sys.argv[1:]

    if not args or args[0] in ['-h', '--help', 'help']:
        print(__doc__)
        sys.exit(0)

    resource = args[0]
    rest = args[1:]

    # ─────────────────────────────────────────────────────────────────────────
    # Projects
    # ─────────────────────────────────────────────────────────────────────────
    if resource == 'projects':
        if not rest or rest[0] == 'list':
            run_management_command('list_projects')
        elif rest[0] == 'create':
            if len(rest) < 2:
                print("Usage: i2 projects create <name>", file=sys.stderr)
                sys.exit(1)
            run_management_command('create_project', rest[1], *rest[2:])
        elif rest[0] == 'show':
            if len(rest) < 2:
                print("Usage: i2 projects show <project_id>", file=sys.stderr)
                sys.exit(1)
            run_management_command('list_project', rest[1])
        elif rest[0] == 'tree':
            if len(rest) < 2:
                print("Usage: i2 projects tree <project_id>", file=sys.stderr)
                sys.exit(1)
            run_management_command('tree_project', rest[1])
        else:
            # Assume it's a project ID for show
            run_management_command('list_project', rest[0])

    # ─────────────────────────────────────────────────────────────────────────
    # Jobs
    # ─────────────────────────────────────────────────────────────────────────
    elif resource == 'jobs':
        if not rest:
            print("Usage: i2 jobs <project_id> [list|create|run|tree|clone]", file=sys.stderr)
            sys.exit(1)

        # Check if first arg is an action or a project ID
        if rest[0] == 'list':
            if len(rest) < 2:
                print("Usage: i2 jobs list <project_id>", file=sys.stderr)
                sys.exit(1)
            run_management_command('list_jobs', rest[1])
        elif rest[0] == 'create':
            if len(rest) < 3:
                print("Usage: i2 jobs create <project_id> <task>", file=sys.stderr)
                sys.exit(1)
            run_management_command('create_job', *rest[1:])
        elif rest[0] == 'run':
            if len(rest) < 2:
                print("Usage: i2 jobs run <job_id>", file=sys.stderr)
                sys.exit(1)
            run_management_command('run_job', rest[1], *rest[2:])
        elif rest[0] == 'tree':
            if len(rest) < 2:
                print("Usage: i2 jobs tree <job_id>", file=sys.stderr)
                sys.exit(1)
            run_management_command('tree_job', rest[1])
        elif rest[0] == 'clone':
            if len(rest) < 2:
                print("Usage: i2 jobs clone <job_id>", file=sys.stderr)
                sys.exit(1)
            run_management_command('clone_job', rest[1], *rest[2:])
        else:
            # Assume first arg is project ID, list jobs in that project
            run_management_command('list_jobs', rest[0])

    # ─────────────────────────────────────────────────────────────────────────
    # Run (shortcut for i2run)
    # ─────────────────────────────────────────────────────────────────────────
    elif resource == 'run':
        if not rest:
            print("Usage: i2 run <task> [options]", file=sys.stderr)
            sys.exit(1)
        run_management_command('i2run', *rest)

    # ─────────────────────────────────────────────────────────────────────────
    # Files
    # ─────────────────────────────────────────────────────────────────────────
    elif resource == 'files':
        if not rest:
            print("Usage: i2 files <job_id> [list|cat]", file=sys.stderr)
            sys.exit(1)

        if rest[0] == 'list':
            if len(rest) < 2:
                print("Usage: i2 files list <job_id>", file=sys.stderr)
                sys.exit(1)
            run_management_command('list_files', rest[1])
        elif rest[0] == 'cat':
            if len(rest) < 3:
                print("Usage: i2 files cat <job_id> <filename>", file=sys.stderr)
                sys.exit(1)
            run_management_command('cat_job_file', rest[1], rest[2])
        else:
            # Assume it's a job ID, list files
            run_management_command('list_files', rest[0])

    # ─────────────────────────────────────────────────────────────────────────
    # Report
    # ─────────────────────────────────────────────────────────────────────────
    elif resource == 'report':
        if not rest:
            print("Usage: i2 report <job_id>", file=sys.stderr)
            sys.exit(1)
        run_management_command('get_job_report', rest[0], *rest[1:])

    # ─────────────────────────────────────────────────────────────────────────
    # Import/Export
    # ─────────────────────────────────────────────────────────────────────────
    elif resource == 'export':
        if not rest:
            print("Usage: i2 export <job|project> <id>", file=sys.stderr)
            sys.exit(1)
        if rest[0] == 'job':
            if len(rest) < 2:
                print("Usage: i2 export job <job_id>", file=sys.stderr)
                sys.exit(1)
            run_management_command('export_job', rest[1], *rest[2:])
        elif rest[0] == 'project':
            if len(rest) < 2:
                print("Usage: i2 export project <project_id>", file=sys.stderr)
                sys.exit(1)
            run_management_command('export_project', rest[1], *rest[2:])
        else:
            print(f"Unknown export target: {rest[0]}", file=sys.stderr)
            print("Usage: i2 export <job|project> <id>", file=sys.stderr)
            sys.exit(1)

    elif resource == 'import':
        if not rest:
            print("Usage: i2 import <zipfile>", file=sys.stderr)
            sys.exit(1)
        run_management_command('import_ccp4_project_zip', rest[0], *rest[1:])

    # ─────────────────────────────────────────────────────────────────────────
    # Unknown
    # ─────────────────────────────────────────────────────────────────────────
    else:
        print(f"Unknown resource: {resource}", file=sys.stderr)
        print("Available: projects, jobs, run, files, report, export, import", file=sys.stderr)
        print("Run 'i2 --help' for usage", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
