#!/usr/bin/env python3
"""
i2run - CCP4i2 command-line job runner

Main entry point for the i2run command-line tool.

Usage:
    ccp4-python -m ccp4i2.cli.i2run <task> [options]

Examples:
    ccp4-python -m ccp4i2.cli.i2run freerflag --project_name gamma \\
        --F_SIGF fullPath=merged_intensities.mtz --FRAC 0.05

This is the entry point the job panel's "i2run command" button renders, so it
has to work in a dev checkout and in a pip-installed (packaged) tree alike --
neither of which can rely on ``manage.py`` being on disk.
"""

import argparse
import os
import sys


def _setup_django() -> None:
    """Configure Django before anything imports a model.

    i2run is a Django application wearing a CLI: the runner, the job models
    and the task registry all read ``django.conf.settings``. Nothing else on
    this path does the bootstrap that ``manage.py`` does, so without this the
    console script and ``-m`` invocation both died on
    ``ImproperlyConfigured: Requested setting INSTALLED_APPS``.

    ``setdefault`` so an explicitly exported settings module still wins (the
    cloud deployments set their own).
    """
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "ccp4i2.config.settings")
    import django

    django.setup()


def main():
    """
    Main entry point for i2run.

    Sets up Django, then delegates to the Django-backed runner.
    """
    _setup_django()

    # Imported after django.setup(): the runner pulls in ccp4i2.db.models.
    from .CCP4i2RunnerDjango import CCP4i2RunnerDjango as Runner

    args = sys.argv[1:]
    if not args or args[0] in ("-h", "--help"):
        print(__doc__.strip())
        print("\nRun with a task name to see that task's own parameters, e.g.")
        print("    ccp4-python -m ccp4i2.cli.i2run freerflag --help")
        return 0

    # sys.argv[1:] is passed through as a list rather than joined into a string
    # and re-split: ' '.join(...) followed by shlex.split() loses the shell's
    # quoting, so `--jobTitle "my first run"` arrived as three arguments.
    runner = Runner(the_args=args, parser=argparse.ArgumentParser(prog="i2run"))

    # getPlugin(arguments_parsed=True), reached from execute(), assumes the
    # task's own arguments have already been added to the parser and parsed.
    runner.parseArgs()

    job_id, exit_code = runner.execute()
    sys.exit(exit_code)


if __name__ == '__main__':
    main()
