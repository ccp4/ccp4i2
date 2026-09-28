import sys
import logging
import argparse
import traceback

from django.core.management.base import BaseCommand
from ccp4i2.cli.i2run.CCP4i2RunnerDjango import CCP4i2RunnerDjango
from xml.etree import ElementTree as ET

# Get an instance of a logger
logger = logging.getLogger("root")
logger.setLevel(logging.WARNING)


class Command(BaseCommand):

    help = "Configure and run a job in the database"
    requires_system_checks = []

    def add_arguments(self, parser):
        """Collect the task name and its parameters verbatim.

        A task's parameters come from its ``.def.xml`` and are only known once
        the task name has been read, so they cannot be declared here. But
        declaring *nothing* does not mean "accept anything": Django parses
        argv before ``handle()`` runs, so every argument was rejected --
        ``manage.py i2run freerflag --project_name gamma`` died with
        ``error: unrecognized arguments: freerflag --project_name gamma``.
        The command was therefore unusable from a shell, and only worked for
        the i2run test tier, which fakes ``sys.argv`` and calls
        ``call_command('i2run')`` with no arguments at all.

        ``REMAINDER`` takes the task name and everything after it untouched,
        for ``handle()`` to hand to the runner.
        """
        parser.add_argument(
            "i2run_args",
            nargs=argparse.REMAINDER,
            metavar="task_name [--PARAM value ...]",
            help="Task name followed by that task's own parameters.",
        )

    def handle(self, *args, **options):
        """
        Use sys.argv directly to bypass Django's argument parsing.
        CCP4i2RunnerDjango will handle all argument parsing.

        Special handling for --i2run_configure flag:
        - If present, configure the job but do not execute it
        - Remove the flag before passing args to CCP4i2RunnerDjango
        """
        # Normally the arguments arrive through the parser, in options. The
        # i2run test tier instead sets sys.argv and calls call_command with no
        # arguments (see tests/i2run/utils.py), so fall back to that shape.
        the_args = list(options.get("i2run_args") or [])
        if not the_args:
            # sys.argv structure: ['manage.py', 'i2run', 'task_name', ...args...]
            the_args = sys.argv[2:]
        logger.info(f"i2run args: {the_args}")

        if not the_args:
            self.stderr.write(
                "i2run needs a task name, e.g. "
                "manage.py i2run freerflag --project_name my_project"
            )
            return

        # Check for --i2run_configure flag and remove it from args
        configure_only = False
        if '--i2run_configure' in the_args:
            configure_only = True
            the_args = [arg for arg in the_args if arg != '--i2run_configure']
            logger.info("--i2run_configure flag detected: will configure but not execute")

        try:
            parser = argparse.ArgumentParser()

            # Modern approach: No Qt parent needed
            self.i2_runner = CCP4i2RunnerDjango(
                the_args=the_args,
                parser=parser,
            )

            self.i2_runner.parseArgs()

            # Execute only if not in configure-only mode
            if not configure_only:
                result = self.i2_runner.execute()
                logger.warning(f"i2run execute() returned: {result}")
            else:
                logger.warning("Skipping execute() due to --i2run_configure flag")
                thePlugin = self.i2_runner.getPlugin(arguments_parsed=True)
                plugin_etree = thePlugin.getEtree(excludeUnset=True)
                plugin_xml_str = ET.tostring(plugin_etree, encoding='unicode')
                logger.warning(f"Configured plugin XML:\n{plugin_xml_str}")
                result = None

        except Exception as e:
            logger.error(f"i2run failed with exception: {e}")
            logger.error(f"Traceback:\n{traceback.format_exc()}")
            print(f"\nERROR: i2run failed")
            print(f"Exception: {e}")
            print(f"\nTraceback:")
            print(traceback.format_exc())
            raise
