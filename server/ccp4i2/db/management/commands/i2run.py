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
        # Don't add any arguments - we'll handle them manually in handle()
        pass

    def handle(self, *args, **options):
        """
        Use sys.argv directly to bypass Django's argument parsing.
        CCP4i2RunnerDjango will handle all argument parsing.

        Special handling for --i2run_configure flag:
        - If present, configure the job but do not execute it
        - Remove the flag before passing args to CCP4i2RunnerDjango
        """
        # sys.argv structure: ['manage.py', 'i2run', 'task_name', ...args...]
        # We want everything after 'i2run'
        the_args = sys.argv[2:]
        logger.info(f"i2run args: {the_args}")

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
                parent=None,  # No Qt spoof needed - using modern async approach
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
