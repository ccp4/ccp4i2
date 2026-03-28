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
from django.core.management.base import BaseCommand

from ccp4i2.db.import_i2xml import import_i2xml_from_file


class Command(BaseCommand):

    help = "Import a project"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument("project_xml")

    def handle(self, *args, **options):
        self.stdout.write(f"{options}")
        import_i2xml_from_file(options["program_xml"])
