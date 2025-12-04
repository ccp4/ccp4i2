from django.core.management.base import BaseCommand

from ccp4x.db.import_i2xml import import_i2xml_from_file


class Command(BaseCommand):

    help = "Import a project"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument("project_xml")

    def handle(self, *args, **options):
        self.stdout.write(f"{options}")
        import_i2xml_from_file(options["program_xml"])
