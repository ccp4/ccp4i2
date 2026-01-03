"""Django's command-line utility for administrative tasks."""

from os import environ
from sys import argv, path
from pathlib import Path
from django.core.management import execute_from_command_line

# Add project root to Python path for imports from core/, pipelines/, etc.
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in path:
    path.insert(0, str(PROJECT_ROOT))

# Import ccp4dll early to set up DLL paths on Windows
try:
    import ccp4dll
except ImportError:
    pass


def main():
    """Run administrative tasks."""
    environ.setdefault("DJANGO_SETTINGS_MODULE", "ccp4i2.config.settings")
    execute_from_command_line(argv)


if __name__ == "__main__":
    main()
