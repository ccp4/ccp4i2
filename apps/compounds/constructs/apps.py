"""
Constructs App Configuration

Django app for managing plasmid constructs, proteins, cassettes,
and sequencing results.
"""

from django.apps import AppConfig


class ConstructsConfig(AppConfig):
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'compounds.constructs'
    label = 'constructs'
    verbose_name = 'Construct Database'

    def ready(self):
        # Import signal handlers
        from . import models  # noqa: F401
