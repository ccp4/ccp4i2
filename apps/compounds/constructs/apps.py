# Copyright (C) 2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
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
