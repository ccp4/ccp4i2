# Copyright (C) 2025-2026 University of York
# Copyright (C) 2025-2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
from django.db.models.signals import pre_save
from django.dispatch import receiver
from .models import Job


@receiver(pre_save, sender=Job)
def job_status_change_handler(sender, instance, **kwargs):
    if not instance.pk:
        # New Job, not an update
        return
    try:
        old_instance = Job.objects.get(pk=instance.pk)
    except Job.DoesNotExist:
        return
    if old_instance.status != instance.status:
        # Status has changed
        print(
            f"Job {instance.pk} status changed from {old_instance.status} to {instance.status}"
        )
        # Place your custom logic here
