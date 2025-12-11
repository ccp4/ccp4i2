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
