from django.apps import AppConfig


class AzureExtensionsConfig(AppConfig):
    default_auto_field = "django.db.models.BigAutoField"
    name = "azure_extensions"
    verbose_name = "Azure Extensions"

    def ready(self):
        # Import signal handlers or perform startup tasks here if needed
        pass
