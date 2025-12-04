from django.apps import AppConfig


class DbConfig(AppConfig):
    name = "ccp4x.db"
    label = "ccp4x"
    default_auto_field = "django.db.models.BigAutoField"

    def ready(self):
        print("status: ready")
        import ccp4x.db.signals  # noqa
