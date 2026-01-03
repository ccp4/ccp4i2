from django.apps import AppConfig


class DbConfig(AppConfig):
    name = "ccp4i2.db"
    label = "ccp4i2"
    default_auto_field = "django.db.models.BigAutoField"

    def ready(self):
        print("status: ready")
        import ccp4i2.db.signals  # noqa
