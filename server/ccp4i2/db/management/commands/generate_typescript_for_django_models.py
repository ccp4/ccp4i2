from django.core.management.base import BaseCommand
from django.apps import apps
from django.db.models.fields import (
    TextField,
    BigAutoField,
    IntegerField,
    FloatField,
    CharField,
    DateTimeField,
    UUIDField,
)
from django.db.models.fields.related import ForeignKey, ManyToManyField
from django.db.models.fields.reverse_related import (
    ManyToManyRel,
    ManyToOneRel,
    OneToOneRel,
)


class Command(BaseCommand):

    help = "Generate start point for template bindings"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument("app_label")

    def handle(self, **kwargs):
        app_label = kwargs.get("app_label", None)
        if app_label is None:
            print("Need to provide app_label")

        the_app = apps.get_app_config(app_label)
        the_models = the_app.get_models()
        for the_model in the_models:
            print(f"export class {the_model.__name__} {{\n    constructor(")
            for the_field in the_model._meta.get_fields():
                type = the_field.__class__
                if isinstance(the_field, TextField):
                    type = "string"
                if isinstance(the_field, CharField):
                    type = "string"
                if isinstance(the_field, UUIDField):
                    type = "string"
                if isinstance(the_field, DateTimeField):
                    type = "string"
                if isinstance(the_field, IntegerField):
                    type = "number"
                if isinstance(the_field, FloatField):
                    type = "number"
                if isinstance(the_field, ForeignKey):
                    type = "number"
                if isinstance(the_field, BigAutoField):
                    type = "number"
                if isinstance(the_field, ManyToOneRel):
                    type = "number[]"
                if isinstance(the_field, ManyToManyRel):
                    type = "number[]"
                if isinstance(the_field, ManyToManyField):
                    type = "number[]"
                if isinstance(the_field, OneToOneRel):
                    type = "number"
                print(f"        public {the_field.name}: {type},")
            print("    ){}\n}")
