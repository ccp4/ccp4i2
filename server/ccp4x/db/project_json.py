import os
from django.core.management import call_command
from django.db import models
from io import StringIO


def project_json(project_instance):
    """
    Exports a JSON file representing a dump of the database contents for a given project instance,
    all related objects through ForeignKey relationships, and recursively all objects related to those.

    Args:
        project_instance (models.Model): The project instance to export.
    """
    if not isinstance(project_instance, models.Model):
        raise ValueError("project_instance must be a Django model instance.")

    def collect_related_objects(model_instance, collected_models, collected_pks):
        """
        Recursively collect related objects through ForeignKey relationships.
        """
        model = model_instance.__class__
        pk = model_instance.pk
        model_identifier = f"{model._meta.app_label}.{model._meta.model_name}"

        if (model_identifier, pk) in collected_pks:
            return  # Avoid circular references

        collected_models.add(model_identifier)
        collected_pks.add((model_identifier, pk))

        for field in model._meta.get_fields():
            if isinstance(field, models.ForeignKey):
                related_instance = getattr(model_instance, field.name, None)
                if related_instance:
                    collect_related_objects(
                        related_instance, collected_models, collected_pks
                    )
            elif field.one_to_many or field.one_to_one:
                related_manager = getattr(model_instance, field.name, None)
                if related_manager:
                    if hasattr(
                        related_manager, "all"
                    ):  # For reverse ForeignKey relationships
                        for related_instance in related_manager.all():
                            collect_related_objects(
                                related_instance, collected_models, collected_pks
                            )
                    else:  # For one-to-one relationships
                        collect_related_objects(
                            related_manager, collected_models, collected_pks
                        )

    # Collect all related objects recursively
    collected_models = set()
    collected_pks = set()
    collect_related_objects(project_instance, collected_models, collected_pks)

    # Use Django's dumpdata management command to export the data to a string
    output = StringIO()
    try:
        for model in collected_models:
            pks = [str(pk) for m, pk in collected_pks if m == model]
            call_command(
                "dumpdata",
                model,
                format="json",
                indent=2,
                stdout=output,
                primary_keys=",".join(pks),
            )
        return output.getvalue()
    except Exception as e:
        print(f"An error occurred while exporting data: {e}")
        return None
