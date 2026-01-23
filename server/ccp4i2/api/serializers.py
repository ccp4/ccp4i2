from pathlib import Path
from django.utils.text import slugify
from django.conf import settings
from rest_framework.serializers import ModelSerializer, ValidationError
from ..db import models


class FileTypeSerializer(ModelSerializer):
    class Meta:
        model = models.FileType
        fields = "__all__"


class FileImportSerializer(ModelSerializer):
    class Meta:
        model = models.FileImport
        fields = "__all__"


class ProjectTagSerializer(ModelSerializer):
    class Meta:
        model = models.ProjectTag
        exclude = []

    def validate(self, attrs):
        """Validate unique constraint on text and parent combination."""
        text = attrs.get("text")
        parent = attrs.get("parent")

        # Check for existing tag with same text and parent
        existing_tag = models.ProjectTag.objects.filter(
            text=text, parent=parent
        ).first()

        if existing_tag and (not self.instance or existing_tag.id != self.instance.id):
            raise ValidationError(
                {"text": "A tag with this text and parent already exists."}
            )

        return attrs


class ProjectListSerializer(ModelSerializer):
    """Lightweight serializer for project lists with prefetched tags.

    Excludes 'directory' field to reduce response size - long paths repeated
    across many projects add significant bytes. Directory is available in
    detail view via ProjectSerializer.
    """

    # Use nested serializer for tags - requires prefetch_related('tags') in ViewSet
    tags = ProjectTagSerializer(many=True, read_only=True)

    class Meta:
        model = models.Project
        fields = [
            "id",
            "uuid",
            "name",
            "creation_time",
            "last_access",
            "tags",
        ]


class ProjectSerializer(ModelSerializer):
    # Include tag details in project serialization
    tags = ProjectTagSerializer(many=True, read_only=True)

    class Meta:
        model = models.Project
        fields = "__all__"

    def validate(self, attrs):
        # Validation will depend on whether this is a new project or an update
        # If this is a new project (i.e. no existing instance), we need to provide a default for the directory
        instance = (
            self.instance
        )  # This is the instance being updated (or None if creating)

        if instance is None:
            if (
                "directory" not in attrs
                or not attrs["directory"]
                or len(attrs["directory"]) == 0
                or attrs["directory"] == "__default__"
            ):
                attrs["directory"] = str(
                    Path(settings.CCP4I2_PROJECTS_DIR) / slugify(attrs["name"])
                )
        return super().validate(attrs)

    def create(self, validated_data):

        Path(validated_data["directory"]).mkdir(parents=True, exist_ok=True)

        for sub_dir in [
            "CCP4_JOBS",
            "CCP4_IMPORTED_FILES",
            "CCP4_COOT",
            "CCP4_TMP",
            "CCP4_PROJECT_FILES",
        ]:
            (Path(validated_data["directory"]) / sub_dir).mkdir(exist_ok=True)

        return models.Project.objects.create(**validated_data)

    def validate_name(self, data: str):
        if any((not c.isalnum() and c not in ["_", "-"]) for c in data):
            raise ValidationError(
                f"Your project name contains whitespace or special characters [{data}]"
            )
        project_names = [
            project.name.upper() for project in models.Project.objects.all()
        ]
        if "uuid" not in self.initial_data and data.upper() in project_names:
            raise ValidationError("A project with this name already exists!")
        if "directory" not in self.initial_data:
            assert Path(settings.CCP4I2_PROJECTS_DIR).is_dir()
            try:
                testWritePath = Path(settings.CCP4I2_PROJECTS_DIR) / "testWrite.txt"
                with open(testWritePath, "w") as testWrite:
                    testWrite.write("test")
                testWritePath.unlink()
            except Exception as err:
                raise ValidationError(
                    f"Failure trying to write to  [{testWritePath}], {err}"
                ) from err
        return data


class FileSerializer(ModelSerializer):
    class Meta:
        model = models.File
        fields = "__all__"


class JobSerializer(ModelSerializer):
    class Meta:
        model = models.Job
        fields = "__all__"


class FileUseSerializer(ModelSerializer):
    class Meta:
        model = models.FileUse
        fields = "__all__"


class FileExportSerializer(ModelSerializer):
    class Meta:
        model = models.FileExport
        fields = "__all__"


class ProjectExportSerializer(ModelSerializer):
    class Meta:
        model = models.ProjectExport
        fields = "__all__"


class XDataSerializer(ModelSerializer):
    class Meta:
        model = models.XData
        fields = "__all__"


class JobFloatValueSerializer(ModelSerializer):
    class Meta:
        model = models.JobFloatValue
        fields = "__all__"


class JobCharValueSerializer(ModelSerializer):
    class Meta:
        model = models.JobCharValue
        fields = "__all__"


class ProjectGroupMembershipSerializer(ModelSerializer):
    """Serializer for project group memberships."""

    class Meta:
        model = models.ProjectGroupMembership
        fields = "__all__"


class ProjectGroupSerializer(ModelSerializer):
    """Basic serializer for project groups."""

    class Meta:
        model = models.ProjectGroup
        fields = "__all__"


class ProjectGroupDetailSerializer(ModelSerializer):
    """Extended serializer with nested memberships for detail view."""

    memberships = ProjectGroupMembershipSerializer(many=True, read_only=True)

    class Meta:
        model = models.ProjectGroup
        fields = "__all__"
