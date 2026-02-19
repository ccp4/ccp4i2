from pathlib import Path
from django.utils.text import slugify
from django.conf import settings
from rest_framework.serializers import ModelSerializer, ValidationError, SerializerMethodField
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
        fields = ['id', 'text']

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
    """Lightweight serializer for project lists with tags.

    Excludes 'directory' to keep response size manageable.
    Tags use simplified serializer (id, text only) to avoid N+1 queries.
    """

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


class ProjectExportSerializer(ModelSerializer):
    class Meta:
        model = models.ProjectExport
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
    """Basic serializer for project groups with member count."""

    member_count = SerializerMethodField()

    class Meta:
        model = models.ProjectGroup
        fields = "__all__"


class ProjectGroupDetailSerializer(ModelSerializer):
    """Extended serializer with nested memberships for detail view."""

    memberships = ProjectGroupMembershipSerializer(many=True, read_only=True)

    class Meta:
        model = models.ProjectGroup
        fields = "__all__"
