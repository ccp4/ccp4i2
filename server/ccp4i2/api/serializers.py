from pathlib import Path
from django.utils.text import slugify
from django.conf import settings
from rest_framework.serializers import (
    ModelSerializer,
    ReadOnlyField,
    SerializerMethodField,
    ValidationError,
)
from ..db import models


class FileTypeSerializer(ModelSerializer):
    class Meta:
        model = models.FileType
        fields = "__all__"


class FileImportSerializer(ModelSerializer):
    class Meta:
        model = models.FileImport
        fields = "__all__"


class ProjectTagListSerializer(ModelSerializer):
    """Lightweight tag serializer for embedding in project lists (no reverse projects)."""

    display_path = ReadOnlyField()

    class Meta:
        model = models.ProjectTag
        fields = ['id', 'text', 'parent', 'path', 'display_path']


class ProjectTagSerializer(ModelSerializer):
    display_path = ReadOnlyField()

    class Meta:
        model = models.ProjectTag
        fields = ['id', 'text', 'parent', 'path', 'display_path', 'projects']

    def validate(self, attrs):
        """Guard the two ways an edit can break the tag forest.

        Renaming and re-parenting both arrive as partial updates, so the value
        not being edited has to come from the instance — reading it out of
        ``attrs`` alone would compare against ``None`` and wave through a
        collision that the unique path then rejects with a 500.
        """
        text = attrs.get("text", getattr(self.instance, "text", None))
        if "parent" in attrs:
            parent = attrs["parent"]
        else:
            parent = getattr(self.instance, "parent", None)

        # A tag cannot be moved beneath itself or its own descendant; the model
        # refuses it too, but only by raising on save, which is a 500.
        if parent is not None and self.instance is not None:
            ancestor = parent
            seen = set()
            while ancestor is not None and ancestor.pk not in seen:
                if ancestor.pk == self.instance.pk:
                    raise ValidationError(
                        {"parent": "A tag cannot be moved beneath itself."}
                    )
                seen.add(ancestor.pk)
                ancestor = ancestor.parent

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

    tags = ProjectTagListSerializer(many=True, read_only=True)

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


def default_project_parent() -> Path:
    """Where a project with no explicit directory should go.

    The parent of the most recently created project, falling back to the
    configured projects directory. So a user who put their last project
    somewhere particular gets offered the same place again, WITHOUT that choice
    being written into a preference: the "default" is derived on demand rather
    than stored and mutated. A stored-and-mutated default is how a one-off
    choice silently became everybody's default, and how a second database could
    appear somewhere unexpected.

    A project whose recorded directory no longer exists is skipped -- an
    unplugged external disk should not send the next project somewhere
    unwritable.
    """
    from ..db.models import Project

    for directory in (
        Project.objects.exclude(directory="")
        .order_by("-creation_time")
        .values_list("directory", flat=True)[:10]
    ):
        parent = Path(directory).parent
        if parent.is_dir():
            return parent
    return Path(settings.CCP4I2_PROJECTS_DIR)


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
                    default_project_parent() / slugify(attrs["name"])
                )
        elif "directory" in attrs and attrs["directory"] != instance.directory:
            # Changing this field alone would leave the record pointing at a
            # directory that is still sitting somewhere else, with every job's
            # baked-in absolute path stale. Relocating a project is a job for
            # POST /projects/{id}/move/, which moves the bytes and rebases the
            # paths as one operation.
            raise ValidationError(
                {
                    "directory": (
                        "A project's directory cannot be changed directly. Use "
                        "the project move endpoint, which relocates the files "
                        "and updates the paths inside them."
                    )
                }
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
        # An update has to be allowed to send back the name it already has: an
        # edit form that PATCHes every field, changed or not, would otherwise
        # collide with itself.
        others = models.Project.objects.all()
        if self.instance is not None:
            others = others.exclude(pk=self.instance.pk)
        project_names = [project.name.upper() for project in others]
        if "uuid" not in self.initial_data and data.upper() in project_names:
            raise ValidationError("A project with this name already exists!")
        # Only a project about to be created needs somewhere under
        # CCP4I2_PROJECTS_DIR to be created in. One that already exists has a
        # directory of its own, which need not be under that root at all, and
        # renaming it does not move it.
        if self.instance is None and "directory" not in self.initial_data:
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
    path = SerializerMethodField()

    class Meta:
        model = models.File
        fields = "__all__"

    def get_path(self, obj):
        try:
            return str(obj.path)
        except Exception:
            return None


class JobSerializer(ModelSerializer):
    float_values = SerializerMethodField()
    char_values = SerializerMethodField()

    class Meta:
        model = models.Job
        fields = "__all__"

    def get_float_values(self, obj):
        # JobValueKey.name is the PK, so kv.key_id is the KPI name string.
        return {kv.key_id: kv.value for kv in obj.float_values.all()}

    def get_char_values(self, obj):
        return {kv.key_id: kv.value for kv in obj.char_values.all()}


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
    """Basic serializer for project groups with member count."""

    member_count = SerializerMethodField()

    class Meta:
        model = models.ProjectGroup
        fields = "__all__"

    def get_member_count(self, obj):
        return obj.memberships.filter(
            type=models.ProjectGroupMembership.MembershipType.MEMBER
        ).count()


class ProjectGroupDetailSerializer(ModelSerializer):
    """Extended serializer with nested memberships for detail view."""

    memberships = ProjectGroupMembershipSerializer(many=True, read_only=True)

    class Meta:
        model = models.ProjectGroup
        fields = "__all__"
