from getpass import getuser
from socket import gethostname
from uuid import uuid4
from pathlib import Path

from django.db.models import (
    CASCADE,
    CharField,
    DateTimeField,
    FloatField,
    ForeignKey,
    IntegerChoices,
    IntegerField,
    JSONField,
    ManyToManyField,
    Model,
    OneToOneField,
    RESTRICT,
    SET_NULL,
    TextField,
    UUIDField,
    TextChoices,
)
from django.core.exceptions import ValidationError
from django.utils import timezone


class Project(Model):
    uuid = UUIDField(default=uuid4, unique=True)
    name = CharField(max_length=100, unique=True)
    description = TextField(blank=True)
    directory = TextField(unique=True)
    creation_time = DateTimeField(default=timezone.now)
    creation_user = TextField(default=getuser)
    creation_host = TextField(default=gethostname)
    last_access = DateTimeField(default=timezone.now)
    last_job_number = IntegerField(default=0)
    follow_from_job = ForeignKey(
        "Job", SET_NULL, blank=True, null=True, related_name="+"
    )
    i1_project_name = TextField(blank=True)
    i1_project_directory = TextField(blank=True)

    def __str__(self):
        return self.name


class ProjectGroup(Model):
    """A named set of projects that owns type-dependent state and behaviour.

    Use a ProjectGroup when the collection is itself a working object (a
    fragment campaign has sites, a summary scene, a PanDDA export). For pure
    organisation of the project list, use ProjectTag instead — see
    docs/organising-projects.md.
    """

    class GroupType(TextChoices):
        GENERAL_SET = "general_set", "General set"
        FRAGMENT_SET = "fragment_set", "Fragment set"

    name = CharField(max_length=100, unique=True)
    type = CharField(
        max_length=32, choices=GroupType.choices, default=GroupType.GENERAL_SET
    )

    # Convenience relation to access projects in a group
    projects = ManyToManyField(
        Project,
        related_name="groups",
        through="ProjectGroupMembership",
        through_fields=("group", "project"),
    )

    # Sites for fragment campaigns - stored view states for quick navigation
    # Schema: [{"name": "Site 1", "origin": [x, y, z], "quat": [x, y, z, w], "zoom": z}, ...]
    sites = JSONField(default=list, blank=True)

    def __str__(self):
        return self.name


class ProjectGroupMembership(Model):
    class MembershipType(TextChoices):
        PARENT = "parent", "Parent"
        MEMBER = "member", "Member"

    group = ForeignKey(ProjectGroup, CASCADE, related_name="memberships")
    project = ForeignKey(Project, CASCADE, related_name="group_memberships")
    type = CharField(max_length=16, choices=MembershipType.choices)

    class Meta:
        unique_together = ["group", "project"]

    def __str__(self):
        return f"{self.project} in {self.group} as {self.type}"


class ProjectTag(Model):
    """A label on a project. Carries no state beyond its text and its place in
    the tag forest (``parent``), which is what makes it safe to apply freely.
    Contrast ProjectGroup — see docs/organising-projects.md.

    ``path`` is the materialised chain of ancestor texts, maintained on save.
    It exists for two reasons: it gives the tree a real database-level
    uniqueness constraint (``unique_together = [parent, text]`` does not
    constrain rows with ``parent = NULL``, because SQL treats NULLs as
    distinct), and it turns descendant lookup — needed whenever a browser rolls
    projects up to an ancestor node — into a prefix query rather than recursive
    descent.

    The separator is an ASCII unit separator rather than "/" so that no
    restriction has to be placed on what a tag may be called: "A/B" as a root
    tag and "B" under a root "A" are distinct nodes and must not collide.
    """

    PATH_SEPARATOR = "\x1f"
    DISPLAY_SEPARATOR = "/"

    parent = ForeignKey("self", CASCADE, blank=True, null=True, related_name="children")
    text = CharField(max_length=50)
    path = CharField(max_length=1024, unique=True, editable=False, default="")
    projects = ManyToManyField(Project, related_name="tags", blank=True)

    class Meta:
        unique_together = ["parent", "text"]

    def __str__(self):
        return self.text

    # -- tree maintenance ------------------------------------------------

    def compute_path(self):
        """The materialised path this tag should have, given its parent."""
        if self.parent_id is None:
            return self.text
        return self.parent.path + self.PATH_SEPARATOR + self.text

    def check_no_cycle(self):
        """Raise if this tag's parent chain would loop back to itself."""
        seen = set()
        ancestor = self.parent
        while ancestor is not None:
            if ancestor.pk == self.pk:
                raise ValidationError(
                    f"Tag '{self.text}' cannot be a descendant of itself."
                )
            if ancestor.pk in seen:
                raise ValidationError("Tag hierarchy contains a cycle.")
            seen.add(ancestor.pk)
            ancestor = ancestor.parent
        return True

    def clean(self):
        super().clean()
        self.check_no_cycle()

    def save(self, *args, **kwargs):
        self.check_no_cycle()
        previous_path = None
        if self.pk is not None:
            previous_path = (
                type(self).objects.filter(pk=self.pk).values_list("path", flat=True).first()
            )
        self.path = self.compute_path()
        super().save(*args, **kwargs)
        # Re-parenting or renaming moves an entire subtree; rewrite it so the
        # descendants' paths stay consistent with their ancestry.
        if previous_path is not None and previous_path != self.path:
            self._rewrite_descendant_paths(previous_path)

    def _rewrite_descendant_paths(self, previous_path):
        prefix = previous_path + self.PATH_SEPARATOR
        descendants = type(self).objects.filter(path__startswith=prefix)
        for descendant in descendants:
            descendant.path = self.path + self.PATH_SEPARATOR + descendant.path[len(prefix):]
            super(ProjectTag, descendant).save(update_fields=["path"])

    # -- tree queries ----------------------------------------------------

    @property
    def display_path(self):
        """The path as a user would read it, e.g. ``SARS-CoV-2/Mpro``."""
        return self.DISPLAY_SEPARATOR.join(self.path.split(self.PATH_SEPARATOR))

    @property
    def depth(self):
        return self.path.count(self.PATH_SEPARATOR)

    def descendants(self):
        """Every tag below this one, at any depth."""
        return type(self).objects.filter(
            path__startswith=self.path + self.PATH_SEPARATOR
        )

    def self_and_descendants(self):
        from django.db.models import Q

        return type(self).objects.filter(
            Q(pk=self.pk) | Q(path__startswith=self.path + self.PATH_SEPARATOR)
        )

    def tagged_projects(self, include_descendants=True):
        """Projects under this node.

        With ``include_descendants`` (the default, and what a nested browser
        wants) a node answers for everything filed beneath it. The stored
        many-to-many is untouched either way: roll-up is a read-time concern,
        never something written back into the database.
        """
        if not include_descendants:
            return self.projects.all()
        return Project.objects.filter(
            tags__in=self.self_and_descendants()
        ).distinct()


class ProjectExport(Model):
    STATUS_PENDING = "pending"
    STATUS_RUNNING = "running"
    STATUS_COMPLETED = "completed"
    STATUS_FAILED = "failed"
    STATUS_CHOICES = [
        (STATUS_PENDING, "Pending"),
        (STATUS_RUNNING, "Running"),
        (STATUS_COMPLETED, "Completed"),
        (STATUS_FAILED, "Failed"),
    ]

    project = ForeignKey(Project, CASCADE, related_name="exports")
    time = DateTimeField(default=timezone.now)
    status = CharField(
        max_length=16, choices=STATUS_CHOICES, default=STATUS_PENDING
    )

    @property
    def file_exists(self):
        """Check if the export file exists on disk"""
        from django.utils.text import slugify
        import os

        project_name = slugify(self.project.name or f"project_{self.project.id}")
        timestamp = self.time.strftime("%Y%m%d_%H%M%S")
        export_file_name = f"{project_name}_export_{timestamp}.ccp4_project.zip"
        export_file_path = os.path.join(
            self.project.directory, "CCP4_EXPORT_FILES", export_file_name
        )
        return os.path.exists(export_file_path)

    def __str__(self):
        return f"{self.project} at {self.time}"


class ProjectImport(Model):
    project = ForeignKey(Project, CASCADE, related_name="imports")
    time = DateTimeField(default=timezone.now)

    def __str__(self):
        return f"{self.project} at {self.time}"


class Job(Model):
    class Status(IntegerChoices):
        UNKNOWN = 0, "Unknown"
        PENDING = 1, "Pending"
        QUEUED = 2, "Queued"
        RUNNING = 3, "Running"
        INTERRUPTED = 4, "Interrupted"
        FAILED = 5, "Failed"
        FINISHED = 6, "Finished"
        RUNNING_REMOTELY = 7, "Running remotely"
        FILE_HOLDER = 8, "File holder"
        TO_DELETE = 9, "To delete"
        UNSATISFACTORY = 10, "Unsatisfactory"

    class Evaluation(IntegerChoices):
        UNKNOWN = 0, "Unknown"
        BEST = 1, "Best"
        GOOD = 2, "Good"
        REJECTED = 3, "Rejected"

    uuid = UUIDField(default=uuid4, unique=True)
    project = ForeignKey(Project, CASCADE, related_name="jobs")
    parent = ForeignKey("self", CASCADE, blank=True, null=True, related_name="children")
    number = CharField(max_length=50)  # 1 or 1.1 etc
    title = CharField(max_length=255)
    status = IntegerField(choices=Status.choices, default=Status.UNKNOWN)
    evaluation = IntegerField(choices=Evaluation.choices, default=Evaluation.UNKNOWN)
    comments = TextField(blank=True)
    creation_time = DateTimeField(default=timezone.now)
    finish_time = DateTimeField(blank=True, null=True)
    task_name = CharField(max_length=100)
    process_id = IntegerField(blank=True, null=True)

    class Meta:
        unique_together = ["project", "number"]

    def __str__(self):
        return f"{self.number} {self.title}"

    @property
    def rel_path(self) -> str:
        path_elements = ["CCP4_JOBS"] + [
            f"job_{element}" for element in self.number.split(".")
        ]
        return "/".join(path_elements)

    @property
    def directory(self):
        path_elements = [f"job_{element}" for element in self.number.split(".")]
        jobs_dir = Path(self.project.directory) / "CCP4_JOBS"
        return jobs_dir.joinpath(*path_elements)


class ServerJob(Model):
    job = OneToOneField(Job, CASCADE, primary_key=True)
    server_process_id = IntegerField(blank=True, null=True)
    machine = CharField(max_length=255, blank=True)
    username = CharField(max_length=100, blank=True)
    mechanism = CharField(max_length=32, blank=True)
    remote_path = CharField(max_length=255, blank=True)
    custom_code_file = CharField(max_length=255, blank=True)
    validate = CharField(max_length=32, blank=True)
    key_file_name = CharField(max_length=255, blank=True)
    server_group = CharField(max_length=32, blank=True)

    def __str__(self):
        return str(self.job)


class JobValueKey(Model):
    name = CharField(max_length=50, primary_key=True)
    description = TextField()

    def __str__(self):
        return self.name


class JobFloatValue(Model):
    # If a job is deleted, all its float values should be deleted
    job = ForeignKey(Job, CASCADE, related_name="float_values")
    # Existence of a JobFloatValue should preclude deletion of the corresponding JobValueKey
    key = ForeignKey(JobValueKey, RESTRICT, related_name="+")
    value = FloatField()

    class Meta:
        unique_together = ["job", "key"]

    def __str__(self):
        return f"{self.job} {self.key} = {self.value}"


class JobCharValue(Model):
    # If a job is deleted, all its char values should be deleted
    job = ForeignKey(Job, CASCADE, related_name="char_values")
    # Existence of a JobFloatValue should preclude deletion of the corresponding JobValueKey
    key = ForeignKey(JobValueKey, RESTRICT, related_name="+")
    value = CharField(max_length=255)

    class Meta:
        unique_together = ["job", "key"]

    def __str__(self):
        return f"{self.job} {self.key} = {self.value}"


class FileType(Model):
    name = CharField(max_length=50, primary_key=True)
    description = TextField()

    def __str__(self):
        return self.name


class File(Model):
    class Directory(IntegerChoices):
        JOB_DIR = 1, "Job"
        IMPORT_DIR = 2, "Import"

    uuid = UUIDField(default=uuid4, unique=True)
    name = TextField()
    directory = IntegerField(choices=Directory.choices)
    type = ForeignKey(FileType, RESTRICT, related_name="files")
    sub_type = IntegerField(blank=True, null=True)
    content = IntegerField(blank=True, null=True)
    annotation = TextField(blank=True)
    job = ForeignKey(Job, CASCADE, blank=True, null=True, related_name="files")
    job_param_name = CharField(max_length=255, blank=True)

    def __str__(self):
        return self.name

    @property
    def path(self) -> Path:
        if self.directory == File.Directory.JOB_DIR:
            return self.job.directory / self.name
        elif self.directory == File.Directory.IMPORT_DIR:
            return Path(self.job.project.directory) / "CCP4_IMPORTED_FILES" / self.name

    @property
    def rel_path(self) -> str:
        if self.directory == File.Directory.JOB_DIR:
            return self.job.rel_path
        if self.directory == File.Directory.IMPORT_DIR:
            return "CCP4_IMPORTED_FILES"
        return ""


class FileExport(Model):
    file = ForeignKey(File, CASCADE, related_name="exports")
    time = DateTimeField(default=timezone.now)
    name = TextField()

    def __str__(self):
        return f"{self.name} at {self.time}"


class FileImport(Model):
    # MN: Note this has dropped fields annotation lastmodifiedtime importnumber
    file = OneToOneField(File, on_delete=CASCADE, primary_key=True)
    time = DateTimeField(default=timezone.now)
    name = TextField()
    # `checksum` is the checksum of `file` as it ended up in the project. For
    # anything derived --- an MTZ split down to the columns one parameter
    # wanted --- that is the checksum of the derivative, not of what the user
    # uploaded, so it cannot answer "have I seen this file before?".
    checksum = CharField(max_length=32)
    # `source_checksum` is the checksum of the bytes as uploaded, taken before
    # any splitting or conversion. Blank on rows written before this field
    # existed, and on rows whose source could not be hashed; blank never
    # matches, so those rows simply take no part in duplicate detection.
    source_checksum = CharField(max_length=32, blank=True, default="", db_index=True)
    last_modified = DateTimeField(blank=True, null=True)

    def __str__(self):
        return f"{self.name} at {self.time}"


class FileUse(Model):
    class Role(IntegerChoices):
        OUT = 0, "OUT"
        IN = 1, "IN"

    file = ForeignKey(File, CASCADE, related_name="file_uses")
    job = ForeignKey(Job, CASCADE, related_name="file_uses")
    role = IntegerField(choices=Role.choices)
    job_param_name = CharField(max_length=255, blank=True)

    class Meta:
        unique_together = ["file", "job", "role", "job_param_name"]

    def __str__(self):
        return self.job_param_name


class XData(Model):
    data_class = TextField(db_column="class")
    xml = TextField()
    job = ForeignKey(Job, RESTRICT, blank=True, null=True, related_name="xdatas")

    def __str__(self):
        return self.id
