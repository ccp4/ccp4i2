"""
Django management command to import legacy CCP4i2 dumpdata fixtures.

This command imports data from the legacy CCP4i2 schema (CCP4i2Docker) into
the new CCP4i2 schema, handling the key differences:
- Project parent relationships -> ProjectGroupMembership
- UUID1 hex strings -> UUID4 UUIDField
- Float timestamps -> DateTimeField
- Lookup tables -> IntegerChoices/string PKs

Usage:
    python manage.py import_legacy_ccp4i2 legacy_dump.json
    python manage.py import_legacy_ccp4i2 legacy_dump.json --dry-run --verbose
    python manage.py import_legacy_ccp4i2 legacy_dump.json --remap-dirs /old/path /new/path
"""

import json
from collections import defaultdict
from datetime import datetime
from uuid import UUID

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction
from django.utils import timezone

from ccp4i2.db.models import (
    File,
    FileExport,
    FileImport,
    FileType,
    FileUse,
    Job,
    JobCharValue,
    JobFloatValue,
    JobValueKey,
    Project,
    ProjectExport,
    ProjectGroup,
    ProjectGroupMembership,
    ProjectImport,
    ProjectTag,
    ServerJob,
    XData,
)


class Command(BaseCommand):
    """Import legacy CCP4i2 dumpdata fixtures into new schema."""

    help = "Import legacy CCP4i2 dumpdata fixtures"
    requires_system_checks = []

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # ID mapping caches (legacy_id -> new_object)
        self.project_map = {}  # legacy projectid -> Project
        self.job_map = {}  # legacy jobid -> Job
        self.file_map = {}  # legacy fileid -> File
        self.filetype_map = {}  # legacy filetypeid -> FileType
        self.keytype_map = {}  # legacy keytypeid -> JobValueKey
        self.tag_map = {}  # legacy tagid -> ProjectTag
        self.filerole_map = {}  # legacy roleid -> FileUse.Role value

        # Track parent relationships for group creation
        self.parent_relationships = []  # [(child_project, parent_project), ...]

        # Statistics
        self.stats = defaultdict(int)

    def add_arguments(self, parser):
        """Add command-line arguments."""
        parser.add_argument(
            "fixture_file",
            type=str,
            help="Path to legacy CCP4i2 dumpdata JSON fixture file",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Validate without committing changes",
        )
        parser.add_argument(
            "--verbose",
            action="store_true",
            help="Show detailed progress",
        )
        parser.add_argument(
            "--remap-dirs",
            nargs=2,
            metavar=("FROM", "TO"),
            help="Remap project directories (e.g., --remap-dirs /old/path /new/path)",
        )
        parser.add_argument(
            "--projects-only",
            action="store_true",
            help="Import only projects (skip jobs/files for testing)",
        )
        parser.add_argument(
            "--continue-on-error",
            action="store_true",
            help="Continue importing remaining records if one fails",
        )

    def handle(self, *args, **options):
        """Handle the import_legacy_ccp4i2 command."""
        fixture_file = options["fixture_file"]
        self.dry_run = options["dry_run"]
        self.verbose = options["verbose"]
        self.remap_dirs = options["remap_dirs"]
        self.projects_only = options["projects_only"]
        self.continue_on_error = options["continue_on_error"]

        try:
            with open(fixture_file, "r") as f:
                content = f.read()

            # Strip any garbage before the JSON array starts
            # Legacy dumpdata may have stdout noise at the beginning
            json_start = content.find("[")
            if json_start == -1:
                raise CommandError(f"No JSON array found in fixture file: {fixture_file}")
            if json_start > 0:
                self.stdout.write(
                    self.style.WARNING(f"Skipping {json_start} bytes of non-JSON content at start of file")
                )
                content = content[json_start:]

            data = json.loads(content)
        except FileNotFoundError:
            raise CommandError(f"Fixture file not found: {fixture_file}")
        except json.JSONDecodeError as e:
            raise CommandError(f"Invalid JSON in fixture file: {e}")

        self.stdout.write(f"\nImporting legacy CCP4i2 fixture: {fixture_file}")
        self.stdout.write("-" * 60)

        # Group records by model
        by_model = defaultdict(list)
        for item in data:
            model_name = item.get("model", "").lower()
            by_model[model_name].append(item)

        if self.verbose:
            self.stdout.write(f"\nModels found in fixture:")
            for model, records in sorted(by_model.items()):
                self.stdout.write(f"  {model}: {len(records)} records")

        try:
            if self.dry_run:
                self.stdout.write(self.style.WARNING("\n[DRY RUN - no changes will be committed]\n"))

            with transaction.atomic():
                # Phase 1: Lookup tables
                self.stdout.write("\nPhase 1: Lookup tables")
                self.import_filetypes(by_model.get("ccp4i2.filetypes", []))
                self.import_keytypes(by_model.get("ccp4i2.keytypes", []))
                self.import_fileroles(by_model.get("ccp4i2.fileroles", []))
                self.import_tags(by_model.get("ccp4i2.tags", []))

                # Phase 2: Projects
                self.stdout.write("\nPhase 2: Projects")
                self.import_projects(by_model.get("ccp4i2.projects", []))
                self.import_projecttags(by_model.get("ccp4i2.projecttags", []))

                # Phase 3: Project Groups (parent relationships)
                self.stdout.write("\nPhase 3: Project Groups (parent relationships)")
                self.create_project_groups()

                if not self.projects_only:
                    # Phase 4: Jobs
                    self.stdout.write("\nPhase 4: Jobs")
                    self.import_jobs(by_model.get("ccp4i2.jobs", []))
                    self.import_serverjobs(by_model.get("ccp4i2.serverjobs", []))
                    self.import_jobkeyvalues(by_model.get("ccp4i2.jobkeyvalues", []))
                    self.import_jobkeycharvalues(by_model.get("ccp4i2.jobkeycharvalues", []))
                    self.import_xdata(by_model.get("ccp4i2.xdata", []))

                    # Phase 5: Files
                    self.stdout.write("\nPhase 5: Files")
                    self.import_files(by_model.get("ccp4i2.files", []))
                    self.import_importfiles(by_model.get("ccp4i2.importfiles", []))
                    self.import_exportfiles(by_model.get("ccp4i2.exportfiles", []))
                    self.import_fileuses(by_model.get("ccp4i2.fileuses", []))

                    # Phase 6: Project history
                    self.stdout.write("\nPhase 6: Project history")
                    self.import_projectexports(by_model.get("ccp4i2.projectexports", []))
                    self.import_projectimports(by_model.get("ccp4i2.projectimports", []))

                if self.dry_run:
                    raise DryRunRollback()

        except DryRunRollback:
            self.stdout.write(self.style.WARNING("\n[DRY RUN - all changes rolled back]\n"))

        # Print summary
        self.print_summary()

    # -------------------------------------------------------------------------
    # Utility functions
    # -------------------------------------------------------------------------

    def convert_uuid(self, legacy_uuid_hex):
        """Convert 32-char UUID1 hex string to UUID object."""
        if not legacy_uuid_hex:
            return None
        try:
            # UUID hex strings can be parsed directly
            return UUID(legacy_uuid_hex)
        except (ValueError, TypeError):
            return None

    def convert_timestamp(self, epoch_float):
        """Convert epoch float to timezone-aware datetime."""
        if epoch_float is None:
            return None
        try:
            dt = datetime.fromtimestamp(float(epoch_float))
            return timezone.make_aware(dt, timezone.get_default_timezone())
        except (ValueError, TypeError, OSError):
            return None

    def remap_directory(self, directory):
        """Apply directory remapping if configured."""
        if not directory or not self.remap_dirs:
            return directory
        from_path, to_path = self.remap_dirs
        if directory.startswith(from_path):
            return directory.replace(from_path, to_path, 1)
        return directory

    def get_unique_project_name(self, name):
        """Get unique project name by adding suffix if needed."""
        if not Project.objects.filter(name=name).exists():
            return name

        # Try _legacy suffix first
        legacy_name = f"{name}_legacy"
        if not Project.objects.filter(name=legacy_name).exists():
            self.stats["projects_renamed"] += 1
            return legacy_name

        # Try numbered suffixes
        for i in range(2, 1000):
            numbered_name = f"{name}_legacy_{i}"
            if not Project.objects.filter(name=numbered_name).exists():
                self.stats["projects_renamed"] += 1
                return numbered_name

        raise CommandError(f"Unable to find unique name for project: {name}")

    def log(self, message):
        """Log message if verbose mode is enabled."""
        if self.verbose:
            self.stdout.write(f"  {message}")

    def log_count(self, entity, count, extra=""):
        """Log import count for an entity."""
        extra_msg = f" ({extra})" if extra else ""
        self.stdout.write(f"  - {entity}: {count} imported{extra_msg}")

    # -------------------------------------------------------------------------
    # Import functions
    # -------------------------------------------------------------------------

    def import_filetypes(self, records):
        """Import Filetypes -> FileType."""
        for record in records:
            pk = record.get("pk")
            fields = record.get("fields", {})
            name = fields.get("filetypename", f"type_{pk}")
            description = fields.get("filetypedescription", "")

            filetype, created = FileType.objects.get_or_create(
                name=name,
                defaults={"description": description or ""}
            )
            self.filetype_map[pk] = filetype
            if created:
                self.stats["filetypes"] += 1

        self.log_count("FileTypes", self.stats["filetypes"])

    def import_keytypes(self, records):
        """Import Keytypes -> JobValueKey."""
        for record in records:
            pk = record.get("pk")
            fields = record.get("fields", {})
            name = fields.get("keytypename", f"key_{pk}")
            description = fields.get("keytypedescription", "")

            keytype, created = JobValueKey.objects.get_or_create(
                name=name,
                defaults={"description": description or ""}
            )
            self.keytype_map[pk] = keytype
            if created:
                self.stats["keytypes"] += 1

        self.log_count("JobValueKeys", self.stats["keytypes"])

    def import_fileroles(self, records):
        """Import Fileroles -> map to FileUse.Role IntegerChoices."""
        # Map legacy role IDs to new Role enum values
        # Legacy: roletext values like "IN", "OUT"
        for record in records:
            pk = record.get("pk")
            fields = record.get("fields", {})
            roletext = fields.get("roletext", "").upper()

            if roletext == "IN":
                self.filerole_map[pk] = FileUse.Role.IN
            elif roletext == "OUT":
                self.filerole_map[pk] = FileUse.Role.OUT
            else:
                # Default to OUT for unknown roles
                self.filerole_map[pk] = FileUse.Role.OUT

        self.log_count("FileRoles mapped", len(self.filerole_map))

    def import_tags(self, records):
        """Import Tags -> ProjectTag (hierarchical, two passes)."""
        # First pass: create all tags without parent references
        tag_parent_map = {}  # new_tag_id -> legacy_parent_tagid
        for record in records:
            pk = record.get("pk")
            fields = record.get("fields", {})
            text = fields.get("text", f"tag_{pk}")
            parent_tagid = fields.get("parenttagid")

            # Create without parent first
            tag, created = ProjectTag.objects.get_or_create(
                text=text,
                parent=None,
                defaults={}
            )
            self.tag_map[pk] = tag
            if parent_tagid:
                tag_parent_map[tag.id] = parent_tagid
            if created:
                self.stats["tags"] += 1

        # Second pass: set parent references
        for tag_id, legacy_parent_id in tag_parent_map.items():
            if legacy_parent_id in self.tag_map:
                tag = ProjectTag.objects.get(id=tag_id)
                tag.parent = self.tag_map[legacy_parent_id]
                tag.save()

        self.log_count("ProjectTags", self.stats["tags"])

    def import_projects(self, records):
        """Import Projects -> Project."""
        for record in records:
            try:
                pk = record.get("pk")
                fields = record.get("fields", {})

                name = fields.get("projectname", f"project_{pk}")
                unique_name = self.get_unique_project_name(name)

                directory = fields.get("projectdirectory", "")
                directory = self.remap_directory(directory)

                project = Project.objects.create(
                    uuid=self.convert_uuid(pk) or None,
                    name=unique_name,
                    description="",  # Legacy doesn't have description
                    directory=directory or "",
                    creation_time=self.convert_timestamp(fields.get("projectcreated")) or timezone.now(),
                    creation_user=fields.get("userid", "legacy_import") or "legacy_import",
                    creation_host="legacy_import",
                    last_access=self.convert_timestamp(fields.get("lastaccess")) or timezone.now(),
                    last_job_number=int(fields.get("lastjobnumber", 0) or 0),
                    i1_project_name=fields.get("i1projectname", "") or "",
                    i1_project_directory=fields.get("i1projectdirectory", "") or "",
                )
                self.project_map[pk] = project
                self.stats["projects"] += 1

                # Track parent relationship for later
                parent_id = fields.get("parentprojectid")
                if parent_id:
                    self.parent_relationships.append((pk, parent_id))

                self.log(f"Project: {unique_name}")

            except Exception as e:
                if self.continue_on_error:
                    self.stderr.write(self.style.ERROR(f"  Error importing project {pk}: {e}"))
                else:
                    raise

        extra = f"{self.stats['projects_renamed']} renamed" if self.stats["projects_renamed"] else ""
        self.log_count("Projects", self.stats["projects"], extra)

    def import_projecttags(self, records):
        """Import Projecttags -> Project.tags M2M."""
        for record in records:
            fields = record.get("fields", {})
            project_id = fields.get("projectid")
            tag_id = fields.get("tagid")

            if project_id in self.project_map and tag_id in self.tag_map:
                project = self.project_map[project_id]
                tag = self.tag_map[tag_id]
                project.tags.add(tag)
                self.stats["projecttags"] += 1

        self.log_count("Project-Tag associations", self.stats["projecttags"])

    def create_project_groups(self):
        """Create ProjectGroups for parent relationships."""
        # Group children by parent
        parent_children = defaultdict(list)
        for child_id, parent_id in self.parent_relationships:
            if child_id in self.project_map and parent_id in self.project_map:
                parent_children[parent_id].append(child_id)

        for parent_id, child_ids in parent_children.items():
            parent_project = self.project_map[parent_id]

            # Create a group for this parent
            group_name = f"Legacy_{parent_project.name}_group"
            # Ensure unique group name
            base_name = group_name
            counter = 1
            while ProjectGroup.objects.filter(name=group_name).exists():
                group_name = f"{base_name}_{counter}"
                counter += 1

            group = ProjectGroup.objects.create(
                name=group_name,
                type=ProjectGroup.GroupType.FRAGMENT_SET  # Default to fragment_set for legacy imports
            )
            self.stats["groups"] += 1

            # Add parent as PARENT type
            ProjectGroupMembership.objects.create(
                group=group,
                project=parent_project,
                type=ProjectGroupMembership.MembershipType.PARENT
            )
            self.stats["memberships"] += 1

            # Add children as MEMBER type
            for child_id in child_ids:
                child_project = self.project_map[child_id]
                ProjectGroupMembership.objects.create(
                    group=group,
                    project=child_project,
                    type=ProjectGroupMembership.MembershipType.MEMBER
                )
                self.stats["memberships"] += 1

        self.log_count("Groups created", self.stats["groups"])
        self.log_count("Memberships created", self.stats["memberships"])

    def import_jobs(self, records):
        """Import Jobs -> Job."""
        # First pass: create jobs without parent references
        job_parent_map = {}  # new_job_id -> legacy_parent_jobid

        for record in records:
            try:
                pk = record.get("pk")
                fields = record.get("fields", {})
                project_id = fields.get("projectid")

                if project_id not in self.project_map:
                    self.log(f"Skipping job {pk}: project {project_id} not found")
                    continue

                project = self.project_map[project_id]

                job = Job.objects.create(
                    uuid=self.convert_uuid(pk) or None,
                    project=project,
                    number=fields.get("jobnumber", "1"),
                    title=fields.get("jobtitle", "") or "",
                    status=fields.get("status", 0) or 0,
                    evaluation=fields.get("evaluation", 0) or 0,
                    comments="",
                    creation_time=self.convert_timestamp(fields.get("creationtime")) or timezone.now(),
                    finish_time=self.convert_timestamp(fields.get("finishtime")),
                    task_name=fields.get("taskname", "unknown"),
                    process_id=fields.get("processid"),
                )
                self.job_map[pk] = job
                self.stats["jobs"] += 1

                parent_jobid = fields.get("parentjobid")
                if parent_jobid:
                    job_parent_map[job.id] = parent_jobid

            except Exception as e:
                if self.continue_on_error:
                    self.stderr.write(self.style.ERROR(f"  Error importing job {pk}: {e}"))
                else:
                    raise

        # Second pass: set parent references
        for job_id, legacy_parent_id in job_parent_map.items():
            if legacy_parent_id in self.job_map:
                job = Job.objects.get(id=job_id)
                job.parent = self.job_map[legacy_parent_id]
                job.save()

        self.log_count("Jobs", self.stats["jobs"])

    def import_serverjobs(self, records):
        """Import Serverjobs -> ServerJob."""
        for record in records:
            try:
                pk = record.get("pk")
                fields = record.get("fields", {})
                job_id = fields.get("jobid") or pk

                if job_id not in self.job_map:
                    continue

                job = self.job_map[job_id]

                ServerJob.objects.create(
                    job=job,
                    server_process_id=fields.get("serverprocessid"),
                    machine=fields.get("machine", "") or "",
                    username=fields.get("username", "") or "",
                    mechanism=fields.get("mechanism", "") or "",
                    remote_path=fields.get("remotepath", "") or "",
                    custom_code_file=fields.get("customcodefile", "") or "",
                    validate=fields.get("validate", "") or "",
                    key_file_name=fields.get("keyfilename", "") or "",
                    server_group=fields.get("servergroup", "") or "",
                )
                self.stats["serverjobs"] += 1

            except Exception as e:
                if self.continue_on_error:
                    self.stderr.write(self.style.ERROR(f"  Error importing serverjob {pk}: {e}"))
                else:
                    raise

        self.log_count("ServerJobs", self.stats["serverjobs"])

    def import_jobkeyvalues(self, records):
        """Import Jobkeyvalues -> JobFloatValue."""
        for record in records:
            try:
                fields = record.get("fields", {})
                job_id = fields.get("jobid")
                keytype_id = fields.get("keytypeid")
                value = fields.get("value")

                if job_id not in self.job_map or keytype_id not in self.keytype_map:
                    continue

                JobFloatValue.objects.create(
                    job=self.job_map[job_id],
                    key=self.keytype_map[keytype_id],
                    value=float(value) if value is not None else 0.0,
                )
                self.stats["jobfloatvalues"] += 1

            except Exception as e:
                if self.continue_on_error:
                    pass  # Silently skip to avoid too much noise
                else:
                    raise

        self.log_count("JobFloatValues", self.stats["jobfloatvalues"])

    def import_jobkeycharvalues(self, records):
        """Import Jobkeycharvalues -> JobCharValue."""
        for record in records:
            try:
                fields = record.get("fields", {})
                job_id = fields.get("jobid")
                keytype_id = fields.get("keytypeid")
                value = fields.get("value", "")

                if job_id not in self.job_map or keytype_id not in self.keytype_map:
                    continue

                JobCharValue.objects.create(
                    job=self.job_map[job_id],
                    key=self.keytype_map[keytype_id],
                    value=value or "",
                )
                self.stats["jobcharvalues"] += 1

            except Exception as e:
                if self.continue_on_error:
                    pass
                else:
                    raise

        self.log_count("JobCharValues", self.stats["jobcharvalues"])

    def import_xdata(self, records):
        """Import Xdata -> XData."""
        for record in records:
            try:
                fields = record.get("fields", {})
                job_id = fields.get("jobid")

                job = self.job_map.get(job_id) if job_id else None

                XData.objects.create(
                    data_class=fields.get("xdataclass", ""),
                    xml=fields.get("xdataxml", ""),
                    job=job,
                )
                self.stats["xdata"] += 1

            except Exception as e:
                if self.continue_on_error:
                    pass
                else:
                    raise

        self.log_count("XData", self.stats["xdata"])

    def import_files(self, records):
        """Import Files -> File."""
        for record in records:
            try:
                pk = record.get("pk")
                fields = record.get("fields", {})
                job_id = fields.get("jobid")
                filetype_id = fields.get("filetypeid")

                if job_id not in self.job_map:
                    self.log(f"Skipping file {pk}: job {job_id} not found")
                    continue

                if filetype_id not in self.filetype_map:
                    # Create a default filetype
                    filetype, _ = FileType.objects.get_or_create(
                        name=f"unknown_{filetype_id}",
                        defaults={"description": "Unknown file type from legacy import"}
                    )
                    self.filetype_map[filetype_id] = filetype

                # Map pathflag to Directory enum
                pathflag = fields.get("pathflag", 1)
                if pathflag == 2:
                    directory = File.Directory.IMPORT_DIR
                else:
                    directory = File.Directory.JOB_DIR

                file_obj = File.objects.create(
                    uuid=self.convert_uuid(pk) or None,
                    name=fields.get("filename", ""),
                    directory=directory,
                    type=self.filetype_map[filetype_id],
                    sub_type=fields.get("filesubtype"),
                    content=fields.get("filecontent"),
                    annotation=fields.get("annotation", "") or "",
                    job=self.job_map[job_id],
                    job_param_name=fields.get("jobparamname", "") or "",
                )
                self.file_map[pk] = file_obj
                self.stats["files"] += 1

            except Exception as e:
                if self.continue_on_error:
                    self.stderr.write(self.style.ERROR(f"  Error importing file {pk}: {e}"))
                else:
                    raise

        self.log_count("Files", self.stats["files"])

    def import_importfiles(self, records):
        """Import Importfiles -> FileImport."""
        for record in records:
            try:
                fields = record.get("fields", {})
                file_id = fields.get("fileid")

                if file_id not in self.file_map:
                    continue

                FileImport.objects.create(
                    file=self.file_map[file_id],
                    time=self.convert_timestamp(fields.get("creationtime")) or timezone.now(),
                    name=fields.get("sourcefilename", "") or "",
                    checksum=fields.get("checksum", "") or "",
                    last_modified=self.convert_timestamp(fields.get("lastmodifiedtime")),
                )
                self.stats["fileimports"] += 1

            except Exception as e:
                if self.continue_on_error:
                    pass
                else:
                    raise

        self.log_count("FileImports", self.stats["fileimports"])

    def import_exportfiles(self, records):
        """Import Exportfiles -> FileExport."""
        for record in records:
            try:
                fields = record.get("fields", {})
                file_id = fields.get("fileid")

                if file_id not in self.file_map:
                    continue

                FileExport.objects.create(
                    file=self.file_map[file_id],
                    time=self.convert_timestamp(fields.get("creationtime")) or timezone.now(),
                    name=fields.get("exportfilename", "") or "",
                )
                self.stats["fileexports"] += 1

            except Exception as e:
                if self.continue_on_error:
                    pass
                else:
                    raise

        self.log_count("FileExports", self.stats["fileexports"])

    def import_fileuses(self, records):
        """Import Fileuses -> FileUse."""
        for record in records:
            try:
                fields = record.get("fields", {})
                file_id = fields.get("fileid")
                job_id = fields.get("jobid")
                role_id = fields.get("roleid")

                if file_id not in self.file_map or job_id not in self.job_map:
                    continue

                role = self.filerole_map.get(role_id, FileUse.Role.OUT)

                FileUse.objects.create(
                    file=self.file_map[file_id],
                    job=self.job_map[job_id],
                    role=role,
                    job_param_name=fields.get("jobparamname", "") or "",
                )
                self.stats["fileuses"] += 1

            except Exception as e:
                if self.continue_on_error:
                    pass
                else:
                    raise

        self.log_count("FileUses", self.stats["fileuses"])

    def import_projectexports(self, records):
        """Import Projectexports -> ProjectExport."""
        for record in records:
            try:
                fields = record.get("fields", {})
                project_id = fields.get("projectid")

                if project_id not in self.project_map:
                    continue

                ProjectExport.objects.create(
                    project=self.project_map[project_id],
                    time=self.convert_timestamp(fields.get("projectexporttime")) or timezone.now(),
                )
                self.stats["projectexports"] += 1

            except Exception as e:
                if self.continue_on_error:
                    pass
                else:
                    raise

        self.log_count("ProjectExports", self.stats["projectexports"])

    def import_projectimports(self, records):
        """Import Projectimports -> ProjectImport."""
        for record in records:
            try:
                fields = record.get("fields", {})
                project_id = fields.get("projectid")

                if project_id not in self.project_map:
                    continue

                ProjectImport.objects.create(
                    project=self.project_map[project_id],
                    time=self.convert_timestamp(fields.get("projectimporttime")) or timezone.now(),
                )
                self.stats["projectimports"] += 1

            except Exception as e:
                if self.continue_on_error:
                    pass
                else:
                    raise

        self.log_count("ProjectImports", self.stats["projectimports"])

    def print_summary(self):
        """Print import summary."""
        self.stdout.write("\n" + "-" * 60)

        if self.dry_run:
            self.stdout.write(self.style.WARNING("DRY RUN completed (no changes saved)"))
        else:
            self.stdout.write(self.style.SUCCESS("Import completed successfully!"))

        total = sum(self.stats.values())
        self.stdout.write(f"  Total records processed: {total}")

        if self.stats["projects_renamed"]:
            self.stdout.write(f"  Projects renamed due to collision: {self.stats['projects_renamed']}")

        self.stdout.write("-" * 60 + "\n")


class DryRunRollback(Exception):
    """Exception to trigger rollback in dry-run mode."""
    pass
