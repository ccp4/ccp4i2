"""
Django management command to list jobs for a CCP4i2 project.

Usage:
    python manage.py list_jobs <project_name_or_uuid>
    python manage.py list_jobs toxd --json
    python manage.py list_jobs toxd --status FINISHED
    python manage.py list_jobs toxd --sort-by creation_time --reverse
    python manage.py list_jobs toxd --filter "demo*"

Example:
    python manage.py list_jobs toxd
    python manage.py list_jobs toxd --status RUNNING --json
    python manage.py list_jobs toxd --sort-by number
    python manage.py list_jobs --all  # List all jobs from all projects
"""

from django.core.management.base import BaseCommand, CommandError
from ccp4x.db import models
import fnmatch


class Command(BaseCommand):
    """
    Django management command to list jobs for a CCP4i2 project.

    Displays job information including:
    - Job number
    - Title
    - Task name
    - Status
    - Evaluation
    - Creation time
    - Finish time
    - Process ID

    Supports multiple output formats and sorting/filtering options.

    Attributes:
        help (str): Description of the command displayed with --help
        requires_system_checks (list): System checks required before running
    """

    help = "List jobs for a CCP4i2 project"
    requires_system_checks = []

    def add_arguments(self, parser):
        """
        Add command-line arguments.

        Args:
            parser: Django argument parser
        """
        # Positional argument for project
        parser.add_argument(
            "project",
            nargs="?",
            type=str,
            help="Project name or UUID (required unless --all is specified)",
        )

        # Option to list all jobs from all projects
        parser.add_argument(
            "--all",
            action="store_true",
            help="List all jobs from all projects",
        )

        # Output format options
        parser.add_argument(
            "--json",
            action="store_true",
            help="Output as JSON array",
        )

        parser.add_argument(
            "--ids-only",
            action="store_true",
            help="Output only job UUIDs (one per line)",
        )

        parser.add_argument(
            "--numbers-only",
            action="store_true",
            help="Output only job numbers (one per line)",
        )

        # Sorting options
        parser.add_argument(
            "--sort-by",
            type=str,
            choices=["number", "title", "task_name", "status", "creation_time", "finish_time"],
            default="number",
            help="Sort jobs by field (default: number)",
        )

        parser.add_argument(
            "--reverse",
            action="store_true",
            help="Reverse sort order",
        )

        # Filtering options
        parser.add_argument(
            "--filter",
            type=str,
            help="Filter jobs by title pattern (supports wildcards like 'demo*')",
        )

        parser.add_argument(
            "--status",
            type=str,
            choices=[
                "UNKNOWN", "PENDING", "QUEUED", "RUNNING", "INTERRUPTED",
                "FAILED", "FINISHED", "RUNNING_REMOTELY", "FILE_HOLDER",
                "TO_DELETE", "UNSATISFACTORY"
            ],
            help="Filter jobs by status",
        )

        parser.add_argument(
            "--evaluation",
            type=str,
            choices=["UNKNOWN", "BEST", "GOOD", "REJECTED"],
            help="Filter jobs by evaluation",
        )

        parser.add_argument(
            "--task",
            type=str,
            help="Filter jobs by task name",
        )

        parser.add_argument(
            "--limit",
            type=int,
            help="Limit number of results",
        )

    def handle(self, *args, **options):
        """
        Handle the list_jobs command.

        Args:
            *args: Positional arguments (unused)
            **options: Command options

        Raises:
            CommandError: If listing fails or project not found
        """
        import json

        try:
            # Determine which jobs to query
            list_all = options.get("all")
            project_identifier = options.get("project")

            if not list_all and not project_identifier:
                raise CommandError(
                    "Either specify a project name/UUID or use --all to list all jobs"
                )

            if list_all and project_identifier:
                raise CommandError(
                    "Cannot specify both a project and --all option"
                )

            # Get jobs
            if list_all:
                jobs = models.Job.objects.select_related('project').all()
                project_name = "All Projects"
            else:
                # Find the project by name or UUID
                project = None
                try:
                    # Try as UUID first
                    from uuid import UUID
                    UUID(project_identifier)  # Validate UUID format
                    project = models.Project.objects.get(uuid=project_identifier)
                except (models.Project.DoesNotExist, ValueError, TypeError):
                    # Not a valid UUID or not found, try as name
                    try:
                        project = models.Project.objects.get(name=project_identifier)
                    except models.Project.DoesNotExist:
                        raise CommandError(
                            f"Project not found: '{project_identifier}'. "
                            "Use 'python manage.py list_projects' to see available projects."
                        )

                project_name = project.name
                jobs = models.Job.objects.filter(project=project)

            # Apply filters
            status_filter = options.get("status")
            if status_filter:
                status_value = getattr(models.Job.Status, status_filter)
                jobs = jobs.filter(status=status_value)

            evaluation_filter = options.get("evaluation")
            if evaluation_filter:
                eval_value = getattr(models.Job.Evaluation, evaluation_filter)
                jobs = jobs.filter(evaluation=eval_value)

            task_filter = options.get("task")
            if task_filter:
                jobs = jobs.filter(task_name=task_filter)

            # Convert to list for filtering and sorting
            jobs = list(jobs)

            # Apply title filter if specified
            title_filter = options.get("filter")
            if title_filter:
                jobs = [j for j in jobs if fnmatch.fnmatch(j.title, title_filter)]

            # Sort jobs
            sort_by = options["sort_by"]
            reverse = options["reverse"]

            if sort_by == "number":
                # Natural sort for job numbers (1, 2, 10 instead of 1, 10, 2)
                def natural_sort_key(job):
                    parts = job.number.split('.')
                    return [int(p) if p.isdigit() else p for p in parts]
                jobs.sort(key=natural_sort_key, reverse=reverse)
            elif sort_by == "title":
                jobs.sort(key=lambda j: j.title, reverse=reverse)
            elif sort_by == "task_name":
                jobs.sort(key=lambda j: j.task_name, reverse=reverse)
            elif sort_by == "status":
                jobs.sort(key=lambda j: j.status, reverse=reverse)
            elif sort_by == "creation_time":
                jobs.sort(key=lambda j: j.creation_time, reverse=reverse)
            elif sort_by == "finish_time":
                jobs.sort(key=lambda j: j.finish_time or "", reverse=reverse)

            # Apply limit if specified
            limit = options.get("limit")
            if limit:
                jobs = jobs[:limit]

            # Handle empty results
            if not jobs:
                if title_filter or status_filter or evaluation_filter or task_filter:
                    self.stdout.write(f"No jobs found matching the filters in project '{project_name}'")
                else:
                    self.stdout.write(f"No jobs found in project '{project_name}'")
                return

            # Output based on format
            json_output = options["json"]
            ids_only = options["ids_only"]
            numbers_only = options["numbers_only"]

            if json_output:
                # JSON output
                output = []
                for job in jobs:
                    job_data = {
                        "uuid": str(job.uuid),
                        "number": job.number,
                        "title": job.title,
                        "task_name": job.task_name,
                        "status": models.Job.Status(job.status).label,
                        "status_code": job.status,
                        "evaluation": models.Job.Evaluation(job.evaluation).label,
                        "evaluation_code": job.evaluation,
                        "creation_time": job.creation_time.isoformat() if job.creation_time else None,
                        "finish_time": job.finish_time.isoformat() if job.finish_time else None,
                        "process_id": job.process_id,
                        "comments": job.comments,
                    }
                    if list_all:
                        job_data["project_name"] = job.project.name
                        job_data["project_uuid"] = str(job.project.uuid)
                    output.append(job_data)
                self.stdout.write(json.dumps(output, indent=2))

            elif ids_only:
                # UUIDs only output (one per line)
                for job in jobs:
                    self.stdout.write(str(job.uuid))

            elif numbers_only:
                # Job numbers only output (one per line)
                for job in jobs:
                    self.stdout.write(job.number)

            else:
                # Verbose table output
                self.stdout.write(self.style.SUCCESS(f"\n{'='*140}"))
                if list_all:
                    self.stdout.write(self.style.SUCCESS(f"Jobs from All Projects ({len(jobs)} found)"))
                else:
                    self.stdout.write(self.style.SUCCESS(f"Jobs in Project: {project_name} ({len(jobs)} found)"))
                self.stdout.write(self.style.SUCCESS(f"{'='*140}\n"))

                # Table header
                if list_all:
                    self.stdout.write(
                        f"{'#':<6} {'Title':<25} {'Task':<20} {'Status':<15} {'Eval':<8} {'Project':<20} {'Created':<20}"
                    )
                else:
                    self.stdout.write(
                        f"{'#':<6} {'Title':<30} {'Task':<25} {'Status':<15} {'Eval':<8} {'Created':<20}"
                    )
                self.stdout.write("-" * 140)

                # Table rows
                for job in jobs:
                    # Format fields
                    number = job.number[:5] if len(job.number) > 5 else job.number

                    if list_all:
                        title = job.title[:24] if len(job.title) > 24 else job.title
                        task = job.task_name[:19] if len(job.task_name) > 19 else job.task_name
                        project = job.project.name[:19] if len(job.project.name) > 19 else job.project.name
                    else:
                        title = job.title[:29] if len(job.title) > 29 else job.title
                        task = job.task_name[:24] if len(job.task_name) > 24 else job.task_name

                    status = models.Job.Status(job.status).label
                    evaluation = models.Job.Evaluation(job.evaluation).label
                    created = job.creation_time.strftime("%Y-%m-%d %H:%M:%S") if job.creation_time else "Unknown"

                    if list_all:
                        self.stdout.write(
                            f"{number:<6} {title:<25} {task:<20} {status:<15} {evaluation:<8} {project:<20} {created:<20}"
                        )
                    else:
                        self.stdout.write(
                            f"{number:<6} {title:<30} {task:<25} {status:<15} {evaluation:<8} {created:<20}"
                        )

                # Summary footer
                self.stdout.write(f"\n{'='*140}")
                self.stdout.write(f"Total: {len(jobs)} jobs")

                # Status breakdown
                status_counts = {}
                for job in jobs:
                    status_label = models.Job.Status(job.status).label
                    status_counts[status_label] = status_counts.get(status_label, 0) + 1

                if status_counts:
                    status_summary = ", ".join([f"{count} {status}" for status, count in sorted(status_counts.items())])
                    self.stdout.write(f"Status: {status_summary}")

                self.stdout.write(f"{'='*140}\n")

                # Show filter info if applied
                filters_applied = []
                if title_filter:
                    filters_applied.append(f"title={title_filter}")
                if status_filter:
                    filters_applied.append(f"status={status_filter}")
                if evaluation_filter:
                    filters_applied.append(f"evaluation={evaluation_filter}")
                if task_filter:
                    filters_applied.append(f"task={task_filter}")

                if filters_applied:
                    self.stdout.write(f"Filters: {', '.join(filters_applied)}")

                if sort_by != "number" or reverse:
                    sort_desc = f"Sorted by: {sort_by}"
                    if reverse:
                        sort_desc += " (reversed)"
                    self.stdout.write(sort_desc)

        except CommandError:
            raise
        except Exception as e:
            raise CommandError(f"Failed to list jobs: {str(e)}")
