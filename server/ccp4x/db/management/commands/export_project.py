import uuid
import os
import subprocess
import platform
from pathlib import Path
from typing import Set
from django.core.management.base import BaseCommand
from ccp4x.db.models import Project, Job
from ccp4x.db.export_project import export_project_to_zip


class Command(BaseCommand):
    """
    A Django management command to export a CCP4 project to a ZIP archive.
    """

    help = "Export a CCP4 project to ZIP archive"
    requires_system_checks = []

    def add_arguments(self, parser):
        # Project identification arguments
        parser.add_argument("-pn", "--projectname", help="Project name", type=str)
        parser.add_argument("-pi", "--projectid", help="Integer project id", type=int)
        parser.add_argument("-pu", "--projectuuid", help="Project uuid", type=str)

        # Job selection arguments
        parser.add_argument(
            "-j",
            "--jobs",
            help="Comma-separated list of top-level job numbers to export (e.g., '1,3,5')",
            type=str,
        )

        # Export options
        parser.add_argument(
            "-o",
            "--output",
            help="Output ZIP file path (default: {project_name}_{timestamp}.zip)",
            type=str,
        )
        parser.add_argument(
            "-d", "--detach", help="Run export in detached process", action="store_true"
        )

    def handle(self, *args, **options):
        try:
            project = self.get_project(options)
        except Project.DoesNotExist as e:
            self.stderr.write(self.style.ERROR(str(e)))
            return

        # Get job selection if specified (as Set[str])
        job_selection = self.get_job_selection(project, options)

        if job_selection is not None:
            # Validate and report on selected jobs
            valid_count = self._count_valid_jobs(project, job_selection)
            self.stdout.write(
                f"Exporting {valid_count} selected top-level jobs and their descendants"
            )

        output_path = self.get_output_path(project, options)

        if options["detach"]:
            self.run_detached_export(project, output_path, options.get("jobs"))
        else:
            self.run_export(project, output_path, job_selection)

    def get_project(self, options):
        """Retrieve project based on provided options."""
        if options["projectname"] is not None:
            return Project.objects.get(name=options["projectname"])
        if options["projectid"] is not None:
            return Project.objects.get(id=options["projectid"])
        if options["projectuuid"] is not None:
            return Project.objects.get(uuid=uuid.UUID(options["projectuuid"]))

        raise Project.DoesNotExist(
            "No project specified. Use --projectname, --projectid, or --projectuuid."
        )

    def get_job_selection(self, project, options) -> Set[str]:
        """Get the set of job numbers to export as strings."""
        if not options.get("jobs"):
            return None

        try:
            # Parse job numbers from comma-separated string and keep as strings
            job_numbers = {num.strip() for num in options["jobs"].split(",")}

            # Validate that the job numbers exist in the project
            valid_jobs = set()
            for job_num_str in job_numbers:
                try:
                    job_num = int(job_num_str)
                    if Job.objects.filter(project=project, number=job_num).exists():
                        valid_jobs.add(job_num_str)
                    else:
                        self.stderr.write(
                            self.style.WARNING(
                                f"Job number {job_num_str} not found in project"
                            )
                        )
                except ValueError:
                    self.stderr.write(
                        self.style.WARNING(f"Invalid job number format: {job_num_str}")
                    )

            if not valid_jobs:
                self.stderr.write(self.style.ERROR("No valid jobs found for selection"))
                return set()

            return valid_jobs

        except ValueError as e:
            self.stderr.write(self.style.ERROR(f"Invalid job numbers format: {e}"))
            return set()

    def _count_valid_jobs(self, project, job_selection: Set[str]) -> int:
        """Count valid jobs for reporting purposes."""
        if not job_selection:
            return 0

        count = 0
        for job_num_str in job_selection:
            try:
                job_num = int(job_num_str)
                if Job.objects.filter(project=project, number=job_num).exists():
                    count += 1
            except ValueError:
                continue
        return count

    def get_output_path(self, project, options):
        """Determine output path for the exported ZIP file."""
        if options["output"]:
            return Path(options["output"])

        # Default output path: {project_name}_{uuid_first_8_chars}.zip
        safe_name = "".join(c for c in project.name if c.isalnum() or c in "._-")
        uuid_prefix = str(project.uuid).replace("-", "")[:8]

        # Add job selection indicator to filename if jobs are selected
        job_suffix = ""
        if options.get("jobs"):
            job_suffix = f"_jobs_{options['jobs'].replace(',', '_')}"

        filename = f"{safe_name}_{uuid_prefix}{job_suffix}.zip"

        return Path.cwd() / filename

    def run_detached_export(self, project, output_path, job_selection_str):
        """Run export in a detached subprocess."""
        # Determine the program name based on the OS
        ccp4_python_program = "ccp4-python"
        if platform.system() == "Windows":
            ccp4_python_program += ".bat"

        # Build command arguments for detached process
        cmd_args = [
            ccp4_python_program,
            "manage.py",
            "export_project",
            "-pu",
            str(project.uuid),
            "-o",
            str(output_path),
        ]

        # Add job selection if specified
        if job_selection_str:
            cmd_args.extend(["-j", job_selection_str])

        # Create log file for detached process
        log_path = output_path.parent / f"{output_path.stem}_export.log"

        try:
            with open(log_path, "w", encoding="utf-8") as log_file:
                process = subprocess.Popen(
                    cmd_args,
                    start_new_session=True,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    cwd=os.getcwd(),
                )

                job_info = f" (Jobs: {job_selection_str})" if job_selection_str else ""
                self.stdout.write(
                    self.style.SUCCESS(
                        f"Export started in detached process (PID: {process.pid})\n"
                        f"Project: {project.name} (UUID: {project.uuid}){job_info}\n"
                        f"Output: {output_path}\n"
                        f"Log file: {log_path}"
                    )
                )

        except Exception as e:
            self.stderr.write(
                self.style.ERROR(f"Failed to start detached export process: {e}")
            )

    def run_export(self, project, output_path, job_selection: Set[str] = None):
        """Run export in the current process."""
        try:
            if job_selection:
                job_info = (
                    f" (Selected job numbers: {', '.join(sorted(job_selection))})"
                )
            else:
                job_info = " (All jobs)"

            self.stdout.write(
                f"Exporting project '{project.name}' (UUID: {project.uuid}){job_info}"
            )
            self.stdout.write(f"Output file: {output_path}")

            # Ensure output directory exists
            output_path.parent.mkdir(parents=True, exist_ok=True)

            # Perform the export with job selection (now passing Set[str])
            result_path = export_project_to_zip(
                project, output_path, job_selection=job_selection
            )

            # Get file size for confirmation
            file_size = result_path.stat().st_size
            file_size_mb = file_size / (1024 * 1024)

            self.stdout.write(
                self.style.SUCCESS(
                    f"Export completed successfully!\n"
                    f"Output file: {result_path}\n"
                    f"File size: {file_size_mb:.2f} MB"
                )
            )

        except Exception as e:
            self.stderr.write(self.style.ERROR(f"Export failed: {e}"))

            # Clean up partial file if it exists
            if output_path.exists():
                try:
                    output_path.unlink()
                    self.stdout.write("Cleaned up partial export file.")
                except Exception:
                    pass
