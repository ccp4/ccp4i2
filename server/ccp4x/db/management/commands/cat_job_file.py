"""
Django management command to view files within a job directory.

Usage:
    python manage.py cat_job_file <project> <job> <relative_path>
    python manage.py cat_job_file toxd 1 input_params.xml
    python manage.py cat_job_file toxd 1 program.log --tail 50
    python manage.py cat_job_file toxd 1 output.mtz --binary
"""

import uuid
from pathlib import Path
from django.core.management.base import BaseCommand, CommandError
from ccp4x.db.models import Job, Project


class Command(BaseCommand):
    """View file contents within a job directory."""

    help = "Display file contents relative to a job directory"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument(
            'project',
            help='Project name or UUID'
        )
        parser.add_argument(
            'job',
            help='Job number or UUID'
        )
        parser.add_argument(
            'relative_path',
            help='Relative path to file within job directory'
        )
        parser.add_argument(
            '--lines',
            type=int,
            help='Number of lines to show from start (like head -n)'
        )
        parser.add_argument(
            '--tail',
            type=int,
            help='Number of lines to show from end (like tail -n)'
        )
        parser.add_argument(
            '--binary',
            action='store_true',
            help='Allow viewing binary files (show hex dump)'
        )
        parser.add_argument(
            '--no-line-numbers',
            action='store_true',
            help='Do not show line numbers'
        )

    def handle(self, *args, **options):
        project_identifier = options['project']
        job_identifier = options['job']
        relative_path = options['relative_path']

        # Find job
        try:
            job = self.get_job(project_identifier, job_identifier)
        except (Job.DoesNotExist, Project.DoesNotExist) as e:
            raise CommandError(str(e))

        # Construct full path
        job_dir = Path(job.directory)
        file_path = job_dir / relative_path

        # Security check: ensure file is within job directory
        try:
            file_path = file_path.resolve()
            job_dir = job_dir.resolve()
            if not str(file_path).startswith(str(job_dir)):
                raise CommandError(
                    f"Security error: Path '{relative_path}' is outside job directory"
                )
        except Exception as e:
            raise CommandError(f"Invalid path: {e}")

        # Check file exists
        if not file_path.exists():
            raise CommandError(
                f"File not found: {relative_path}\n"
                f"Full path: {file_path}\n"
                f"Use 'ccp4i2 jobs tree {job.project.name} {job.number}' to see available files."
            )

        if not file_path.is_file():
            raise CommandError(
                f"Not a file: {relative_path} (it's a directory)\n"
                f"Use 'ccp4i2 jobs tree {job.project.name} {job.number}' to see directory contents."
            )

        # Display file
        self.display_file(
            file_path,
            relative_path,
            job,
            options['lines'],
            options['tail'],
            options['binary'],
            not options['no_line_numbers']
        )

    def display_file(self, file_path: Path, relative_path: str, job,
                     head_lines: int, tail_lines: int, allow_binary: bool, show_line_numbers: bool):
        """Display file contents with formatting."""

        # Get file info
        file_size = file_path.stat().st_size

        # Check if binary
        is_binary = self.is_binary_file(file_path)

        if is_binary and not allow_binary:
            raise CommandError(
                f"File appears to be binary: {relative_path}\n"
                f"Use --binary flag to view as hex dump, or use an appropriate viewer.\n"
                f"File size: {self.format_size(file_size)}"
            )

        # Print header
        self.stdout.write(self.style.SUCCESS(f"\n{'='*80}"))
        self.stdout.write(self.style.SUCCESS(f"Project: {job.project.name}"))
        self.stdout.write(self.style.SUCCESS(f"Job: #{job.number} - {job.title} ({job.task_name})"))
        self.stdout.write(self.style.SUCCESS(f"File: {relative_path}"))
        self.stdout.write(self.style.SUCCESS(f"Size: {self.format_size(file_size)}"))
        self.stdout.write(self.style.SUCCESS(f"{'='*80}\n"))

        try:
            if is_binary:
                self.display_binary_file(file_path, head_lines, tail_lines)
            else:
                self.display_text_file(file_path, head_lines, tail_lines, show_line_numbers)
        except Exception as e:
            raise CommandError(f"Error reading file: {e}")

        # Print footer
        self.stdout.write(self.style.SUCCESS(f"\n{'='*80}"))
        self.stdout.write(self.style.SUCCESS(f"End of file: {relative_path}"))
        self.stdout.write(self.style.SUCCESS(f"{'='*80}\n"))

    def display_text_file(self, file_path: Path, head_lines: int, tail_lines: int, show_line_numbers: bool):
        """Display text file contents."""
        with open(file_path, 'r', encoding='utf-8', errors='replace') as f:
            lines = f.readlines()

        total_lines = len(lines)

        # Apply head/tail filtering
        if head_lines is not None:
            lines = lines[:head_lines]
            if total_lines > head_lines:
                truncated_msg = f"... (showing first {head_lines} of {total_lines} lines)"
        elif tail_lines is not None:
            lines = lines[-tail_lines:]
            if total_lines > tail_lines:
                truncated_msg = f"... (showing last {tail_lines} of {total_lines} lines)"
        else:
            truncated_msg = None

        # Display lines
        if show_line_numbers:
            # Calculate width for line numbers
            max_line_num = total_lines if tail_lines is None else total_lines
            width = len(str(max_line_num))

            start_line = 1 if tail_lines is None else (total_lines - len(lines) + 1)

            for idx, line in enumerate(lines, start=start_line):
                # Remove trailing newline for display
                line_content = line.rstrip('\n')
                self.stdout.write(f"{idx:>{width}} | {line_content}")
        else:
            for line in lines:
                self.stdout.write(line.rstrip('\n'))

        # Show truncation message
        if truncated_msg:
            self.stdout.write(f"\n{truncated_msg}")

    def display_binary_file(self, file_path: Path, head_bytes: int, tail_bytes: int):
        """Display binary file as hex dump."""
        with open(file_path, 'rb') as f:
            data = f.read()

        # Apply head/tail filtering
        total_bytes = len(data)
        if head_bytes is not None:
            data = data[:head_bytes]
            truncated = total_bytes > head_bytes
        elif tail_bytes is not None:
            data = data[-tail_bytes:]
            truncated = total_bytes > tail_bytes
        else:
            # Limit binary display to first 1KB by default
            if len(data) > 1024:
                data = data[:1024]
                truncated = True
            else:
                truncated = False

        # Display hex dump
        self.stdout.write("Hex dump (offset | hex | ASCII):\n")
        for i in range(0, len(data), 16):
            chunk = data[i:i+16]
            hex_part = ' '.join(f'{b:02x}' for b in chunk)
            ascii_part = ''.join(chr(b) if 32 <= b < 127 else '.' for b in chunk)
            self.stdout.write(f"{i:08x} | {hex_part:<48} | {ascii_part}")

        if truncated:
            self.stdout.write(f"\n... (binary file truncated, showing first {len(data)} bytes of {total_bytes})")

    def is_binary_file(self, file_path: Path) -> bool:
        """Check if file is binary by reading first chunk."""
        try:
            with open(file_path, 'rb') as f:
                chunk = f.read(8192)
                if b'\0' in chunk:
                    return True
                # Check if mostly text characters
                text_chars = sum(1 for b in chunk if 32 <= b < 127 or b in (9, 10, 13))
                return text_chars / len(chunk) < 0.7 if chunk else False
        except Exception:
            return True

    def format_size(self, size_bytes: int) -> str:
        """Format file size in human-readable format."""
        for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
            if size_bytes < 1024.0:
                if unit == 'B':
                    return f"{size_bytes} {unit}"
                else:
                    return f"{size_bytes:.1f} {unit}"
            size_bytes /= 1024.0
        return f"{size_bytes:.1f} PB"

    def get_job(self, project_identifier: str, job_identifier: str) -> Job:
        """Get job by project and job identifiers."""
        # First, get the project
        project = None

        # Try project UUID
        try:
            uuid.UUID(project_identifier)
            project = Project.objects.get(uuid=project_identifier)
        except (ValueError, Project.DoesNotExist):
            # Try project name
            try:
                project = Project.objects.get(name=project_identifier)
            except Project.DoesNotExist:
                raise Project.DoesNotExist(
                    f"Project not found: '{project_identifier}'"
                )

        # Now get the job
        # Try job UUID
        try:
            uuid.UUID(job_identifier)
            job = Job.objects.get(uuid=job_identifier, project=project)
            return job
        except (ValueError, Job.DoesNotExist):
            pass

        # Try job number
        try:
            job = Job.objects.get(number=job_identifier, project=project)
            return job
        except Job.DoesNotExist:
            raise Job.DoesNotExist(
                f"Job not found: '{job_identifier}' in project '{project.name}'"
            )
