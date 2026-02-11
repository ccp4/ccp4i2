import faulthandler
import uuid
import os
import subprocess
import platform
from django.core.management.base import BaseCommand
from asgiref.sync import async_to_sync
from ccp4i2.db.models import Job, Project
from ccp4i2.lib.async_run_job import run_job_async


class Command(BaseCommand):
    """
    A Django management command to run a job.

    If the job is not in PENDING status, prompts for confirmation before running
    (unless -y/--yes is specified).

    Usage:
        python manage.py run_job --projectname toxd --jobnumber 5
        python manage.py run_job --projectname toxd --jobnumber 5 -y  # Skip confirmation
        python manage.py run_job --jobid 42 --detach
    """

    help = "Run a job"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument("-ji", "--jobid", help="Integer job id", type=int)
        parser.add_argument("-ju", "--jobuuid", help="Job UUID", type=str)
        parser.add_argument("-pn", "--projectname", help="Project name", type=str)
        parser.add_argument("-pi", "--projectid", help="Integer project id", type=int)
        parser.add_argument("-pu", "--projectuuid", help="Project uuid", type=str)
        parser.add_argument("-jn", "--jobnumber", help="Job number", type=str)
        parser.add_argument("-d", "--detach", help="Detach job", action="store_true")
        parser.add_argument(
            "-y", "--yes",
            help="Skip confirmation prompt for non-pending jobs",
            action="store_true"
        )

    def handle(self, *args, **options):
        try:
            the_job = self.get_job(options)
        except (Job.DoesNotExist, Project.DoesNotExist) as e:
            self.stderr.write(self.style.ERROR(str(e)))
            return

        # Check if job is not in PENDING status
        if the_job.status not in [Job.Status.PENDING, Job.Status.QUEUED]:
            status_label = Job.Status(the_job.status).label
            if not options.get("yes"):
                self.stdout.write(
                    self.style.WARNING(
                        f"Job #{the_job.number} '{the_job.title}' has status '{status_label}' (not Pending)."
                    )
                )
                self.stdout.write(
                    f"Re-running this job will overwrite existing results."
                )
                try:
                    response = input("Continue? [y/N] ").strip().lower()
                    if response not in ('y', 'yes'):
                        self.stdout.write("Aborted.")
                        return
                except (EOFError, KeyboardInterrupt):
                    self.stdout.write("\nAborted.")
                    return

        if options["detach"]:
            # Determine the program name based on the OS
            ccp4_python_program = "ccp4-python"
            if platform.system() == "Windows":
                ccp4_python_program += ".bat"

            # Open file for capturing stdout
            with open(
                the_job.directory / "cplusplus_stdout.txt", "w", encoding="utf-8"
            ) as stdout_file:
                process = subprocess.Popen(
                    [
                        ccp4_python_program,
                        "manage.py",
                        "run_job",
                        "-ju",
                        f"{str(the_job.uuid)}",
                    ],
                    start_new_session=True,
                    stdout=stdout_file,  # Capture stdout
                    stderr=subprocess.STDOUT,  # Merge stderr into stdout
                )
                the_job.process_id = process.pid
                the_job.save()
        else:
            # Enable faulthandler to produce a Python traceback on SIGSEGV/SIGABRT.
            # Without this, a segfault in a C extension (e.g. iris_validation/mmdb2)
            # kills the process silently with no diagnostic output.
            crash_log_path = the_job.directory / "crash_traceback.txt"
            self._crash_log = open(crash_log_path, 'w')
            faulthandler.enable(file=self._crash_log)

            with open(
                the_job.directory / "cplusplus_stdout.txt", "w", encoding="utf-8"
            ) as stdout_file:
                # Redirect file descriptors
                stdout_fd = stdout_file.fileno()
                os.dup2(stdout_fd, 1)  # Redirect stdout
                os.dup2(stdout_fd, 2)  # Redirect stderr
                # Use modern async run_job with proper async/sync bridging
                async_to_sync(run_job_async)(the_job.uuid)

    def get_job(self, options):
        if options["jobid"] is not None:
            return Job.objects.get(id=options["jobid"])
        if options["jobuuid"] is not None:
            return Job.objects.get(uuid=uuid.UUID(options["jobuuid"]))
        if options["projectname"] is not None and options["jobnumber"] is not None:
            the_project = Project.objects.get(name=options["projectname"])
            return Job.objects.get(number=options["jobnumber"], project=the_project)
        if options["projectid"] is not None and options["jobnumber"] is not None:
            the_project = Project.objects.get(id=options["projectid"])
            return Job.objects.get(number=options["jobnumber"], project=the_project)
        if options["projectuuid"] is not None and options["jobnumber"] is not None:
            the_project = Project.objects.get(uuid=uuid.UUID(options["projectuuid"]))
            return Job.objects.get(number=options["jobnumber"], project=the_project)
        raise Job.DoesNotExist("No job found with the provided criteria.")
