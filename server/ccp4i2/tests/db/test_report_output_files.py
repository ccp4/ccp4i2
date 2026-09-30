"""A job's report lists as output what it made, not what it imported.

A job that imports a file (a model given by path) is recorded as the file's
job, with the import directory; the report's Output Data listed it too, a
second time and under its parameter's name ("XYZIN"), beside the files the
job actually produced.
"""
from django.test import TestCase

from ...db.models import File, FileType, Job, Project
from ...lib.utils.reporting.i2_report import _output_files


class OutputFilesTests(TestCase):
    def test_imported_files_are_not_outputs(self):
        file_type, _ = FileType.objects.get_or_create(
            name="chemical/x-pdb", defaults={"description": "Model coordinates"})
        project = Project.objects.create(name="Links", directory="/tmp/Links")
        job = Job.objects.create(project=project, number="1", task_name="MakeLink",
                                 title="Make link", status=Job.Status.FINISHED)
        File.objects.create(name="6ndn.pdb", directory=File.Directory.IMPORT_DIR,
                            type=file_type, job=job, job_param_name="XYZIN")
        File.objects.create(name="ModelWithLinks.pdb", directory=File.Directory.JOB_DIR,
                            type=file_type, job=job, job_param_name="XYZOUT",
                            annotation="Model with 1 link applied")

        names = [f["filename"] for f in _output_files(job)]
        self.assertEqual(names, ["ModelWithLinks.pdb"])
