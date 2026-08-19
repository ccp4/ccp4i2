"""The project export must not issue queries in proportion to project size.

Exporting walked ``file.job``, ``file.type`` and ``file_use.file`` row by row,
which cost several queries per file: a 250-job project ran 3757 queries and took
over half a second. That is slow for an export and prohibitive for anything
that wants to snapshot a project routinely.

These tests pin the shape rather than a magic number -- doubling the project
must not change the query count -- so the regression they guard against is
caught however the export is rewritten.
"""

from xml.etree import ElementTree as ET

from django.db import connection
from django.test import TestCase
from django.test.utils import CaptureQueriesContext

from ...db.export_project import generate_project_xml_tree
from ...db.models import (
    File,
    FileImport,
    FileType,
    FileUse,
    Job,
    JobCharValue,
    JobFloatValue,
    JobValueKey,
    Project,
)


def build_project(name: str, n_jobs: int, files_per_job: int = 3) -> Project:
    """A project of a given size, with the relations the export walks."""
    file_type, _ = FileType.objects.get_or_create(
        name="chemical/x-pdb", defaults={"description": "Model coordinates"}
    )
    project = Project.objects.create(name=name, directory=f"/tmp/{name}")

    Job.objects.bulk_create(
        [
            Job(
                project=project,
                number=str(i + 1),
                task_name="prosmart_refmac",
                title=f"Refinement {i + 1}",
                status=Job.Status.FINISHED,
            )
            for i in range(n_jobs)
        ]
    )
    jobs = list(Job.objects.filter(project=project))

    File.objects.bulk_create(
        [
            File(
                name=f"out{k}.pdb",
                # Alternate so both branches of the path handling are covered.
                directory=(
                    File.Directory.JOB_DIR if k % 2 == 0 else File.Directory.IMPORT_DIR
                ),
                type=file_type,
                job=job,
                annotation=f"Output {k} of job {job.number}",
                sub_type=k,
                content=k,
                job_param_name=f"XYZOUT{k}",
            )
            for job in jobs
            for k in range(files_per_job)
        ]
    )
    files = list(File.objects.filter(job__project=project).select_related("job"))

    FileUse.objects.bulk_create(
        [
            FileUse(
                file=file, job=file.job, role=FileUse.Role.OUT, job_param_name="XYZOUT"
            )
            for file in files
        ]
    )
    FileImport.objects.bulk_create(
        [
            FileImport(file=file, name="/somewhere/original.pdb", checksum="deadbeef")
            for file in files
            if file.directory == File.Directory.IMPORT_DIR
        ]
    )

    r_free, _ = JobValueKey.objects.get_or_create(
        name="RFree", defaults={"description": "Free R factor"}
    )
    space_group, _ = JobValueKey.objects.get_or_create(
        name="SpaceGroup", defaults={"description": "Space group"}
    )
    JobFloatValue.objects.bulk_create(
        [JobFloatValue(job=job, key=r_free, value=0.21) for job in jobs]
    )
    JobCharValue.objects.bulk_create(
        [JobCharValue(job=job, key=space_group, value="P 21 21 21") for job in jobs]
    )
    return project


def count_queries(project: Project) -> int:
    with CaptureQueriesContext(connection) as captured:
        generate_project_xml_tree(project)
    return len(captured.captured_queries)


class ExportQueryCountTest(TestCase):
    def test_query_count_does_not_grow_with_the_project(self):
        small = build_project("Small", n_jobs=5)
        large = build_project("Large", n_jobs=40)

        self.assertEqual(count_queries(small), count_queries(large))

    def test_query_count_is_a_handful(self):
        """A per-table query each, not a per-row one."""
        project = build_project("Modest", n_jobs=20)
        self.assertLess(count_queries(project), 20)

    def test_every_row_still_appears(self):
        project = build_project("Complete", n_jobs=7, files_per_job=4)
        root = generate_project_xml_tree(project)

        self.assertEqual(len(root.findall(".//jobTable/job")), 7)
        self.assertEqual(len(root.findall(".//fileTable/file")), 28)
        self.assertEqual(len(root.findall(".//fileuseTable/fileuse")), 28)
        self.assertEqual(len(root.findall(".//importfileTable/importfile")), 14)
        self.assertEqual(len(root.findall(".//jobkeyvalueTable/jobkeyvalue")), 7)
        self.assertEqual(len(root.findall(".//jobkeyvalueTable/jobkeycharvalue")), 7)

    def test_file_attributes_survive(self):
        project = build_project("Attributes", n_jobs=1, files_per_job=2)
        root = generate_project_xml_tree(project)
        files = {
            element.get("filename"): element
            for element in root.findall(".//fileTable/file")
        }

        self.assertEqual(files["out0.pdb"].get("pathflag"), str(File.Directory.JOB_DIR))
        self.assertEqual(
            files["out1.pdb"].get("pathflag"), str(File.Directory.IMPORT_DIR)
        )
        self.assertEqual(files["out0.pdb"].get("annotation"), "Output 0 of job 1")
        self.assertEqual(files["out1.pdb"].get("filesubtype"), "1")
        self.assertEqual(files["out1.pdb"].get("filecontent"), "1")
        self.assertEqual(files["out0.pdb"].get("jobparamname"), "XYZOUT0")

    def test_export_survives_a_file_with_no_job(self):
        """``pathflag`` used to be derived by building the file's full path,
        which reaches through ``job`` -- and ``File.job`` is nullable."""
        project = build_project("Orphan", n_jobs=1, files_per_job=1)
        file_type = FileType.objects.get(name="chemical/x-pdb")
        File.objects.create(
            name="orphan.pdb",
            directory=File.Directory.JOB_DIR,
            type=file_type,
            job=None,
        )

        root = generate_project_xml_tree(project)
        self.assertEqual(len(root.findall(".//fileTable/file")), 1)
