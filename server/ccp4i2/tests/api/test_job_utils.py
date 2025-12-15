from pathlib import Path
from shutil import rmtree
from typing import List
from django.test import Client
from django.conf import settings
from django.test import TestCase, override_settings
from ccp4i2.core import CCP4Container
from ccp4i2.core.CCP4ModelData import CDictDataFile
from ccp4i2.core.base_object.cdata_file import CDataFile
from ...lib.utils.jobs.i2run import i2run_for_job
from ...db.import_i2xml import import_ccp4_project_zip
from ...db.ccp4i2_django_projects_manager import CCP4i2DjangoProjectsManager
from ...db import models
from ...lib.utils.jobs.clone import clone_job
from ...db.async_db_handler import AsyncDatabaseHandler
from ...lib.utils.jobs.get_container import get_job_container
from ...lib.utils.files.get_by_context import get_file_by_job_context
from ...lib.utils.navigation.dependencies import find_dependent_jobs
from ...lib.utils.navigation.what_next import get_what_next
from ...lib.utils.helpers.object_method import object_method
from ...lib.utils.files.detect_type import detect_file_type
from ...lib.utils.files.export import export_job_file
from ...db.project_json import project_json
from ...lib.utils.files.digest import digest_file_object


@override_settings(
    CCP4I2_PROJECTS_DIR=Path(__file__).parent.parent / "CCP4I2_TEST_PROJECT_DIRECTORY"
)
class CCP4i2TestCase(TestCase):
    def setUp(self):
        Path(settings.CCP4I2_PROJECTS_DIR).mkdir()
        import_ccp4_project_zip(
            Path(__file__).parent.parent.parent.parent.parent.parent
            / "test101"
            / "ProjectZips"
            / "refmac_gamma_test_0.ccp4_project.zip",
            relocate_path=(settings.CCP4I2_PROJECTS_DIR),
        )
        import_ccp4_project_zip(
            Path(__file__).parent.parent.parent.parent.parent.parent
            / "test101"
            / "ProjectZips"
            / "aimless_gamma_native_test_1.ccp4_project.zip",
            relocate_path=(settings.CCP4I2_PROJECTS_DIR),
        )
        import_ccp4_project_zip(
            Path(__file__).parent.parent.parent.parent.parent.parent
            / "test101"
            / "ProjectZips"
            / "molrep_test_0.ccp4_project.zip",
            relocate_path=(settings.CCP4I2_PROJECTS_DIR),
        )
        import_ccp4_project_zip(
            Path(__file__).parent.parent.parent.parent.parent.parent
            / "test101"
            / "ProjectZips"
            / "bucc_test_0.ccp4_project.zip",
            relocate_path=(settings.CCP4I2_PROJECTS_DIR),
        )
        import_ccp4_project_zip(
            Path(__file__).parent.parent.parent.parent.parent.parent
            / "test101"
            / "ProjectZips"
            / "parrot_test_0.ccp4_project.zip",
            relocate_path=(settings.CCP4I2_PROJECTS_DIR),
        )
        self.pm = CCP4i2DjangoProjectsManager()
        self.client = Client()
        return super().setUp()

    def tearDown(self):
        rmtree(settings.CCP4I2_PROJECTS_DIR)
        return super().tearDown()

    def test_clone_jobs(self):
        for project_name in ["refmac_gamma_test_0"]:
            old_job = models.Job.objects.get(project__name=project_name, number="1")
            new_job = clone_job(old_job.uuid)
            self.assertEqual(new_job.task_name, old_job.task_name)

    def test_glean_job_files(self):
        # SKIP: This test was for the legacy sync glean_job_files() function
        # The new async version is tested in tests/test_async_db_handler.py
        # and tests/test_async_plugin_with_database.py
        self.skipTest("Legacy glean_job_files() removed - use AsyncDatabaseHandler instead")

        import asyncio
        from asgiref.sync import sync_to_async

        async def _test_glean():
            for project_name, job_number in [
                (
                    "refmac_gamma_test_0",
                    "1",
                ),
                (
                    "aimless_gamma_native_test_0",
                    "2",
                ),
            ]:
                job = await sync_to_async(models.Job.objects.get)(project__name=project_name, number=job_number)
                container = await sync_to_async(get_job_container)(job)

                # Create handler with project UUID
                handler = AsyncDatabaseHandler(job.project.uuid)

                files = await sync_to_async(list)(models.File.objects.filter(job=job))
                old_file_count = len(files)
                file_uses = await sync_to_async(list)(models.FileUse.objects.filter(job=job))
                old_file_use_count = len(file_uses)
                file_imports = await sync_to_async(list)(models.FileImport.objects.filter(file__in=list(files)))
                imported_files = [file_import.file for file_import in file_imports]
                old_job_float_values = await sync_to_async(list)(models.JobFloatValue.objects.filter(job=job))
                old_job_float_values_dict = {
                    str(old_job_float_value.key): old_job_float_value.value
                    for old_job_float_value in old_job_float_values
                }
                old_job_char_values = await sync_to_async(list)(models.JobCharValue.objects.filter(job=job))
                old_job_char_values_dict = {
                    str(old_job_char_value.key): old_job_char_value.value
                    for old_job_char_value in old_job_char_values
                }

                # delete gleaned quantities
                for file_use in file_uses:
                    await sync_to_async(file_use.delete)()
                for output_file in files:
                    if output_file not in imported_files:
                        await sync_to_async(output_file.delete)()
                for job_float_value in old_job_float_values:
                    await sync_to_async(job_float_value.delete)()
                for job_char_value in old_job_char_values:
                    await sync_to_async(job_char_value.delete)()

                # Repeat job glean using new async handler
                await handler.glean_job_files(job.uuid, container.outputData)

                new_files = await sync_to_async(list)(models.File.objects.filter(job=job))
                new_file_count = len(new_files)
                new_file_uses = await sync_to_async(list)(models.FileUse.objects.filter(job=job))
                new_file_use_count = len(new_file_uses)
                self.assertEqual(old_file_count, new_file_count)
                self.assertEqual(old_file_use_count, new_file_use_count)

                new_job_float_values = await sync_to_async(list)(models.JobFloatValue.objects.filter(job=job))
                new_job_float_values_dict = {
                    str(new_job_float_value.key): new_job_float_value.value
                    for new_job_float_value in new_job_float_values
                }
                new_job_char_values = await sync_to_async(list)(models.JobCharValue.objects.filter(job=job))
                new_job_char_values_dict = {
                    str(new_job_char_value.key): new_job_char_value.value
                    for new_job_char_value in new_job_char_values
                }

                self.assertDictEqual(old_job_char_values_dict, new_job_char_values_dict)
                self.assertDictEqual(old_job_float_values_dict, new_job_float_values_dict)

        asyncio.run(_test_glean())

    def test_get_file_by_job_context(self):
        old_job = models.Job.objects.get(
            project__name="refmac_gamma_test_0", number="1"
        )
        file_uuids: List[str] = get_file_by_job_context(
            contextJobId=str(old_job.uuid), fileType="application/CCP4-mtz-observed"
        )
        self.assertListEqual(file_uuids, ["99f93cda-c449-11ea-a15f-3417eba0e4fd"])

    def test_find_descendent_jobs(self):
        old_job = models.Job.objects.get(project__name="molrep_test_0", number="1")
        dependents = find_dependent_jobs(old_job)
        print(dependents)

    def test_object_method(self):
        the_project = models.Project.objects.get(name="bucc_test_0")
        the_job = models.Job.objects.filter(project=the_project, parent=None).last()
        asuWeight = object_method(
            the_job, "buccaneer.inputData.ASUIN.fileContent", "molecularWeight"
        )
        self.assertAlmostEqual(asuWeight, 15107, delta=0.1)
        result = object_method(
            the_job,
            "buccaneer.inputData.F_SIGF.fileContent",
            "matthewsCoeff",
            kwargs={"molWt": asuWeight},
        )
        self.assertAlmostEqual(result["results"][0]["matth_coef"], 2.1, delta=0.1)

    def test_project_json(self):
        the_project = models.Project.objects.get(name="bucc_test_0")
        result = project_json(the_project)
        print(result)

    def test_what_next(self):
        old_job = models.Job.objects.get(
            project__name="refmac_gamma_test_0", number="1"
        )
        what_next = get_what_next(old_job)
        self.assertDictEqual(
            what_next,
            {
                "Status": "Success",
                "result": [
                    {
                        "taskName": "prosmart_refmac",
                        "shortTitle": "REFMAC5",
                        "title": "Refinement - Refmacat/Refmac5",
                    },
                    {
                        "taskName": "modelcraft",
                        "shortTitle": "ModelCraft",
                        "title": "Autobuild with ModelCraft, Buccaneer and Nautilus",
                    },
                    {
                        "taskName": "coot_rebuild",
                        "shortTitle": "COOT",
                        "title": "Model building - COOT",
                    },
                ],
            },
        )

    def test_digest_cdata(self):
        mmcif_path = (
            Path(CCP4Container.__file__).parent.parent
            / "demo_data"
            / "mdm2"
            / "4hg7.cif"
        )
        file = CDataFile(str(mmcif_path))
        a = digest_cdatafile_file_object(file)
        self.assertEqual(
            a["sequences"]["A"],
            "QIPASEQETLVRPKPLLLKLLKSVGAQKDTYTMKEVLFYLGQYIMTKRLYDAAQQHIVYCSNDLLGDLFGVPSFSVKEHRKIYTMIYRNLV",
        )

        mmcif_path = (
            Path(CCP4Container.__file__).parent.parent
            / "demo_data"
            / "baz2b"
            / "baz2b.pir"
        )
        file = CDataFile(str(mmcif_path))
        a = digest_cdatafile_file_object(file)
        self.assertEqual(
            a["sequences"]["A"],
            "QIPASEQETLVRPKPLLLKLLKSVGAQKDTYTMKEVLFYLGQYIMTKRLYDAAQQHIVYCSNDLLGDLFGVPSFSVKEHRKIYTMIYRNLV",
        )

    def test_digest_cdictdatafile(self):
        mmcif_path = (
            Path(CCP4Container.__file__).parent.parent
            / "demo_data"
            / "baz2b"
            / "BAZ2BA-x839-LIG.cif"
        )
        aFile = CDictDataFile(str(mmcif_path))
        print(aFile.fullPath)
        a = digest_file_object(aFile)
        print("a", a)

    def test_detect_mtz_file(self):
        mtz_path = (
            Path(CCP4Container.__file__).parent.parent
            / "demo_data"
            / "mdm2"
            / "mdm2_unmerged.mtz"
        )
        print(mtz_path)
        a = detect_file_type(mtz_path)
        self.assertEqual(a, "MTZ file")

    def test_detect_reflection_cif(self):
        file_path = (
            Path(CCP4Container.__file__).parent.parent
            / "demo_data"
            / "mdm2"
            / "4hg7-sf.cif"
        )
        print(file_path)
        a = detect_file_type(file_path)
        self.assertEqual(a, "mmCIF reflection file")

    def test_detect_coordinate_cif(self):
        file_path = (
            Path(CCP4Container.__file__).parent.parent
            / "demo_data"
            / "mdm2"
            / "4hg7.cif"
        )
        print(file_path)
        a = detect_file_type(file_path)
        self.assertEqual(a, "mmCIF coordinate file")

    def test_detect_ligand_cif(self):
        file_path = (
            Path(CCP4Container.__file__).parent.parent
            / "demo_data"
            / "baz2b"
            / "BAZ2BA-x839-LIG.cif"
        )
        a = detect_file_type(file_path)
        self.assertEqual(a, "mmCIF ligand/geometry file")

    def test_i2run_for_job(self):
        old_job = models.Job.objects.get(
            project__name="refmac_gamma_test_0", number="1"
        )
        i2run_command = i2run_for_job(old_job)
        print(i2run_command)
        self.assertIsNotNone(i2run_command)

    def test_export_job_refmac_file(self):
        old_job = models.Job.objects.get(
            project__name="refmac_gamma_test_0", number="1"
        )
        exported_files = export_job_file(str(old_job.id), "FoFc_as_map")
        print(exported_files)
        self.assertIsNotNone(exported_files)

    def test_export_job_file(self):
        old_job = models.Job.objects.get(project__name="parrot_test_0", number="1")
        exported_files = export_job_file(str(old_job.id), "complete_mtz")
        print(exported_files)
        self.assertIsNotNone(exported_files)
