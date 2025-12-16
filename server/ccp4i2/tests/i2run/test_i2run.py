import shlex
import os
from pathlib import Path
from shutil import rmtree
import uuid
from argparse import ArgumentParser
from django.test import TestCase, override_settings
from django.conf import settings
from ccp4i2.db.models import Job, Project, File, JobCharValue, JobFloatValue
from ccp4i2.db.import_i2xml import import_ccp4_project_zip

from ccp4i2.cli.i2run.CCP4i2RunnerDjango import CCP4i2RunnerDjango
from ccp4i2.core import CCP4Modules

# Set CCP4I2_TOP to the project root
os.environ.setdefault(
    "CCP4I2_TOP", str(Path(__file__).parent.parent.parent.parent.parent.parent)
)

case1 = """aimless_pipe \
    --UNMERGEDFILES \
        crystalName=hg7 \
        dataset=DS1 \
        file=$CCP4I2_TOP/demo_data/mdm2/mdm2_unmerged.mtz \
    --project_name refmac_gamma_test_0"""


case2a = """aimless_pipe \
    --UNMERGEDFILES \
        crystalName=hg7 \
        dataset=DS1 \
        file=$CCP4I2_TOP/demo_data/mdm2/mdm2_unmerged.mtz \
    --XYZIN_REF fullPath=$CCP4I2_TOP/demo_data/mdm2/4hg7.pdb \
	--MODE MATCH \
	--REFERENCE_DATASET XYZ \
    --project_name refmac_gamma_test_0"""

case2b = """prosmart_refmac \
    --F_SIGF fileUse="SubstituteLigand[-1].F_SIGF_OUT" \
    --XYZIN \
        fullPath=$CCP4I2_TOP/demo_data/mdm2/4hg7.pdb \
        selection/text="not (HOH)" \
    --prosmartProtein.REFERENCE_MODELS \
        fullPath=$CCP4I2_TOP/demo_data/mdm2/4qo4.cif \
    --project_name SubstituteLigand_test_0"""

case3 = """phaser_simple \
    --F_SIGF \
        fullPath=$CCP4I2_TOP/demo_data/beta_blip/beta_blip_P3221.mtz \
        columnLabels="/*/*/[Fobs,Sigma]" \
    --F_OR_I F \
    --XYZIN \
        $CCP4I2_TOP/demo_data/beta_blip/beta.pdb \
    --project_name refmac_gamma_test_0"""

OLD_CCP4I2_PROJECTS_DIR = settings.CCP4I2_PROJECTS_DIR


@override_settings(
    CCP4I2_PROJECTS_DIR=Path(__file__).parent.parent.parent.parent.parent
    / "CCP4I2_TEST_PROJECT_DIRECTORY"
)
class CCP4i2TestCase(TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.app = CCP4Modules.QTAPPLICATION(graphical=False)

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
            / "SubstituteLigand_test_0.ccp4_project.zip",
            relocate_path=(settings.CCP4I2_PROJECTS_DIR),
        )
        return super().setUp()

    def tearDown(self):
        rmtree(settings.CCP4I2_PROJECTS_DIR)
        return super().tearDown()

    def test_shlex(self):
        self.assertEqual(len(shlex.split(case1, comments=True)), 7)

    def test_case1(self):
        args = shlex.split(os.path.expandvars(case1), comments=True)
        print(args)
        i2Runner = CCP4i2RunnerDjango(
            the_args=args,
            parser=ArgumentParser(),
            parent=self.app,
        )
        # i2Runner.parseArgs()
        print("Initial file count", File.objects.count())
        i2Runner.execute()
        # self.app.exit()
        self.assertEqual(Job.objects.last().project.name, "refmac_gamma_test_0")
        self.assertEqual(Job.objects.last().number, "2.5")
        the_job = Job.objects.get(uuid=uuid.UUID(i2Runner.jobId))
        self.assertEqual(JobCharValue.objects.filter(job=the_job)[0].value, "P 61 2 2")

    def test_case2(self):
        args = shlex.split(os.path.expandvars(case2b), comments=True)
        print(args)
        self.assertEqual(
            Project.objects.filter(name="SubstituteLigand_test_0").count(), 1
        )
        i2Runner = CCP4i2RunnerDjango(
            the_args=args,
            parser=ArgumentParser(),
            parent=self.app,
        )
        i2Runner.execute()
        # self.app.exit()
        self.assertEqual(Job.objects.last().project.name, "SubstituteLigand_test_0")
        self.assertEqual(Job.objects.filter(parent__isnull=True).last().number, "2")
        the_job = Job.objects.get(uuid=uuid.UUID(i2Runner.jobId))
        # print(glob.glob(str(the_job.directory / "*")))
        self.assertEqual(the_job.status, Job.Status.FINISHED)
        self.assertAlmostEqual(
            JobFloatValue.objects.filter(job=the_job)[0].value, 0.253
        )

    def test_case3(self):
        args = shlex.split(os.path.expandvars(case3), comments=True)
        i2Runner = CCP4i2RunnerDjango(
            the_args=args,
            parser=ArgumentParser(),
            parent=self.app,
        )
        try:
            i2Runner.execute()
        except Exception as e:
            with open(Path(i2Runner.work_directory) / "input_params.xml") as f:
                log = f.read()
                print(log)
