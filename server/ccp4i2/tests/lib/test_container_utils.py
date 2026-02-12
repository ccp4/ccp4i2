import os
from pathlib import Path
from shutil import rmtree
import subprocess
from xml.etree import ElementTree as ET
from django.test import TestCase, override_settings
from django.conf import settings
import gemmi
from ccp4i2.core import CCP4PerformanceData
from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core import CCP4Data
from ccp4i2.core import CCP4Container
from ccp4i2.core import CCP4TaskManager
from ...db.models import Job, File
from ...db.import_i2xml import import_ccp4_project_zip

from ...lib.utils.plugins.get_plugin import get_job_plugin
from ...lib.utils.formats.mtz import mtz_as_dict
from ...lib.utils.parameters.unset_output_data import unset_output_data
from ...lib.utils.containers.remove_defaults import (
    remove_container_default_values,
)
from ...lib.utils.containers.find_objects import find_objects
from ...lib.utils.parameters.load_xml import load_nested_xml
from ...lib.utils.containers.validate import validate_container
from ...lib.utils.jobs.clone import clone_job
from ...lib.utils.jobs.create import create_job
from ...lib.utils.reporting.i2_report import get_report_job_info
from ...lib.utils.formats.gemmi_split_mtz import gemmi_split_mtz


@override_settings(
    CCP4I2_PROJECTS_DIR=Path(__file__).parent.parent.parent.parent.parent
    / "CCP4I2_TEST_PROJECT_DIRECTORY"
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
        return super().setUp()

    def tearDown(self):
        rmtree(settings.CCP4I2_PROJECTS_DIR)
        return super().tearDown()

    def test_unset_output_data(self):
        job = Job.objects.get(pk=1)
        the_job_plugin = get_job_plugin(job)
        unset_output_data(the_job_plugin)
        et = the_job_plugin.container.getEtree()
        ET.indent(et, "\t", 0)
        self.assertEqual(1, 1)

    def test_remove_container_default_values(self):
        job = Job.objects.get(pk=1)
        the_job_plugin = get_job_plugin(job)
        remove_container_default_values(the_job_plugin.container)
        container: CCP4Container.CContainer = the_job_plugin.container
        et = container.getEtree()
        ET.indent(et, "\t", 0)
        self.assertEqual(1, 1)

    def test_find_objects(self):
        job = Job.objects.get(pk=1)
        the_job_plugin = get_job_plugin(job)
        container: CCP4Container.CContainer = the_job_plugin.container
        kpis = find_objects(
            container,
            lambda a: isinstance(a, CCP4PerformanceData.CPerformanceIndicator),
            True,
        )
        for kpi in kpis:
            for kpi_param_name in kpi.dataOrder():
                value = getattr(kpi, kpi_param_name)
                self.assertTrue(isinstance(value, CCP4Data.CString))
                self.assertTreu(isinstance(value, CCP4Data.CFloat))

    def test_load_nested_xml(self):
        startXML = ET.fromstring(prosmart_defmac_xml)
        result = load_nested_xml(startXML)
        ET.indent(result, " ")
        print(ET.tostring(result).decode("utf-8"))

    def test_validate_container(self):
        job = Job.objects.get(pk=1)
        the_job_plugin = get_job_plugin(job)
        container: CCP4Container.CContainer = the_job_plugin.container
        error_etree: ET.Element = validate_container(container)

        self.assertEqual(len(error_etree.findall(".//errorReport")), 47)
        ET.indent(error_etree, " ")
        # print(ET.tostring(err))

    def test_validate_container_no_XYZIN(self):
        job = Job.objects.get(project__name="refmac_gamma_test_0", number="1")
        cloned_job = clone_job(job.uuid)
        the_job_plugin = get_job_plugin(cloned_job)
        container: CCP4Container.CContainer = the_job_plugin.container
        XYZIN = container.find("XYZIN")
        XYZIN.unSet()
        error_etree: ET.Element = validate_container(container)
        ET.indent(error_etree, " ")
        print(ET.tostring(error_etree).decode("utf-8"))
        self.assertEqual(len(error_etree.findall(".//errorReport")), 45)

    def test_validate_container_NCYCLES_MINUS_1(self):
        job = Job.objects.get(project__name="refmac_gamma_test_0", number="1")
        cloned_job = clone_job(job.uuid)
        the_job_plugin = get_job_plugin(cloned_job)
        container: CCP4Container.CContainer = the_job_plugin.container
        NCYCLES = container.find("NCYCLES")
        with self.assertRaises(CCP4ErrorHandling.CException):
            NCYCLES.set(-1)

    def test_def_xml_container(self):
        taskManager: CCP4TaskManager.CTaskManager = CCP4TaskManager.CTaskManager()
        defFile = taskManager.locate_def_xml(task_name="prosmart_refmac")
        container: CCP4Container.CContainer = CCP4Container.CContainer(
            definitionFile=defFile
        )
        def_etree = container.getEtree()
        # ET.indent(def_etree, " ")
        print(ET.tostring(def_etree).decode("utf-8"))

    def test_get_task_tree(self):
        result = CCP4TaskManager.get_task_tree()
        self.assertEqual(len(result["lookup"].items()), 135)
        self.assertEqual(len(result["tree"]), 17)

    def test_get_report_job_info(self):
        job = Job.objects.get(project__name="refmac_gamma_test_0", number="1")
        result = get_report_job_info(job.uuid)
        self.assertEqual(
            result["fileroot"],
            str(settings.CCP4I2_PROJECTS_DIR / "refmac_gamma_test_0/CCP4_JOBS/job_1/"),
        )

    def test_mtz_as_dict(self):
        mtz_file = File.objects.filter(type="application/CCP4-mtz-observed").first()
        result = mtz_as_dict(mtz_file.path, with_reflections=True)
        print(result)
        self.assertEqual(result["header_info"]["title"], "From Clipper CCP4MTZfile")
        self.assertEqual(result["header_info"]["spacegroup"], "P 21 21 21")

    def test_gemmi_split_mtz(self):
        beta_mtz = os.path.expandvars(
            "$CCP4I2_TOP/demo_data/beta_blip/beta_blip_P3221.mtz"
        )
        unsplit_mtz = gemmi.read_mtz_file(beta_mtz)
        uniques = gemmi.make_miller_array(
            unsplit_mtz.cell,
            unsplit_mtz.spacegroup,
            unsplit_mtz.resolution_high(),
            unsplit_mtz.resolution_low(),
            unique=True,
        )
        print("size", uniques.size / 3)
        print(unsplit_mtz.nreflections)
        split_mtz_path = gemmi_split_mtz(
            Path(beta_mtz), "/*/*/[Fobs,Sigma]", Path("/tmp/beta_blip.mtz")
        )
        split_mtz = gemmi.read_mtz_file(str(split_mtz_path))
        subprocess.call(["mtzdmp", str(beta_mtz), "10"])
        subprocess.call(["mtzdmp", str(split_mtz_path), "10"])

    def test_create_coot_patch_file_paths(self):
        job = Job.objects.get(project__name="refmac_gamma_test_0", number="1")
        coot_job_uuid = create_job(
            projectId=job.project.uuid, taskName="coot_rebuild", saveParams=True
        )
        the_job = Job.objects.get(uuid=coot_job_uuid)
        the_job_plugin = get_job_plugin(the_job)

        # Verify the output path ends with the expected filename pattern
        output_path = str(the_job_plugin.container.outputData.XYZOUT[-1])
        self.assertTrue(
            output_path.endswith("coot_rebuild_9.pdb"),
            f"Output path {output_path} should end with coot_rebuild_9.pdb",
        )


prosmart_defmac_xml = """<ns0:ccp4i2 xmlns:ns0="http://www.ccp4.ac.uk/ccp4ns">
  <ccp4i2_header>
    <function>DEF</function>
    <comment/>
    <creationTime>14:00 19/Jul/12</creationTime>
    <userId>cowtan</userId>
    <jobId/>
    <project/>
    <pluginName>prosmart_refmac</pluginName>
    <pluginVersion/>
    <jobNumber/>
  </ccp4i2_header>
  <ccp4i2_body id="prosmart_refmac">
    <file>
      <CI2XmlDataFile>
        <project>CCP4I2_TOP</project>
        <relPath>wrappers/refmac/script</relPath>
        <baseName>refmac.def.xml</baseName>
      </CI2XmlDataFile>
    </file>
    <container id="inputData">
      <content id="AMINOACID_CHAINS">
        <className>CList</className>
        <qualifiers>
          <allowUndefined>True</allowUndefined>
        </qualifiers>
      </content>
      <content id="NUCLEOTIDE_CHAINS">
        <className>CList</className>
        <qualifiers>
          <allowUndefined>True</allowUndefined>
        </qualifiers>
      </content>
    </container>
  </ccp4i2_body>
</ns0:ccp4i2>"""
