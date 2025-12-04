"""Test database-aware path handling for CDataFile objects."""

import shutil
import pathlib
import os
from xml.etree import ElementTree as ET

def test_plugin_script_initialization():
    """Test that plugin scripts can be initialized with database-aware path handling."""
    # Import Django models inside the test function after Django is configured
    from django.core.management import call_command
    from ccp4x.db.models import Project, Job
    from core.CCP4TaskManager import TASKMANAGER
    from core.CCP4PluginScript import CPluginScript
    from ccp4x.db.async_db_handler import AsyncDatabaseHandler
    print("[TEST] Testing plugin script initialization for 'parrot'")
    call_command("create_project", "test_project", "--directory", "./CCP4_PROJECTS/test_project")
    call_command("create_job", "-pn", "test_project", "-tn", "prosmart_refmac")
    prosmart_refmac_job = Job.objects.get(project__name="test_project", task_name="prosmart_refmac")
    prosmart_refmac_plugin:CPluginScript = TASKMANAGER().get_plugin_class("prosmart_refmac")(work_directory=str(prosmart_refmac_job.directory))
    src = pathlib.Path(os.environ["CCP4I2_ROOT"]) / "demo_data" / "gamma" /"gamma_model.pdb"
    dest = pathlib.Path(prosmart_refmac_job.directory) / "input_coords.pdb"
    shutil.copyfile(src, dest)
    db_handler = AsyncDatabaseHandler(project_uuid=prosmart_refmac_job.project.uuid)
    prosmart_refmac_plugin.setDbData(
        handler=db_handler,
        projectName=prosmart_refmac_job.project.name,
        projectId=str(prosmart_refmac_job.project.uuid),
        jobNumber=prosmart_refmac_job.number,
        jobId=str(prosmart_refmac_job.uuid)
    )
    prosmart_refmac_plugin.container.inputData.XYZIN.setFullPath(str(dest))
    
    XYZIN_etree = prosmart_refmac_plugin.container.inputData.XYZIN.getEtree()
    ET.indent(XYZIN_etree, space="  ")
    print(ET.tostring(XYZIN_etree, encoding="unicode"))
    assert prosmart_refmac_plugin is  None