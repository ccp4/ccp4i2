"""Debug test for JSON encoder to track down the inputData.XYZIN issue."""
import json
import os
from pathlib import Path
from shutil import rmtree
from django.test import TestCase, override_settings
from django.conf import settings

from ccp4x.db.import_i2xml import import_ccp4_project_zip
from ccp4x.db.models import Job
from ccp4x.lib.utils.plugins.get_plugin import get_job_plugin
from ccp4x.lib.utils.containers.json_encoder import CCP4i2JsonEncoder


@override_settings(
    CCP4I2_PROJECTS_DIR=Path(__file__).parent.parent.parent.parent.parent
    / "CCP4I2_TEST_PROJECT_DIRECTORY"
)
class JsonEncoderDebugTest(TestCase):
    def setUp(self):
        Path(settings.CCP4I2_PROJECTS_DIR).mkdir(exist_ok=True)
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

    def test_xyzin_serialization(self):
        """Test that inputData.XYZIN children are properly serialized."""
        job = Job.objects.get(project__name="refmac_gamma_test_0", number="1")
        plugin = get_job_plugin(job)
        container = plugin.container

        # Access XYZIN directly
        input_data = container.inputData
        xyzin = input_data.XYZIN

        print("\n=== XYZIN Debug ===")
        print(f"XYZIN type: {type(xyzin)}")
        print(f"XYZIN class name: {xyzin.__class__.__name__}")

        # Check children
        children = xyzin.children()
        print(f"\nXYZIN has {len(children)} children:")
        for child in children:
            child_name = child.objectName() if hasattr(child, 'objectName') else str(child)
            child_value = getattr(child, '_value', 'NO_VALUE_ATTR')
            print(f"  - {child_name}: {type(child).__name__} = {repr(child_value)}")

        # Check baseName specifically
        if hasattr(xyzin, 'baseName'):
            base_name = xyzin.baseName
            print(f"\nxyzin.baseName direct access:")
            print(f"  type: {type(base_name)}")
            print(f"  _value: {repr(getattr(base_name, '_value', 'NO ATTR'))}")
            print(f"  str(): {str(base_name)}")

            # Check if baseName is in children
            base_name_in_children = any(
                c.objectName() == 'baseName' for c in children if hasattr(c, 'objectName')
            )
            print(f"  baseName in children(): {base_name_in_children}")

        # Check project
        if hasattr(xyzin, 'project'):
            project = xyzin.project
            print(f"\nxyzin.project direct access:")
            print(f"  type: {type(project)}")
            print(f"  _value: {repr(getattr(project, '_value', 'NO ATTR'))}")

        # Check relPath
        if hasattr(xyzin, 'relPath'):
            rel_path = xyzin.relPath
            print(f"\nxyzin.relPath direct access:")
            print(f"  type: {type(rel_path)}")
            print(f"  _value: {repr(getattr(rel_path, '_value', 'NO ATTR'))}")

        # Now serialize and check
        print("\n=== JSON Serialization ===")
        json_str = json.dumps(xyzin, cls=CCP4i2JsonEncoder, indent=2)
        json_obj = json.loads(json_str)

        print(f"Serialized XYZIN keys: {list(json_obj.get('_value', {}).keys())}")

        # Check specific values
        xyzin_value = json_obj.get('_value', {})
        if 'baseName' in xyzin_value:
            print(f"baseName._value in JSON: {repr(xyzin_value['baseName'].get('_value'))}")
        else:
            print("baseName NOT in _value dict!")

        if 'project' in xyzin_value:
            print(f"project._value in JSON: {repr(xyzin_value['project'].get('_value'))}")
        else:
            print("project NOT in _value dict!")

        # Assert that values are present
        self.assertIn('baseName', xyzin_value, "baseName should be in serialized _value")
        self.assertIsNotNone(xyzin_value.get('baseName', {}).get('_value'),
                            "baseName._value should not be None/empty")
