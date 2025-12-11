"""
Test DEF XML workflow: loading, modifying, and verifying CData structures.

Note: This test focuses on DEF XML parsing and parameter modification.
The params XML export/import functionality will be tested separately once implemented.
"""

import sys
import os
import tempfile
from pathlib import Path

import pytest

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from ccp4i2.core.task_manager.def_xml_handler import parse_def_xml_file
from ccp4i2.core.base_object.base_classes import ValueState

# Sample DEF XML for testing workflow
SAMPLE_DEF_XML = """<?xml version="1.0" encoding="UTF-8"?>
<ns0:ccp4i2 xmlns:ns0="http://www.ccp4.ac.uk/ccp4ns">
  <ccp4i2_header>
    <pluginName>servalcat_pipe</pluginName>
  </ccp4i2_header>
  <ccp4i2_body id="servalcat_pipe">
    <ccp4i2_body>
      <container id="inputData">
        <content id="XYZIN">
          <className>CPdbDataFile</className>
          <qualifiers>
            <mustExist>True</mustExist>
            <toolTip>Input coordinate file</toolTip>
          </qualifiers>
        </content>
      </container>
      <container id="outputData">
        <content id="XYZOUT">
          <className>CPdbDataFile</className>
          <qualifiers>
            <toolTip>Output coordinate file</toolTip>
          </qualifiers>
        </content>
      </container>
      <container id="controlParameters">
        <content id="DATA_METHOD">
          <className>CString</className>
          <qualifiers>
            <onlyEnumerators>True</onlyEnumerators>
            <enumerators>xtal,spa</enumerators>
            <default>xtal</default>
          </qualifiers>
        </content>
        <content id="ADD_WATERS">
          <className>CBoolean</className>
          <qualifiers>
            <default>False</default>
            <toolTip>Add water molecules</toolTip>
          </qualifiers>
        </content>
        <content id="NCYCLES">
          <className>CInt</className>
          <qualifiers>
            <default>10</default>
            <min>1</min>
            <max>100</max>
            <toolTip>Number of refinement cycles</toolTip>
          </qualifiers>
        </content>
        <content id="WEIGHT">
          <className>CFloat</className>
          <qualifiers>
            <min>0.0</min>
            <toolTip>Refinement weight</toolTip>
          </qualifiers>
        </content>
        <content id="B_REFINEMENT_MODE">
          <className>CString</className>
          <qualifiers>
            <onlyEnumerators>True</onlyEnumerators>
            <enumerators>iso,aniso,fix</enumerators>
            <default>iso</default>
            <toolTip>B-factor refinement mode</toolTip>
          </qualifiers>
        </content>
        <content id="OCCUPANCY_REFINEMENT">
          <className>CBoolean</className>
          <qualifiers>
            <default>True</default>
            <toolTip>Refine occupancies</toolTip>
          </qualifiers>
        </content>
      </container>
    </ccp4i2_body>
  </ccp4i2_body>
  <container id="metalCoordPipeline">
    <content id="RUN_METALCOORD">
      <className>CBoolean</className>
      <qualifiers>
        <default>False</default>
      </qualifiers>
    </content>
    <content id="LINKS">
      <className>CString</className>
      <qualifiers>
        <onlyEnumerators>True</onlyEnumerators>
        <enumerators>UPDATE,KEEP,NOTTOUCH</enumerators>
        <default>UPDATE</default>
      </qualifiers>
    </content>
  </container>
</ns0:ccp4i2>"""


@pytest.fixture
def def_xml_file():
    """Create a temporary DEF XML file."""
    with tempfile.NamedTemporaryFile(
        mode='w', suffix='.def.xml', delete=False
    ) as f:
        f.write(SAMPLE_DEF_XML)
        temp_path = f.name

    yield temp_path

    # Cleanup
    Path(temp_path).unlink(missing_ok=True)


def test_load_task_definition(def_xml_file):
    """Test loading task definition from DEF XML."""
    task = parse_def_xml_file(def_xml_file)

    assert task is not None
    assert task.name == 'servalcat_pipe'
    assert type(task).__name__ == 'CContainer'


def test_initial_default_states(def_xml_file):
    """Test that initial parameters have default values and DEFAULT state."""
    task = parse_def_xml_file(def_xml_file)

    ctrl = task.controlParameters

    # Check default values
    assert ctrl.DATA_METHOD.value == 'xtal'
    assert ctrl.ADD_WATERS.value is False
    assert ctrl.NCYCLES.value == 10
    assert ctrl.B_REFINEMENT_MODE.value == 'iso'
    assert ctrl.OCCUPANCY_REFINEMENT.value is True

    # Check value states (defaults should be DEFAULT state)
    # Note: This depends on implementation of getValueState
    if hasattr(ctrl.NCYCLES, 'get_value_state'):
        state = ctrl.NCYCLES.get_value_state('value')
        assert state in [ValueState.DEFAULT, ValueState.NOT_SET]


def test_parameter_modification_direct_assignment(def_xml_file):
    """Test modifying parameters using direct assignment."""
    task = parse_def_xml_file(def_xml_file)

    ctrl = task.controlParameters

    # Modify using direct assignment
    ctrl.ADD_WATERS.value = True
    ctrl.NCYCLES.value = 25
    ctrl.B_REFINEMENT_MODE.value = 'aniso'

    # Verify modifications
    assert ctrl.ADD_WATERS.value is True
    assert ctrl.NCYCLES.value == 25
    assert ctrl.B_REFINEMENT_MODE.value == 'aniso'


def test_parameter_modification_multiple_types(def_xml_file):
    """Test modifying parameters of different types."""
    task = parse_def_xml_file(def_xml_file)

    ctrl = task.controlParameters

    # Test CInt modification
    original_cycles = ctrl.NCYCLES.value
    assert original_cycles == 10
    ctrl.NCYCLES.value = 25
    assert ctrl.NCYCLES.value == 25

    # Test CFloat modification
    ctrl.WEIGHT.value = 0.15
    assert ctrl.WEIGHT.value == 0.15

    # Test CBoolean modification
    ctrl.ADD_WATERS.value = True
    assert ctrl.ADD_WATERS.value is True

    # Test CString modification
    ctrl.B_REFINEMENT_MODE.value = 'aniso'
    assert ctrl.B_REFINEMENT_MODE.value == 'aniso'


def test_nested_container_modification(def_xml_file):
    """Test modifying parameters in nested containers."""
    task = parse_def_xml_file(def_xml_file)

    # Modify in metalCoordPipeline container
    task.metalCoordPipeline.RUN_METALCOORD.value = True
    task.metalCoordPipeline.LINKS.value = 'KEEP'

    # Verify modifications
    assert task.metalCoordPipeline.RUN_METALCOORD.value is True
    assert task.metalCoordPipeline.LINKS.value == 'KEEP'


def test_constraint_validation(def_xml_file):
    """Test that constraints (min/max) are validated."""
    task = parse_def_xml_file(def_xml_file)

    ncycles = task.controlParameters.NCYCLES

    # Valid values should work
    ncycles.value = 50
    assert ncycles.value == 50

    ncycles.value = 1  # min value
    assert ncycles.value == 1

    ncycles.value = 100  # max value
    assert ncycles.value == 100

    # Invalid values should raise ValueError
    with pytest.raises(ValueError):
        ncycles.value = 0  # Below min

    with pytest.raises(ValueError):
        ncycles.value = 101  # Above max

    with pytest.raises(ValueError):
        ncycles.value = -5  # Way below min


def test_float_constraint_validation(def_xml_file):
    """Test that float constraints are validated."""
    task = parse_def_xml_file(def_xml_file)

    weight = task.controlParameters.WEIGHT

    # Valid values should work
    weight.value = 0.15
    assert weight.value == 0.15

    weight.value = 0.0  # min value
    assert weight.value == 0.0

    # Below min should raise ValueError
    with pytest.raises(ValueError):
        weight.value = -0.1


def test_fresh_task_loading(def_xml_file):
    """Test that loading a fresh task gives clean defaults."""
    # Load first task and modify it
    task1 = parse_def_xml_file(def_xml_file)
    task1.controlParameters.NCYCLES.value = 25
    task1.controlParameters.ADD_WATERS.value = True

    # Load second task - should have defaults
    task2 = parse_def_xml_file(def_xml_file)

    assert task2.controlParameters.NCYCLES.value == 10  # default
    assert task2.controlParameters.ADD_WATERS.value is False  # default


def test_multiple_containers_structure(def_xml_file):
    """Test that multiple containers are properly structured."""
    task = parse_def_xml_file(def_xml_file)

    # Check all expected containers exist
    assert hasattr(task, 'inputData')
    assert hasattr(task, 'outputData')
    assert hasattr(task, 'controlParameters')
    assert hasattr(task, 'metalCoordPipeline')

    # Check container types
    assert type(task.inputData).__name__ == 'CContainer'
    assert type(task.outputData).__name__ == 'CContainer'
    assert type(task.controlParameters).__name__ == 'CContainer'
    assert type(task.metalCoordPipeline).__name__ == 'CContainer'


def test_file_objects_present(def_xml_file):
    """Test that file objects are present and correctly typed."""
    task = parse_def_xml_file(def_xml_file)

    # Check input files
    assert hasattr(task.inputData, 'XYZIN')
    assert type(task.inputData.XYZIN).__name__ == 'CPdbDataFile'

    # Check output files
    assert hasattr(task.outputData, 'XYZOUT')
    assert type(task.outputData.XYZOUT).__name__ == 'CPdbDataFile'


def test_enumerator_constraints(def_xml_file):
    """Test parameters with enumerator constraints."""
    task = parse_def_xml_file(def_xml_file)

    data_method = task.controlParameters.DATA_METHOD

    # Test valid enumerator values
    assert data_method.value == 'xtal'  # default

    data_method.value = 'spa'
    assert data_method.value == 'spa'

    data_method.value = 'xtal'
    assert data_method.value == 'xtal'

    # Test B_REFINEMENT_MODE enumerators
    b_mode = task.controlParameters.B_REFINEMENT_MODE
    assert b_mode.value == 'iso'  # default

    b_mode.value = 'aniso'
    assert b_mode.value == 'aniso'

    b_mode.value = 'fix'
    assert b_mode.value == 'fix'


def test_boolean_values(def_xml_file):
    """Test boolean parameter handling."""
    task = parse_def_xml_file(def_xml_file)

    # Test ADD_WATERS
    add_waters = task.controlParameters.ADD_WATERS
    assert add_waters.value is False  # default

    add_waters.value = True
    assert add_waters.value is True

    add_waters.value = False
    assert add_waters.value is False

    # Test OCCUPANCY_REFINEMENT
    occ_ref = task.controlParameters.OCCUPANCY_REFINEMENT
    assert occ_ref.value is True  # default

    occ_ref.value = False
    assert occ_ref.value is False


def test_complex_workflow_scenario(def_xml_file):
    """Test a complex workflow with multiple modifications."""
    # Step 1: Load task
    task = parse_def_xml_file(def_xml_file)

    # Step 2: Verify initial state
    assert task.controlParameters.NCYCLES.value == 10
    assert task.controlParameters.ADD_WATERS.value is False
    assert task.metalCoordPipeline.RUN_METALCOORD.value is False

    # Step 3: Make multiple modifications
    task.controlParameters.NCYCLES.value = 25
    task.controlParameters.ADD_WATERS.value = True
    task.controlParameters.WEIGHT.value = 0.15
    task.controlParameters.B_REFINEMENT_MODE.value = 'aniso'
    task.controlParameters.OCCUPANCY_REFINEMENT.value = False
    task.metalCoordPipeline.RUN_METALCOORD.value = True
    task.metalCoordPipeline.LINKS.value = 'KEEP'

    # Step 4: Verify all modifications
    assert task.controlParameters.NCYCLES.value == 25
    assert task.controlParameters.ADD_WATERS.value is True
    assert task.controlParameters.WEIGHT.value == 0.15
    assert task.controlParameters.B_REFINEMENT_MODE.value == 'aniso'
    assert task.controlParameters.OCCUPANCY_REFINEMENT.value is False
    assert task.metalCoordPipeline.RUN_METALCOORD.value is True
    assert task.metalCoordPipeline.LINKS.value == 'KEEP'

    # Step 5: Load fresh task and verify it has defaults
    fresh_task = parse_def_xml_file(def_xml_file)
    assert fresh_task.controlParameters.NCYCLES.value == 10
    assert fresh_task.controlParameters.ADD_WATERS.value is False
    assert fresh_task.metalCoordPipeline.RUN_METALCOORD.value is False
