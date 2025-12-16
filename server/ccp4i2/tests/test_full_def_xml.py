"""Test DEF XML parser with full servalcat_pipe task definition."""

import tempfile
from pathlib import Path

import pytest

from ccp4i2.core.task_manager.def_xml_handler import parse_def_xml_file

# Full XML from the servalcat_pipe task
FULL_SERVALCAT_XML = """<?xml version="1.0" encoding="UTF-8"?>
<ns0:ccp4i2 xmlns:ns0="http://www.ccp4.ac.uk/ccp4ns">
  <ccp4i2_header>
    <function>DEF</function>
    <comment/>
    <creationTime>14:00 19/Jul/12</creationTime>
    <userId>martinmaly</userId>
    <ccp4iVersion>0.0.1</ccp4iVersion>
    <jobId/>
    <project/>
    <pluginName>servalcat_pipe</pluginName>
    <pluginVersion/>
    <jobNumber/>
  </ccp4i2_header>
  <ccp4i2_body id="servalcat_pipe">
    <ccp4i2_body>
      <container id="inputData">
        <content id="XYZIN">
          <className>CPdbDataFile</className>
          <qualifiers>
            <ifAtomSelection>True</ifAtomSelection>
            <mustExist>True</mustExist>
            <allowUndefined>False</allowUndefined>
            <fromPreviousJob>True</fromPreviousJob>
            <requiredSubType>1,0</requiredSubType>
            <toolTip>File containing model coordinates (PDB/mmCIF).</toolTip>
          </qualifiers>
        </content>
        <content id="HKLIN">
          <className>CObsDataFile</className>
          <qualifiers>
            <mustExist>True</mustExist>
            <allowUndefined>True</allowUndefined>
            <fromPreviousJob>True</fromPreviousJob>
            <toolTip>File containing structure factor amplitudes/intensities.</toolTip>
          </qualifiers>
        </content>
        <content id="DICT_LIST">
          <className>CList</className>
          <qualifiers>
            <listMinLength>0</listMinLength>
          </qualifiers>
          <subItem>
            <className>CDictDataFile</className>
            <qualifiers>
              <default>
                <contentFlag>1</contentFlag>
              </default>
              <toolTip>Restraint dictionary (mmCIF file).</toolTip>
              <mustExist>True</mustExist>
            </qualifiers>
          </subItem>
        </content>
      </container>
      <container id="outputData">
        <content id="XYZOUT">
          <className>CPdbDataFile</className>
          <qualifiers>
            <default>
              <subType>1</subType>
              <contentFlag>1</contentFlag>
            </default>
            <saveToDb>True</saveToDb>
          </qualifiers>
        </content>
        <content id="FPHIOUT">
          <className>CMapCoeffsDataFile</className>
          <qualifiers>
            <default>
              <subType>1</subType>
              <contentFlag>1</contentFlag>
            </default>
            <saveToDb>True</saveToDb>
          </qualifiers>
        </content>
      </container>
      <container id="controlParameters">
        <content id="DATA_METHOD">
          <className>CString</className>
          <qualifiers>
            <onlyEnumerators>True</onlyEnumerators>
            <menuText>Diffraction data,SPA maps</menuText>
            <enumerators>xtal,spa</enumerators>
            <default>xtal</default>
            <allowUndefined>False</allowUndefined>
          </qualifiers>
        </content>
        <content id="ADD_WATERS">
          <className>CBoolean</className>
          <qualifiers>
            <default>False</default>
            <toolTip>Add waters and perform further refinement.</toolTip>
          </qualifiers>
        </content>
        <content id="NCYCLES">
          <className>CInt</className>
          <qualifiers>
            <default>10</default>
            <min>0</min>
            <toolTip>Number of refinement cycles to perform.</toolTip>
          </qualifiers>
        </content>
        <content id="WEIGHT">
          <className>CFloat</className>
          <qualifiers>
            <min>0.0</min>
            <toolTip>Weight controlling data and restraint terms.</toolTip>
          </qualifiers>
        </content>
        <content id="B_REFINEMENT_MODE">
          <className>CString</className>
          <qualifiers>
            <onlyEnumerators>True</onlyEnumerators>
            <menuText>isotropic,anisotropic,fixed</menuText>
            <enumerators>iso,aniso,fix</enumerators>
            <default>iso</default>
            <toolTip>Specifies ADP parameterisation.</toolTip>
          </qualifiers>
        </content>
        <content id="OCCUPANCY_REFINEMENT">
          <className>CBoolean</className>
          <qualifiers>
            <default>True</default>
            <toolTip>Refine subunitary atomic occupancy parameters.</toolTip>
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
        <menuText>update,keep,keep as they are</menuText>
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
        f.write(FULL_SERVALCAT_XML)
        temp_path = f.name

    yield temp_path

    # Cleanup
    Path(temp_path).unlink(missing_ok=True)


def test_parse_full_def_xml(def_xml_file):
    """Test parsing the complete servalcat_pipe DEF XML."""
    task = parse_def_xml_file(def_xml_file)

    assert task is not None
    assert task.name == 'servalcat_pipe'
    assert type(task).__name__ == 'CContainer'


def test_input_data_container(def_xml_file):
    """Test inputData container structure."""
    task = parse_def_xml_file(def_xml_file)

    assert hasattr(task, 'inputData')
    input_data = task.inputData

    # Test XYZIN
    assert hasattr(input_data, 'XYZIN')
    assert type(input_data.XYZIN).__name__ == 'CPdbDataFile'

    # Test HKLIN
    assert hasattr(input_data, 'HKLIN')
    assert type(input_data.HKLIN).__name__ == 'CObsDataFile'

    # Test DICT_LIST
    assert hasattr(input_data, 'DICT_LIST')
    assert type(input_data.DICT_LIST).__name__ == 'CList'


def test_output_data_container(def_xml_file):
    """Test outputData container structure."""
    task = parse_def_xml_file(def_xml_file)

    assert hasattr(task, 'outputData')
    output_data = task.outputData

    # Test XYZOUT
    assert hasattr(output_data, 'XYZOUT')
    assert type(output_data.XYZOUT).__name__ == 'CPdbDataFile'

    # Test FPHIOUT
    assert hasattr(output_data, 'FPHIOUT')
    assert type(output_data.FPHIOUT).__name__ == 'CMapCoeffsDataFile'


def test_control_parameters_container(def_xml_file):
    """Test controlParameters container with various types."""
    task = parse_def_xml_file(def_xml_file)

    assert hasattr(task, 'controlParameters')
    ctrl = task.controlParameters

    # Test CString with enumerators
    assert hasattr(ctrl, 'DATA_METHOD')
    assert type(ctrl.DATA_METHOD).__name__ == 'CString'
    assert ctrl.DATA_METHOD.value == 'xtal'

    # Test CBoolean
    assert hasattr(ctrl, 'ADD_WATERS')
    assert type(ctrl.ADD_WATERS).__name__ == 'CBoolean'
    assert ctrl.ADD_WATERS.value is False

    # Test CInt with constraints
    assert hasattr(ctrl, 'NCYCLES')
    assert type(ctrl.NCYCLES).__name__ == 'CInt'
    assert ctrl.NCYCLES.value == 10

    # Test CFloat
    assert hasattr(ctrl, 'WEIGHT')
    assert type(ctrl.WEIGHT).__name__ == 'CFloat'

    # Test another CString
    assert hasattr(ctrl, 'B_REFINEMENT_MODE')
    assert ctrl.B_REFINEMENT_MODE.value == 'iso'

    # Test another CBoolean
    assert hasattr(ctrl, 'OCCUPANCY_REFINEMENT')
    assert ctrl.OCCUPANCY_REFINEMENT.value is True


def test_metal_coord_pipeline_container(def_xml_file):
    """Test metalCoordPipeline container at root level."""
    task = parse_def_xml_file(def_xml_file)

    assert hasattr(task, 'metalCoordPipeline')
    metal = task.metalCoordPipeline

    # Test RUN_METALCOORD
    assert hasattr(metal, 'RUN_METALCOORD')
    assert type(metal.RUN_METALCOORD).__name__ == 'CBoolean'
    assert metal.RUN_METALCOORD.value is False

    # Test LINKS
    assert hasattr(metal, 'LINKS')
    assert type(metal.LINKS).__name__ == 'CString'
    assert metal.LINKS.value == 'UPDATE'


def test_default_values_applied(def_xml_file):
    """Test that default values from qualifiers are applied."""
    task = parse_def_xml_file(def_xml_file)

    # Check various defaults
    assert task.controlParameters.DATA_METHOD.value == 'xtal'
    assert task.controlParameters.ADD_WATERS.value is False
    assert task.controlParameters.NCYCLES.value == 10
    assert task.controlParameters.B_REFINEMENT_MODE.value == 'iso'
    assert task.controlParameters.OCCUPANCY_REFINEMENT.value is True
    assert task.metalCoordPipeline.RUN_METALCOORD.value is False
    assert task.metalCoordPipeline.LINKS.value == 'UPDATE'


def test_value_modification(def_xml_file):
    """Test modifying parameter values."""
    task = parse_def_xml_file(def_xml_file)

    ctrl = task.controlParameters

    # Modify CInt
    original_cycles = ctrl.NCYCLES.value
    assert original_cycles == 10

    ctrl.NCYCLES.value = 25
    assert ctrl.NCYCLES.value == 25

    # Modify CBoolean
    original_waters = ctrl.ADD_WATERS.value
    assert original_waters is False

    ctrl.ADD_WATERS.value = True
    assert ctrl.ADD_WATERS.value is True

    # Modify CString
    original_mode = ctrl.B_REFINEMENT_MODE.value
    assert original_mode == 'iso'

    ctrl.B_REFINEMENT_MODE.value = 'aniso'
    assert ctrl.B_REFINEMENT_MODE.value == 'aniso'


def test_min_max_constraints(def_xml_file):
    """Test that min/max constraints are validated."""
    task = parse_def_xml_file(def_xml_file)

    ncycles = task.controlParameters.NCYCLES

    # Should have min=0
    min_val = ncycles.get_qualifier('min')
    assert min_val == 0

    # Setting below min should raise ValueError
    with pytest.raises(ValueError):
        ncycles.value = -5


def test_structure_statistics(def_xml_file):
    """Test that task has expected containers and content."""
    task = parse_def_xml_file(def_xml_file)

    # Test containers exist
    assert hasattr(task, 'inputData')
    assert hasattr(task, 'outputData')
    assert hasattr(task, 'controlParameters')
    assert hasattr(task, 'metalCoordPipeline')

    # Test inputData contents
    assert hasattr(task.inputData, 'XYZIN')
    assert hasattr(task.inputData, 'HKLIN')
    assert hasattr(task.inputData, 'DICT_LIST')

    # Test outputData contents
    assert hasattr(task.outputData, 'XYZOUT')
    assert hasattr(task.outputData, 'FPHIOUT')

    # Test controlParameters contents
    assert hasattr(task.controlParameters, 'DATA_METHOD')
    assert hasattr(task.controlParameters, 'ADD_WATERS')
    assert hasattr(task.controlParameters, 'NCYCLES')
    assert hasattr(task.controlParameters, 'WEIGHT')
    assert hasattr(task.controlParameters, 'B_REFINEMENT_MODE')
    assert hasattr(task.controlParameters, 'OCCUPANCY_REFINEMENT')

    # Test metalCoordPipeline contents
    assert hasattr(task.metalCoordPipeline, 'RUN_METALCOORD')
    assert hasattr(task.metalCoordPipeline, 'LINKS')
