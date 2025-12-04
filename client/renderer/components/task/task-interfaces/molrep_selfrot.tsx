import React, { useCallback, useMemo } from "react";
import { Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";

/**
 * Task interface component for MOLREP Self Rotation Function calculation.
 *
 * Provides functionality for:
 * - Self rotation function calculation using structure factor data
 * - Patterson function analysis for non-crystallographic symmetry detection
 * - Resolution and angular sampling parameter control
 * - Advanced options for rotation function optimization
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem, fetchDigest } = useJob(job.id);

  // Get task items for file handling and parameter updates
  const { item: F_SIGFItem } = useTaskItem("F_SIGF");
  const { update: updateWAVELENGTH } = useTaskItem("WAVELENGTH");

  // Handle F_SIGF file change - extract wavelength from digest
  const handleF_SIGFChange = useCallback(async () => {
    if (!F_SIGFItem?._objectPath) return;

    const digestData = await fetchDigest(F_SIGFItem._objectPath);
    const wavelength = digestData?.wavelengths?.at(-1);
    if (wavelength && wavelength > 0 && wavelength < 9) {
      await updateWAVELENGTH(wavelength);
    }
  }, [F_SIGFItem?._objectPath, fetchDigest, updateWAVELENGTH]);

  // Task values for visibility conditions
  const taskValues = useMemo(
    () => ({
      sculptorMode: useTaskItem("SEQ").value,
      useAdvanced: useTaskItem("USE_ADVANCED").value,
    }),
    [useTaskItem]
  );

  // Visibility conditions
  const visibility = useMemo(
    () => ({
      hasAdvancedOptions: () => taskValues.useAdvanced,
      isSculptorEnabled: () => taskValues.sculptorMode !== "n",
    }),
    [taskValues]
  );

  // Element configurations
  const elementConfigs = useMemo(
    () => ({
      inputData: [{ key: "F_SIGF", label: "Reflections", onChange: handleF_SIGFChange }],
      basicOptions: [
        {
          key: "SEQ",
          label:
            "Perform alignment and use it to rename residues and trim side chains",
          guiMode: "radio",
          menuText: ["always", "only for sequence identity > 20%", "never"],
          enumerators: ["y", "d", "n"],
        },
      ],
      rotationFunction: [
        { key: "RESO_LOW", label: "Low resolution limit" },
        { key: "RESO_HIGH", label: "High resolution limit" },
        { key: "ANGLE_STEP", label: "Angular step size" },
        { key: "RADIUS", label: "Integration radius" },
      ],
      advancedOptions: [
        { key: "USE_ADVANCED", label: "Use advanced options" },
        {
          key: "NPEAKS",
          label: "Number of peaks to find",
          visible: visibility.hasAdvancedOptions,
        },
        {
          key: "THETA_MIN",
          label: "Minimum theta angle",
          visible: visibility.hasAdvancedOptions,
        },
        {
          key: "THETA_MAX",
          label: "Maximum theta angle",
          visible: visibility.hasAdvancedOptions,
        },
        {
          key: "PHI_STEP",
          label: "Phi step size",
          visible: visibility.hasAdvancedOptions,
        },
      ],
      keywords: [{ key: "keywords", label: "Additional keywords" }],
    }),
    [visibility, handleF_SIGFChange]
  );

  // Render helper function
  const renderElements = useCallback(
    (elements: any[]) =>
      elements.map(({ key, label, visible = () => true, onChange, ...extraProps }) => (
        <CCP4i2TaskElement
          {...props}
          key={key}
          itemName={key}
          qualifiers={{ guiLabel: label, ...extraProps }}
          visibility={visible}
          onChange={onChange}
        />
      )),
    [props]
  );

  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab label="Main inputs" key="main">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Input data",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
            key="Input data"
          >
            {renderElements(elementConfigs.inputData)}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="Basic options" key="basic">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel:
                "Perform alignment and use it to rename residues and trim side chains",
            }}
            containerHint="BlockLevel"
            key="Sculptor"
          >
            {renderElements(elementConfigs.basicOptions)}
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Rotation function parameters" }}
            containerHint="FolderLevel"
            key="Rotation function"
          >
            {renderElements(elementConfigs.rotationFunction)}
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Advanced options" }}
            containerHint="FolderLevel"
            key="Advanced"
          >
            {renderElements(elementConfigs.advancedOptions)}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="Keywords" key="keywords">
          {renderElements(elementConfigs.keywords)}
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
