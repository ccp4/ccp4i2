import React, { useCallback, useEffect, useMemo, useRef } from "react";
import { Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";
import {
  CCP4i2ErrorReport,
  useRunCheck,
} from "../../../providers/run-check-provider";

/**
 * Task interface component for SHELX - Experimental Phasing Pipeline.
 *
 * SHELX is used for:
 * - Single wavelength anomalous dispersion (SAD) phasing
 * - Multi-wavelength anomalous dispersion (MAD) phasing
 * - Heavy atom substructure determination
 * - Phase calculation and density modification
 * - Automated model building with Buccaneer integration
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem, useFileDigest, mutateContainer, validation } = useJob(
    job.id
  );
  const { setProcessedErrors } = useRunCheck();

  // Refs for preventing cycles
  const initializationDone = useRef(false);

  // Get task items for file handling and parameter updates
  const { value: ATOM_TYPEValue } = useTaskItem("ATOM_TYPE");
  const { item: F_SIGFanomItem } = useTaskItem("F_SIGFanom");
  const { updateNoMutate: updateWAVELENGTH } = useTaskItem("WAVELENGTH");
  const { updateNoMutate: updateSHELXCDE } = useTaskItem("SHELXCDE");
  const { updateNoMutate: updateUSE_COMB } = useTaskItem("USE_COMB");
  const { updateNoMutate: updateSHELX_SEPAR } = useTaskItem("SHELX_SEPAR");
  const { updateNoMutate: updateMB_PROGRAM } = useTaskItem("MB_PROGRAM");

  // Get current values for initial setup
  const taskValues = useMemo(
    () => ({
      SHELXCDE: useTaskItem("SHELXCDE").value,
      USE_COMB: useTaskItem("USE_COMB").value,
      SHELX_SEPAR: useTaskItem("SHELX_SEPAR").value,
      MB_PROGRAM: useTaskItem("MB_PROGRAM").value,
    }),
    [useTaskItem]
  );

  // File digest for wavelength extraction
  const { data: F_SIGFanomDigest } = useFileDigest(F_SIGFanomItem?._objectPath);

  // Wavelength extraction handler for F_SIGFanom onChange
  const handleF_SIGFanomChange = useCallback(async () => {
    if (!updateWAVELENGTH || !F_SIGFanomDigest || !job || job.status !== 1)
      return;

    // F_SIGFanomDigest is now unwrapped by useFileDigest
    const digestData = F_SIGFanomDigest;

    // Extract wavelength from digest
    if (digestData?.wavelengths?.length > 0) {
      const newWavelength =
        digestData.wavelengths[digestData.wavelengths.length - 1];

      // Only update if wavelength is valid
      if (newWavelength && newWavelength < 9) {
        try {
          console.log(
            `Extracting wavelength from F_SIGFanom: ${newWavelength}`
          );
          await updateWAVELENGTH(newWavelength);
          await mutateContainer();
        } catch (error) {
          console.error("Error updating wavelength:", error);
        }
      }
    }
  }, [updateWAVELENGTH, F_SIGFanomDigest, job?.status, mutateContainer]);

  // Element configurations
  const elementConfigs = useMemo(
    () => ({
      keyFiles: [
        {
          key: "F_SIGFanom",
          label: "Reflections",
          toolTip: "Anomalous reflection data for phasing",
          onChange: handleF_SIGFanomChange,
        },
        {
          key: "WAVELENGTH",
          label: "Wavelength",
          toolTip: "X-ray wavelength used for data collection",
        },
        {
          key: "SEQIN",
          label: "Asymmetric unit content",
          toolTip: "Sequence file defining the protein content",
        },
        {
          key: "FREERFLAG",
          label: "Free R flags",
          toolTip: "Test set flags for cross-validation",
        },
      ],
      parameters: [
        {
          key: "ATOM_TYPE",
          label: "Anomalous atom type",
          toolTip: "Type of heavy atom providing anomalous signal",
        },
        {
          key: "START_PIPELINE",
          label: "First step for analysis",
          toolTip: "Starting point in the SHELX pipeline",
        },
        {
          key: "END_PIPELINE",
          label: "Last step for analysis",
          toolTip: "Ending point in the SHELX pipeline",
        },
      ],
    }),
    [handleF_SIGFanomChange]
  );

  // Reset initialization flag when job changes
  useEffect(() => {
    initializationDone.current = false;
  }, [job?.id]);

  // Process validation errors - add custom ATOM_TYPE validation
  // TODO: This client-side validation should be moved to Python's validity() method
  // in pipelines/crank2/script/crank2.py (SHELX uses crank2 under the hood).
  // Python validity() can check ATOM_TYPE and add appropriate error messages server-side.
  useEffect(() => {
    if (!validation || !setProcessedErrors) return;
    if (!ATOM_TYPEValue || ATOM_TYPEValue.trim() === "") {
      const processedErrors = {
        ...validation,
        "crank2.container.inputData.ATOM_TYPE": {
          maxSeverity: 2,
          messages: ["ATOM_TYPE is required"],
        },
      };
      setProcessedErrors(processedErrors);
    }
  }, [job, validation, setProcessedErrors, ATOM_TYPEValue]);

  // Initialize defaults once when job becomes editable
  // Use a direct effect with minimal dependencies to avoid cycles
  useEffect(() => {
    const initializeDefaults = async () => {
      if (initializationDone.current || !job || job.status !== 1) return;

      const updates: Promise<any>[] = [];

      // Set default values only if they're not already set
      if (taskValues.SHELXCDE === undefined || taskValues.SHELXCDE === null) {
        updates.push(updateSHELXCDE(true));
      }
      if (taskValues.USE_COMB === true) {
        updates.push(updateUSE_COMB(false));
      }
      if (
        taskValues.SHELX_SEPAR === undefined ||
        taskValues.SHELX_SEPAR === null
      ) {
        updates.push(updateSHELX_SEPAR(true));
      }
      if (taskValues.MB_PROGRAM !== "buccaneer") {
        updates.push(updateMB_PROGRAM("buccaneer"));
      }

      if (updates.length > 0) {
        try {
          await Promise.all(updates);
          await mutateContainer();
        } catch (error) {
          console.error("Error initializing SHELX defaults:", error);
        }
      }

      initializationDone.current = true;
    };

    // Only run when job status becomes editable and we haven't initialized yet
    if (job?.status === 1 && !initializationDone.current) {
      initializeDefaults();
    }
  }, [
    job?.status,
    // Only depend on the actual values, not the update functions
    taskValues.SHELXCDE,
    taskValues.USE_COMB,
    taskValues.SHELX_SEPAR,
    taskValues.MB_PROGRAM,
  ]);

  // Render helper function
  const renderElements = useCallback(
    (elements: any[]) =>
      elements.map(({ key, label, onChange, ...extraProps }) => (
        <CCP4i2TaskElement
          {...props}
          key={key}
          itemName={key}
          qualifiers={{ guiLabel: label, ...extraProps }}
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
              guiLabel: "Key files",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            {renderElements(elementConfigs.keyFiles)}
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Parameters",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            {renderElements(elementConfigs.parameters)}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
