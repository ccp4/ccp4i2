import React, { useCallback, useEffect, useMemo, useRef } from "react";
import { Box, IconButton, Paper, Tooltip } from "@mui/material";
import { Download } from "@mui/icons-material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";
import {
  useRunCheck,
  useSequenceWarning,
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
  const { useTaskItem, useFileDigest, fetchDigest, mutateContainer, validation, createPeerTask } = useJob(
    job.id
  );
  const { setProcessedErrors } = useRunCheck();

  // Get SEQIN for sequence warning
  const { value: SEQIN } = useTaskItem("SEQIN");

  // Use centralized sequence warning hook
  useSequenceWarning({
    job,
    taskName: "shelx",
    sequence: SEQIN,
    validation,
    createPeerTask,
  });

  // Refs for preventing cycles
  const initializationDone = useRef(false);

  // Get task items for file handling and parameter updates
  const { item: ATOM_TYPEItem, value: ATOM_TYPEValue } = useTaskItem("ATOM_TYPE");
  const { item: F_SIGFanomItem, value: F_SIGFanomValue } = useTaskItem("F_SIGFanom");
  const { forceUpdate: forceUpdateWAVELENGTH } = useTaskItem("WAVELENGTH");
  const { forceUpdate: forceUpdateUSER_WAVELENGTH } = useTaskItem("USER_WAVELENGTH");
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

  // File digest for wavelength extraction - only fetch when file has been uploaded
  const hasF_SIGFanom = Boolean(F_SIGFanomValue?.dbFileId);
  const f_sigfanomDigestPath = hasF_SIGFanom && F_SIGFanomItem?._objectPath ? F_SIGFanomItem._objectPath : "";
  const { data: F_SIGFanomDigest } = useFileDigest(f_sigfanomDigestPath);

  // Wavelength extraction handler for F_SIGFanom onChange (imperative)
  // Uses fetchDigest to get the digest for the newly selected file
  const handleF_SIGFanomChange = useCallback(async () => {
    if (!forceUpdateWAVELENGTH || !F_SIGFanomItem?._objectPath || !job || job.status !== 1)
      return;

    // Fetch digest imperatively for the new file
    const digestData = await fetchDigest(F_SIGFanomItem._objectPath);

    // Extract wavelength from digest
    if (digestData?.wavelengths?.length > 0) {
      const newWavelength =
        digestData.wavelengths[digestData.wavelengths.length - 1];

      // Only update if wavelength is valid
      if (newWavelength && newWavelength < 9) {
        try {
          // Use forceUpdate to ensure UI reflects the new value immediately
          await forceUpdateWAVELENGTH(newWavelength);
          // Set USER_WAVELENGTH=true so crank2 knows to use it
          if (forceUpdateUSER_WAVELENGTH) {
            await forceUpdateUSER_WAVELENGTH(true);
          }
          mutateContainer();
        } catch (error) {
          console.error("Error updating wavelength:", error);
        }
      }
    }
  }, [forceUpdateWAVELENGTH, forceUpdateUSER_WAVELENGTH, F_SIGFanomItem?._objectPath, fetchDigest, job?.status, mutateContainer]);

  // Handler for fetching wavelength from reflection file (button click)
  const handleFetchWavelength = useCallback(async () => {
    if (!forceUpdateWAVELENGTH || !F_SIGFanomDigest || !job || job.status !== 1) {
      return;
    }

    const digestData = F_SIGFanomDigest;

    if (digestData?.wavelengths?.length > 0) {
      const newWavelength =
        digestData.wavelengths[digestData.wavelengths.length - 1];

      if (newWavelength && newWavelength < 9) {
        try {
          // Use forceUpdate to ensure UI reflects the new value immediately
          await forceUpdateWAVELENGTH(newWavelength);
          if (forceUpdateUSER_WAVELENGTH) {
            await forceUpdateUSER_WAVELENGTH(true);
          }
          mutateContainer();
        } catch (error) {
          console.error("Error updating wavelength:", error);
        }
      }
    }
  }, [forceUpdateWAVELENGTH, forceUpdateUSER_WAVELENGTH, F_SIGFanomDigest, job?.status, mutateContainer]);

  // Check if we can fetch wavelength (file is set and has wavelength data)
  const canFetchWavelength = hasF_SIGFanom && F_SIGFanomDigest?.wavelengths?.length > 0 && job.status === 1;

  // Element configurations for parameters section
  const parameterConfigs = useMemo(
    () => [
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
    []
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
    if (!setProcessedErrors || !ATOM_TYPEItem?._objectPath) return;
    if (!ATOM_TYPEValue || ATOM_TYPEValue.trim() === "") {
      // ATOM_TYPE is empty - add error using the actual item's objectPath
      const errors = {
        ...(validation || {}),
        [ATOM_TYPEItem._objectPath]: {
          maxSeverity: 2,
          messages: ["ATOM_TYPE is required"],
        },
      };
      setProcessedErrors(errors);
    } else {
      // ATOM_TYPE is valid - clear custom errors (revert to server validation)
      setProcessedErrors(null);
    }
  }, [validation, setProcessedErrors, ATOM_TYPEItem?._objectPath, ATOM_TYPEValue]);

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
            <CCP4i2TaskElement
              {...props}
              itemName="F_SIGFanom"
              qualifiers={{
                guiLabel: "Reflections",
                toolTip: "Anomalous reflection data for phasing",
              }}
              onChange={handleF_SIGFanomChange}
            />
            <Box sx={{ display: "flex", alignItems: "flex-start", gap: 0.5 }}>
              {canFetchWavelength && (
                <Tooltip title="Fetch wavelength from reflection file">
                  <IconButton
                    size="small"
                    onClick={handleFetchWavelength}
                    color="primary"
                    sx={{ mt: 1 }}
                  >
                    <Download fontSize="small" />
                  </IconButton>
                </Tooltip>
              )}
              <Box sx={{ flex: 1 }}>
                <CCP4i2TaskElement
                  {...props}
                  itemName="WAVELENGTH"
                  qualifiers={{
                    guiLabel: "Wavelength",
                    toolTip: "X-ray wavelength used for data collection",
                  }}
                />
              </Box>
            </Box>
            <CCP4i2TaskElement
              {...props}
              itemName="SEQIN"
              qualifiers={{
                guiLabel: "Asymmetric unit content",
                toolTip: "Sequence file defining the protein content",
              }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="FREERFLAG"
              qualifiers={{
                guiLabel: "Free R flags",
                toolTip: "Test set flags for cross-validation",
              }}
            />
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
            {renderElements(parameterConfigs)}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
