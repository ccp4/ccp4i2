import React, {
  useCallback,
  useContext,
  useEffect,
  useMemo,
  useRef,
} from "react";
import { Button, Paper, Typography } from "@mui/material";
import { useRouter } from "next/navigation";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { FieldRow } from "../task-elements/field-row";
import { useJob } from "../../../utils";
import {
  CCP4i2ErrorReport,
  useRunCheck,
} from "../../../providers/run-check-provider";
import { Job } from "../../../types/models";

/**
 * Task interface component for ServalCat Pipe - Macromolecular Refinement Pipeline.
 *
 * ServalCat Pipe is used for:
 * - Automated macromolecular structure refinement
 * - Integration with geometry validation and analysis
 * - Map calculation with optional sharpening
 * - NCS and twinning handling
 * - Comprehensive validation reporting (MolProbity, Ramachandran, B-factors)
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const router = useRouter();
  const { useTaskItem, createPeerTask, validation } = useJob(job.id);

  // Use refs to track processed states and prevent cycles
  const lastProcessedValidation = useRef<any>(null);
  const lastProcessedFreeRFlag = useRef<any>(null);

  // Get task items
  const { value: HKLINValue } = useTaskItem("servalcat_pipe.container.inputData.HKLIN");
  const { value: MAP_SHARP } = useTaskItem("MAP_SHARP");
  const { value: MAP_SHARP_CUSTOM } = useTaskItem("MAP_SHARP_CUSTOM");
  const { value: freeRFlag } = useTaskItem("FREERFLAG");

  // Context for error handling
  const {
    processedErrors,
    setProcessedErrors,
    extraDialogActions,
    setExtraDialogActions,
    setRunTaskRequested,
  } = useRunCheck();

  // Derived state (memoized for performance)
  const intensitiesAvailable = useMemo(() => {
    return [1, 3].includes(HKLINValue?.contentFlag);
  }, [HKLINValue?.contentFlag]);

  // Visibility conditions (stable references)
  const visibility = useMemo(
    () => ({
      showMapSharpening: () => MAP_SHARP === true,
      showCustomBFactor: () => MAP_SHARP === true && MAP_SHARP_CUSTOM === true,
      showAmplitudesWarning: () => !intensitiesAvailable,
      showRefinementChoice: () => intensitiesAvailable,
    }),
    [MAP_SHARP, MAP_SHARP_CUSTOM, intensitiesAvailable]
  );

  // Stable task creation function
  const createFreeRTask = useCallback(async () => {
    try {
      const createdJob: Job | undefined = await createPeerTask("freerflag");
      if (createdJob) {
        router.push(`/project/${job.project}/job/${createdJob.id}`);
        setRunTaskRequested(null);
      }
    } catch (error) {
      console.error("Error creating FreeR task:", error);
    }
  }, [createPeerTask, job.project, router, setRunTaskRequested]);

  // Process validation errors with cycle prevention
  const processErrors = useCallback(() => {
    if (!validation) return;

    // Filter out specific validation errors
    const newProcessedErrors = Object.fromEntries(
      Object.entries(validation as CCP4i2ErrorReport).filter(
        ([key, _]) => key !== "servalcat_pipe.metalCoordWrapper.inputData.XYZIN"
      )
    );

    // Add FreeR flag warning if not set
    if (!freeRFlag?.dbFileId?.length) {
      newProcessedErrors["servalcat_pipe.container.inputData.FREERFLAG"] = {
        messages: [
          "Setting the Free R flag file is strongly recommended for refinement",
          "You are advised to select an existing set or create a new one",
        ],
        maxSeverity: 3, // Allows execution but shows warning
      };
    }

    // Only update if errors have actually changed
    const newErrorsKey = JSON.stringify(newProcessedErrors);
    const currentErrorsKey = JSON.stringify(processedErrors);

    if (newErrorsKey !== currentErrorsKey) {
      setProcessedErrors(newProcessedErrors);
    }
  }, [validation, freeRFlag, processedErrors, setProcessedErrors]);

  // Handle extra dialog actions for FreeR flag
  const updateExtraDialogActions = useCallback(() => {
    if (!freeRFlag?.dbFileId?.length) {
      // Only add action if it doesn't already exist
      if (!extraDialogActions?.FREERFLAG) {
        const newExtraDialogActions = {
          ...extraDialogActions,
          FREERFLAG: (
            <Button variant="contained" onClick={createFreeRTask}>
              Create FreeR task
            </Button>
          ),
        };
        setExtraDialogActions(newExtraDialogActions);
      }
    } else {
      // Remove action if FreeR flag is now set
      if (extraDialogActions?.FREERFLAG) {
        const { FREERFLAG, ...remainingActions } = extraDialogActions;
        setExtraDialogActions(
          Object.keys(remainingActions).length > 0 ? remainingActions : null
        );
      }
    }
  }, [freeRFlag, extraDialogActions, setExtraDialogActions, createFreeRTask]);

  // Effect for error processing with minimal dependencies
  useEffect(() => {
    if (freeRFlag !== undefined) {
      processErrors();
    }
  }, [freeRFlag, processErrors]);

  // Effect for handling extra dialog actions
  useEffect(() => {
    updateExtraDialogActions();
  }, [updateExtraDialogActions]);

  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab label="Input data" key="input">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Main inputs",
              initiallyOpen: true,
            }}
            containerHint="BlockLevel"
          >
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{
                guiLabel: "Structure coordinates",
                initiallyOpen: true,
              }}
              containerHint="BlockLevel"
              sx={{
                border: "3px solid grey",
                borderRadius: "0.5rem",
                p: 2,
                mb: 2,
              }}
            >
              <CCP4i2TaskElement
                {...props}
                itemName="servalcat_pipe.container.inputData.XYZIN"
                qualifiers={{
                  guiLabel: "Input coordinates",
                  toolTip: "Macromolecular coordinates for refinement",
                }}
              />
            </CCP4i2ContainerElement>

            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{
                guiLabel: "Reflection data",
                initiallyOpen: true,
              }}
              containerHint="BlockLevel"
              sx={{
                border: "3px solid grey",
                borderRadius: "0.5rem",
                p: 2,
                mb: 2,
              }}
            >
              <CCP4i2TaskElement
                {...props}
                itemName="HKLIN"
                qualifiers={{
                  guiLabel: "Reflection data",
                  toolTip: "Observed reflection data for refinement",
                }}
              />

              {visibility.showRefinementChoice() ? (
                <CCP4i2TaskElement
                  {...props}
                  itemName="F_SIGF_OR_I_SIGI"
                  qualifiers={{
                    guiLabel: "Refinement against",
                    toolTip: "Choose between structure factors or intensities",
                  }}
                />
              ) : (
                <Typography
                  variant="body1"
                  sx={{ fontWeight: "medium", mt: 1 }}
                >
                  Using <strong>amplitudes</strong>
                </Typography>
              )}

              <CCP4i2TaskElement
                {...props}
                itemName="FREERFLAG"
                qualifiers={{
                  guiLabel: "Free R flags",
                  toolTip: "Test set flags for cross-validation",
                }}
              />
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Additional geometry dictionaries",
              initiallyOpen: true,
            }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="DICT_LIST"
              qualifiers={{
                guiLabel: "Dictionaries",
                toolTip:
                  "Additional geometry dictionaries for non-standard residues",
              }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="Output" key="output">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Refinement options",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="USE_NCS"
              qualifiers={{
                guiLabel: "Use NCS if present",
                toolTip: "Apply non-crystallographic symmetry restraints",
              }}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="USE_TWIN"
              qualifiers={{
                guiLabel: "Use twinning",
                toolTip: "Handle crystal twinning during refinement",
              }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Output options",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="OUTPUT_HYDROGENS"
              qualifiers={{
                guiLabel: "Output calculated riding hydrogens to file",
                toolTip:
                  "Include calculated hydrogen atoms in output coordinates",
              }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Map calculation",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="MAP_SHARP"
              qualifiers={{
                guiLabel: "Perform map sharpening when calculating maps",
                toolTip: "Apply B-factor sharpening to improve map quality",
              }}
            />

            <FieldRow sx={{ mt: 1 }}>
              <CCP4i2TaskElement
                {...props}
                itemName="MAP_SHARP_CUSTOM"
                qualifiers={{
                  guiLabel: "Use custom sharpening parameter (B-factor)",
                  toolTip: "Specify custom B-factor for map sharpening",
                }}
                visibility={visibility.showMapSharpening}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="BSHARP"
                qualifiers={{
                  guiLabel: "B factor to use",
                  toolTip: "B-factor value for custom map sharpening",
                }}
                visibility={visibility.showCustomBFactor}
              />
            </FieldRow>
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Validation and analysis",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            <FieldRow>
              <CCP4i2TaskElement
                {...props}
                itemName="VALIDATE_BAVERAGE"
                qualifiers={{
                  guiLabel: "Analyse B-factor distributions",
                  toolTip: "Generate B-factor distribution analysis",
                }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="VALIDATE_RAMACHANDRAN"
                qualifiers={{
                  guiLabel: "Calculate Ramachandran plots",
                  toolTip: "Generate Ramachandran plot validation",
                }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="VALIDATE_MOLPROBITY"
                qualifiers={{
                  guiLabel: "Run MolProbity to analyse geometry",
                  toolTip:
                    "Perform comprehensive geometry validation with MolProbity",
                }}
              />
            </FieldRow>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
