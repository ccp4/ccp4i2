import React, { useCallback, useMemo } from "react";
import { Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { FieldRow } from "../task-elements/field-row";
import { useJob } from "../../../utils";
import {
  CCP4i2ErrorReport,
  useFreeRWarning,
} from "../../../providers/run-check-provider";

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
  const { useTaskItem, createPeerTask, validation } = useJob(job.id);

  // Get task items
  const { value: HKLINValue } = useTaskItem("servalcat_pipe.container.inputData.HKLIN");
  const { value: MAP_SHARP } = useTaskItem("MAP_SHARP");
  const { value: MAP_SHARP_CUSTOM } = useTaskItem("MAP_SHARP_CUSTOM");
  const { value: freeRFlag } = useTaskItem("FREERFLAG");

  // Filter out metalCoordWrapper.inputData.XYZIN validation errors
  // TODO: This should be moved to Python's validity() method in servalcat_pipe.py
  const filterServalcatErrors = useCallback(
    (errors: CCP4i2ErrorReport): CCP4i2ErrorReport => {
      return Object.fromEntries(
        Object.entries(errors).filter(
          ([key, _]) => key !== "servalcat_pipe.metalCoordWrapper.inputData.XYZIN"
        )
      );
    },
    []
  );

  // Use centralized FreeR warning hook
  useFreeRWarning({
    job,
    taskName: "servalcat_pipe",
    freeRFlag,
    validation,
    createPeerTask,
    filterErrors: filterServalcatErrors,
  });

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
                itemName="VALIDATE_IRIS"
                qualifiers={{
                  guiLabel: "Analyse using iris",
                  toolTip: "Generate using iris validation tool",
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
