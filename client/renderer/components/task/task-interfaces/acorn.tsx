import React from "react";
import { Box, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";
import { useBoolToggle } from "../task-elements/shared-hooks";
import { InlineField } from "../task-elements/inline-field";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);

  // Input Data tab
  const { value: ACORN_PHSIN_TYPE } = useTaskItem("ACORN_PHSIN_TYPE");

  // Advanced tab - boolean toggles
  const customPhase = useBoolToggle(useTaskItem, "ACOPH_CUSTOM");
  const custddm = useBoolToggle(useTaskItem, "ACOPH_CUSTDDM");
  const becut = useBoolToggle(useTaskItem, "ACORN_BECUT");
  const bresol = useBoolToggle(useTaskItem, "ACORN_BRESOL");
  const bexclude = useBoolToggle(useTaskItem, "ACORN_BEXCLUDE");
  const bgrid = useBoolToggle(useTaskItem, "ACORN_BGRID");

  const { value: ACOPH_TRIALS } = useTaskItem("ACOPH_TRIALS");
  const numTrials =
    typeof ACOPH_TRIALS === "number"
      ? ACOPH_TRIALS
      : parseInt(ACOPH_TRIALS as string) || 1;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs {...props}>
        {/* ===== Tab 1: Input Data ===== */}
        <CCP4i2Tab label="Input Data">
          <InlineField label="Run ACORN with" width="auto" sx={{ mb: 1 }}>
            <CCP4i2TaskElement
              itemName="ACORN_PHSIN_TYPE"
              {...props}
              qualifiers={{ guiLabel: " ", guiMode: "radio" }}
            />
          </InlineField>

          {/* Reflection Data */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Reflection Data" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="F_SIGF"
              {...props}
              qualifiers={{ guiLabel: "Reflections" }}
            />
            <CCP4i2TaskElement
              itemName="ABCD"
              {...props}
              qualifiers={{ guiLabel: "Phases" }}
              visibility={() => ACORN_PHSIN_TYPE === "phases"}
            />
          </CCP4i2ContainerElement>

          {/* Model for approximate co-ordinates (model mode only) */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Model for approximate co-ordinates" }}
            containerHint="FolderLevel"
            visibility={() => ACORN_PHSIN_TYPE === "model"}
          >
            <CCP4i2TaskElement
              itemName="XYZIN"
              {...props}
              qualifiers={{ guiLabel: "Atomic model" }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ===== Tab 2: Advanced Acorn Parameters ===== */}
        <CCP4i2Tab label="Advanced Acorn Parameters">
          {/* General ACORN phase improvement parameters */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "General ACORN phase improvement parameters" }}
            containerHint="BlockLevel"
          >
            <InlineField label="Number of trials:">
              <CCP4i2TaskElement
                itemName="ACOPH_TRIALS"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <CCP4i2TaskElement
              itemName="ACOPH_CUSTOM"
              {...props}
              qualifiers={{ guiLabel: "Define custom parameters for each trial" }}
              onChange={customPhase.onChange}
            />
            {customPhase.value && (
              <Box sx={{ pl: 3, display: "flex", flexDirection: "column", gap: 1 }}>
                <CCP4i2TaskElement
                  itemName="ACOPH_CUSTDDM"
                  {...props}
                  qualifiers={{ guiLabel: "Define DDM type for each trial" }}
                  onChange={custddm.onChange}
                />
                {Array.from({ length: numTrials }, (_, i) => i + 1).map((n) => (
                  <Box
                    key={n}
                    sx={{ display: "flex", alignItems: "center", gap: 1, flexWrap: "wrap" }}
                  >
                    <Typography variant="body2" sx={{ fontWeight: "bold", minWidth: "4rem" }}>
                      Trial {n}
                    </Typography>
                    <Box sx={{ width: "6rem" }}>
                      <CCP4i2TaskElement itemName={`ACOPH_NCDDM_${n}`} {...props} qualifiers={{ guiLabel: " " }} />
                    </Box>
                    <Typography variant="body2">cycles of</Typography>
                    <Box sx={{ width: "6rem" }}>
                      <CCP4i2TaskElement itemName={`ACOPH_DDMK_${n}`} {...props} qualifiers={{ guiLabel: " " }} />
                    </Box>
                    <Box sx={{ width: "14rem" }}>
                      <CCP4i2TaskElement itemName={`ACOPH_REFINE_${n}`} {...props} qualifiers={{ guiLabel: " " }} />
                    </Box>
                  </Box>
                ))}
              </Box>
            )}
            <InlineField label="Cease DDM cycling if phase shift between consecutive cycles is less than">
              <CCP4i2TaskElement
                itemName="ACOPH_PSFINISH"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
          </CCP4i2ContainerElement>

          {/* Selection of Reflection Data */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Selection of Reflection Data" }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              itemName="ACORN_BRESOL"
              {...props}
              qualifiers={{ guiLabel: "User defined resolution range for reflection data" }}
              onChange={bresol.onChange}
            />
            {bresol.value && (
              <Box sx={{ pl: 3, display: "flex", flexDirection: "column", gap: 1 }}>
                <InlineField label="Use reflections within resolution limit of" hint="Angstrom">
                  <CCP4i2TaskElement itemName="ACOREF_RESOLL" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
                <InlineField label="to a limit of" hint="Angstrom">
                  <CCP4i2TaskElement itemName="ACOREF_RESOLU" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              </Box>
            )}
            <CCP4i2TaskElement
              itemName="ACORN_BEXCLUDE"
              {...props}
              qualifiers={{ guiLabel: "Exclude low SIGFP reflections" }}
              onChange={bexclude.onChange}
            />
            {bexclude.value && (
              <InlineField label="Reject reflections that are less than" hint="* the standard deviation (sigma)" sx={{ pl: 3 }}>
                <CCP4i2TaskElement itemName="ACOREF_EXCLUDE" {...props} qualifiers={{ guiLabel: " " }} />
              </InlineField>
            )}
            <CCP4i2TaskElement
              itemName="ACORN_BECUT"
              {...props}
              qualifiers={{ guiLabel: "Exclude reflections that are high E-value outliers" }}
              onChange={becut.onChange}
            />
            {becut.value && (
              <InlineField label="Reject observed reflections with E-values greater than" sx={{ pl: 3 }}>
                <CCP4i2TaskElement itemName="ACOREF_ECUT" {...props} qualifiers={{ guiLabel: " " }} />
              </InlineField>
            )}
          </CCP4i2ContainerElement>

          {/* Advanced Settings */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Advanced Settings" }}
            containerHint="BlockLevel"
          >
            <InlineField label="Upper density limit for Dynamic Density Modification (DDM)">
              <CCP4i2TaskElement
                itemName="ACOPH_CUTDDM"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <CCP4i2TaskElement
              itemName="ACORN_BGRID"
              {...props}
              qualifiers={{ guiLabel: "Choose a user defined grid size" }}
              onChange={bgrid.onChange}
            />
            {bgrid.value && (
              <InlineField label="Grid Size (Angstroms)" sx={{ pl: 3 }}>
                <CCP4i2TaskElement itemName="ACOGEN_GRID" {...props} qualifiers={{ guiLabel: " " }} />
              </InlineField>
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
