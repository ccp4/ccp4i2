import { useMemo } from "react";
import { Alert, Paper, Stack, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useJob } from "../../../utils";
import { useBoolToggle } from "../task-elements/shared-hooks";
import { InlineField } from "../task-elements/inline-field";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);

  // Watch GEN_MODE to conditionally show FREERFLAG input
  const { value: GEN_MODE_RAW } = useTaskItem("GEN_MODE");
  const isComplete = useMemo(() => GEN_MODE_RAW === "COMPLETE", [GEN_MODE_RAW]);

  // Watch CUTRESOLUTION to conditionally show RESMAX
  const cutResolution = useBoolToggle(useTaskItem, "CUTRESOLUTION");

  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab label="Input" key="input">
          <Alert severity="info" sx={{ mb: 2 }}>
            By default, this task will create a new set of freeR flags for
            cross-validation during refinement. Potential twinning operations
            will be taken into account.
          </Alert>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Reflection data" }}
            containerHint="FolderLevel"
            initiallyOpen={true}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="F_SIGF"
              qualifiers={{ guiLabel: "Reflections" }}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="GEN_MODE"
              qualifiers={{ guiLabel: "Mode", guiMode: "radio" }}
            />

            {isComplete && (
              <CCP4i2TaskElement
                {...props}
                itemName="FREERFLAG"
                qualifiers={{ guiLabel: "Existing freeR flags to complete" }}
              />
            )}
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Options" }}
            containerHint="FolderLevel"
            initiallyOpen={true}
          >
            <InlineField
              label="Fraction of reflections in freeR set"
              width="120px"
              hint="(default 0.05)"
            >
              <CCP4i2TaskElement
                {...props}
                itemName="FRAC"
                qualifiers={{ guiLabel: "" }}
              />
            </InlineField>

            <Stack direction="row" alignItems="center" spacing={1} sx={{ mt: 1 }}>
              <CCP4i2TaskElement
                {...props}
                itemName="CUTRESOLUTION"
                qualifiers={{ guiLabel: "Set high resolution limit at" }}
                onChange={cutResolution.onChange}
              />
              {cutResolution.value && (
                <InlineField hint={"\u00C5"} width="100px">
                  <CCP4i2TaskElement
                    {...props}
                    itemName="RESMAX"
                    qualifiers={{ guiLabel: "" }}
                  />
                </InlineField>
              )}
            </Stack>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
