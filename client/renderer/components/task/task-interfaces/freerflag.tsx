import { useCallback, useEffect, useMemo, useState } from "react";
import { Alert, Box, Paper, Stack, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useJob } from "../../../utils";

const isTruthy = (val: any): boolean =>
  val === true || val === "True" || val === "true";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);

  // Watch GEN_MODE to conditionally show FREERFLAG input
  const { value: GEN_MODE_RAW } = useTaskItem("GEN_MODE");
  const isComplete = useMemo(() => GEN_MODE_RAW === "COMPLETE", [GEN_MODE_RAW]);

  // Watch CUTRESOLUTION to conditionally show RESMAX
  const { value: CUTRESOLUTION_RAW } = useTaskItem("CUTRESOLUTION");
  const [cutResolution, setCutResolution] = useState(() =>
    isTruthy(CUTRESOLUTION_RAW)
  );

  useEffect(() => {
    setCutResolution(isTruthy(CUTRESOLUTION_RAW));
  }, [CUTRESOLUTION_RAW]);

  const handleCutResolutionChange = useCallback((updatedItem: any) => {
    setCutResolution(isTruthy(updatedItem._value));
  }, []);

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
            <Stack direction="row" alignItems="center" spacing={1}>
              <Typography variant="body1" sx={{ whiteSpace: "nowrap" }}>
                Fraction of reflections in freeR set
              </Typography>
              <Box sx={{ width: 120 }}>
                <CCP4i2TaskElement
                  {...props}
                  itemName="FRAC"
                  qualifiers={{ guiLabel: "" }}
                />
              </Box>
              <Typography variant="body1" color="text.secondary">
                (default 0.05)
              </Typography>
            </Stack>

            <Stack direction="row" alignItems="center" spacing={1} sx={{ mt: 1 }}>
              <Box sx={{ flexShrink: 0 }}>
                <CCP4i2TaskElement
                  {...props}
                  itemName="CUTRESOLUTION"
                  qualifiers={{ guiLabel: "Set high resolution limit at" }}
                  onChange={handleCutResolutionChange}
                />
              </Box>
              {cutResolution && (
                <>
                  <Box sx={{ width: 100 }}>
                    <CCP4i2TaskElement
                      {...props}
                      itemName="RESMAX"
                      qualifiers={{ guiLabel: "" }}
                    />
                  </Box>
                  <Typography variant="body1">{"\u00C5"}</Typography>
                </>
              )}
            </Stack>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
