import { Paper, Stack, Box, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab label="Input" key="input">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Select input data" }}
            containerHint="FolderLevel"
            initiallyOpen={true}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="F_SIGF"
              qualifiers={{ guiLabel: "Reflections" }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Plot generation parameters" }}
            containerHint="FolderLevel"
            initiallyOpen={true}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="YLIM"
              qualifiers={{ guiLabel: "Range along y axis" }}
            />

            <Stack direction="row" alignItems="center" spacing={1} sx={{ mt: 1 }}>
              <Box sx={{ flexShrink: 0 }}>
                <CCP4i2TaskElement
                  {...props}
                  itemName="DLIM"
                  qualifiers={{ guiLabel: "High resolution cut-off" }}
                />
              </Box>
              <Typography variant="body1" color="text.secondary" sx={{ fontStyle: "italic" }}>
                No cut-off includes all data.
              </Typography>
            </Stack>

            <CCP4i2TaskElement
              {...props}
              itemName="SINGFIG"
              qualifiers={{ guiLabel: "Put all plots in one figure" }}
            />

            <Stack direction="row" alignItems="center" spacing={1}>
              <Box sx={{ flexShrink: 0 }}>
                <CCP4i2TaskElement
                  {...props}
                  itemName="FLAGICE"
                  qualifiers={{ guiLabel: "Flag suspected ice rings red" }}
                />
              </Box>
              <Typography variant="body1" color="text.secondary" sx={{ fontStyle: "italic" }}>
                Automatic ice ring detection.
              </Typography>
            </Stack>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
