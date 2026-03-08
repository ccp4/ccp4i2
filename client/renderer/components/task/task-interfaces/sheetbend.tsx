import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { InlineField } from "../task-elements/inline-field";
import { useJob } from "../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { container } = useJob(props.job.id);

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input Data">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Input reflections" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="F_SIGF"
              {...props}
              qualifiers={{ toolTip: "Input reflections" }}
            />
            <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
            <CCP4i2TaskElement itemName="XYZIN" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab key="controlParameters" label="Basic Options">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Parameters to refine" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="REFINE_COORD"
              {...props}
              qualifiers={{ guiLabel: "coordinates" }}
            />
            <CCP4i2TaskElement
              itemName="REFINE_U_ISO"
              {...props}
              qualifiers={{ guiLabel: "isotropic B factors" }}
            />
            <CCP4i2TaskElement
              itemName="REFINE_U_ANISO"
              {...props}
              qualifiers={{ guiLabel: "anisotropic B factors" }}
            />
          </CCP4i2ContainerElement>

          <InlineField label="Number of cycles">
            <CCP4i2TaskElement
              itemName="CYCLES"
              {...props}
              qualifiers={{ guiLabel: " " }}
            />
          </InlineField>
          <InlineField
            label="Resolution"
            hint="Angstroms"
          >
            <CCP4i2TaskElement
              itemName="RESOLUTION"
              {...props}
              qualifiers={{
                guiLabel: " ",
                toolTip: "Enter a single number or a comma separated list",
              }}
            />
          </InlineField>
        </CCP4i2Tab>

        <CCP4i2Tab key="advancedOptions" label="Advanced Options">
          <InlineField label="Sphere radius" hint="times resolution">
            <CCP4i2TaskElement
              itemName="RADIUS_SCALE"
              {...props}
              qualifiers={{ guiLabel: " " }}
            />
          </InlineField>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Refine-regularize macro cycles" }}
            containerHint="FolderLevel"
          >
            <InlineField label="Perform pseudo-regularization">
              <CCP4i2TaskElement
                itemName="PSEUDO_REGULARIZE"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="Number of refine-regularize macro-cycles">
              <CCP4i2TaskElement
                itemName="REFINE_REGULARIZE_CYCLES"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <Typography
              variant="subtitle2"
              sx={{ fontWeight: "bold", mt: 1, mb: 0.5 }}
            >
              Additional parameters to refine before regularization
            </Typography>
            <CCP4i2TaskElement
              itemName="POSTREFINE_COORD"
              {...props}
              qualifiers={{ guiLabel: "coordinates" }}
            />
            <CCP4i2TaskElement
              itemName="POSTREFINE_U_ISO"
              {...props}
              qualifiers={{ guiLabel: "isotropic B factors" }}
            />
            <CCP4i2TaskElement
              itemName="POSTREFINE_U_ANISO"
              {...props}
              qualifiers={{ guiLabel: "anisotropic B factors" }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
