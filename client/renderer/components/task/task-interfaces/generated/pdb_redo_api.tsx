import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { useJob } from "../../../../utils";

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
            qualifiers={{ guiLabel: "Main Input data" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="XYZIN"
              {...props}
              qualifiers={{
                toolTip: "This is the structure which will be refined.",
              }}
            />
            <CCP4i2TaskElement itemName="F_SIGF" {...props} />
            <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Optional input data" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="SEQIN" {...props} />
            <CCP4i2TaskElement itemName="DICT" {...props} />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Paired refinement" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="PAIRED"
              {...props}
              qualifiers={{ guiLabel: "Perform paired refinement" }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Rebuilding Options" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="NOLOOPS"
              {...props}
              qualifiers={{ guiLabel: "Do not try to complete loops" }}
            />
            <CCP4i2TaskElement
              itemName="NOPEPFLIP"
              {...props}
              qualifiers={{ guiLabel: "Do not perform peptide flips" }}
            />
            <CCP4i2TaskElement
              itemName="NOSCBUILD"
              {...props}
              qualifiers={{
                guiLabel: "Do not rebuild or complete side-chain",
              }}
            />
            <CCP4i2TaskElement
              itemName="NOCENTRIFUGE"
              {...props}
              qualifiers={{ guiLabel: "Do not delete poor waters" }}
            />
            <CCP4i2TaskElement
              itemName="NOSUGARBUILD"
              {...props}
              qualifiers={{ guiLabel: "Do not (re)build carbohydrates" }}
            />
            <CCP4i2TaskElement
              itemName="NOREBUILD"
              {...props}
              qualifiers={{ guiLabel: "Skip model rebuilding" }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab key="controlParameters" label="Advanced Options">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Advanced Options" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="NEWMODEL"
              {...props}
              qualifiers={{
                guiLabel:
                  "Update the model even if the initial refinement is poor",
              }}
            />
            <CCP4i2TaskElement
              itemName="ISOTROPIC"
              {...props}
              qualifiers={{ guiLabel: "Force isotropic B-factors" }}
            />
            <CCP4i2TaskElement
              itemName="ANISOTROPIC"
              {...props}
              qualifiers={{
                guiLabel: "Force anisotropic B-factors (within limits)",
              }}
            />
            <CCP4i2TaskElement
              itemName="NOTLS"
              {...props}
              qualifiers={{ guiLabel: "Do not perform TLS refinement" }}
            />
            <CCP4i2TaskElement
              itemName="TIGHTER"
              {...props}
              qualifiers={{
                guiLabel: "Use tighter restraints in refinement",
              }}
            />
            <CCP4i2TaskElement
              itemName="LOOSER"
              {...props}
              qualifiers={{
                guiLabel: "Use looser restraints in refinement",
              }}
            />
            <Typography variant="body2" color="text.secondary">
              Note: The tighter and looser values are subtracted so the net
              effect of having both set to 2 will be none
            </Typography>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
