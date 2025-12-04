import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input Data">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            You do not have valid PDB-REDO token/secret set. Please set these using:  Utilities -&gt; System administration tools -&gt; Configure login tokens for PDB-REDO.  You can obtain a token/secret pair from https://services.pdb-redo.eu/
          </Typography>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Main Input data
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="XYZIN" {...props} qualifiers={{ toolTip: "This is the structure which will be refined." }} />
            <CCP4i2TaskElement itemName="F_SIGF" {...props} />
            <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Optional input data
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="SEQIN" {...props} />
            <CCP4i2TaskElement itemName="DICT" {...props} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Paired refinement
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="PAIRED" {...props} qualifiers={{ guiLabel: "Perform paired refinement" }} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Rebuilding Options
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="NOLOOPS" {...props} qualifiers={{ guiLabel: "Do not try to complete loops" }} />
            <CCP4i2TaskElement itemName="NOPEPFLIP" {...props} qualifiers={{ guiLabel: "Do not perform peptide flips" }} />
            <CCP4i2TaskElement itemName="NOSCBUILD" {...props} qualifiers={{ guiLabel: "Do not rebuild or complete side-chain" }} />
            <CCP4i2TaskElement itemName="NOCENTRIFUGE" {...props} qualifiers={{ guiLabel: "Do not delete poor waters" }} />
            <CCP4i2TaskElement itemName="NOSUGARBUILD" {...props} qualifiers={{ guiLabel: "Do not (re)build carbohydrates" }} />
            <CCP4i2TaskElement itemName="NOREBUILD" {...props} qualifiers={{ guiLabel: "Skip model rebuilding" }} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
        <CCP4i2Tab key="controlParameters" label="Advanced Options">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Advanced Options
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="NEWMODEL" {...props} qualifiers={{ guiLabel: "Update the model even if the initial refinement is poor" }} />
            <CCP4i2TaskElement itemName="ISOTROPIC" {...props} qualifiers={{ guiLabel: "Force isotropic B-factors" }} />
            <CCP4i2TaskElement itemName="ANISOTROPIC" {...props} qualifiers={{ guiLabel: "Force anisotropic B-factors (within limits)" }} />
            <CCP4i2TaskElement itemName="NOTLS" {...props} qualifiers={{ guiLabel: "Do not perform TLS refinement" }} />
            <CCP4i2TaskElement itemName="TIGHTER" {...props} qualifiers={{ guiLabel: "Use tighter restraints in refinement" }} />
            <CCP4i2TaskElement itemName="LOOSER" {...props} qualifiers={{ guiLabel: "Use looser restraints in refinement" }} />
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Note: The tighter and looser values are subtracted so the net effect of having both set to 2 will be none
            </Typography>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;