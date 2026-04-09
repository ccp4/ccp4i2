import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: MOLSMILESORSKETCH } = useTaskItem("MOLSMILESORSKETCH");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            
          </Typography>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            
          </Typography>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Start point
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="MOLSMILESORSKETCH" {...props} qualifiers={{ guiLabel: "Start with molecular structure from" }} />
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Will launch Lidia to sketch molecule. Click Apply and Close in Lidia when sketch is ready.
            </Typography>
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Optionally can provide a starting monomer for the Lidia sketch:
            </Typography>
            {(MOLSMILESORSKETCH === "SKETCH" || MOLSMILESORSKETCH === "MOL") && (
              <CCP4i2TaskElement itemName="MOLIN" {...props} />
            )}
            {(MOLSMILESORSKETCH === "MOL2") && (
              <CCP4i2TaskElement itemName="MOL2IN" {...props} />
            )}
            {(MOLSMILESORSKETCH === "SMILES") && (
              <CCP4i2TaskElement itemName="SMILESIN" {...props} />
            )}
            {(MOLSMILESORSKETCH === "SMILESFILE") && (
              <CCP4i2TaskElement itemName="SMILESFILEIN" {...props} qualifiers={{ guiLabel: "SMILES file" }} />
            )}
            {(MOLSMILESORSKETCH === "DICT") && (
              <CCP4i2TaskElement itemName="DICTIN2" {...props} />
            )}
            {(MOLSMILESORSKETCH === "PDBMMCIF") && (
              <CCP4i2TaskElement itemName="PDBMMCIFIN" {...props} />
            )}
            <CCP4i2TaskElement itemName="TOGGLE_METAL" {...props} qualifiers={{ guiLabel: "This monomer contains a metal atom." }} />
            <CCP4i2TaskElement itemName="METAL_STRUCTURE" {...props} qualifiers={{ guiLabel: "indent" }} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Output monomer
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="TLC" {...props} qualifiers={{ guiLabel: "Monomer code (usually 3 or 5 characters)" }} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Advanced
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            {(MOLSMILESORSKETCH === "MOL" || MOLSMILESORSKETCH === "MOL2" || MOLSMILESORSKETCH === "DICT" || MOLSMILESORSKETCH === "PDBMMCIF") && (
              <CCP4i2TaskElement itemName="USE_COORD" {...props} qualifiers={{ guiLabel: "Use the coordinates from the input file for further optimisation (requires a reasonable input conformation)" }} />
            )}
            <CCP4i2TaskElement itemName="NRANDOM" {...props} qualifiers={{ guiLabel: "Set number of initial conformers to try:" }} />
            <CCP4i2TaskElement itemName="NOPROT" {...props} qualifiers={{ guiLabel: "No further protonation/deprotonation to be done by AceDRG" }} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;