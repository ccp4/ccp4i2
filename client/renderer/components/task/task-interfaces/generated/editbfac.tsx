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
            Atomic model
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="XYZIN" {...props} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Options
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Select B-factor treatment option - it is important this is set correctly
            </Typography>
            <CCP4i2TaskElement itemName="BTREATMENT" {...props} qualifiers={{ guiLabel: "" }} />
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Options for residues and regions
            </Typography>
            <CCP4i2TaskElement itemName="CONFCUT" {...props} qualifiers={{ guiLabel: "Remove low confidence residues (lddt or rmsd) (recommended)" }} />
            <CCP4i2TaskElement itemName="COMPACTREG" {...props} qualifiers={{ guiLabel: "Split model into compact regions (recommended)" }} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Optional Files
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              PAE file
            </Typography>
            <CCP4i2TaskElement itemName="PAEIN" {...props} />
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Distance model file
            </Typography>
            <CCP4i2TaskElement itemName="XYZDISTMOD" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
        <CCP4i2Tab key="addSettings" label="Additional Settings">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Options
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="MAXDOM" {...props} qualifiers={{ guiLabel: "Maximum domains to obtain." }} />
            <CCP4i2TaskElement itemName="DOMAINSIZE" {...props} qualifiers={{ guiLabel: "Approximate size of domains to be found (Angstroms)" }} />
            <CCP4i2TaskElement itemName="MINDOML" {...props} qualifiers={{ guiLabel: "Minimum domain length (residues)" }} />
            <CCP4i2TaskElement itemName="MAXFRACCL" {...props} qualifiers={{ guiLabel: "Maximum fraction close" }} />
            <CCP4i2TaskElement itemName="MINSEQRESI" {...props} qualifiers={{ guiLabel: "Minimum sequential residues" }} />
            <CCP4i2TaskElement itemName="MINREMSEQL" {...props} qualifiers={{ guiLabel: "Minimum remainder sequence length" }} />
            <CCP4i2TaskElement itemName="MINLDDT" {...props} qualifiers={{ guiLabel: "Minimum LDDT for removing residues" }} />
            <CCP4i2TaskElement itemName="MAXRMSD" {...props} qualifiers={{ guiLabel: "Maximum RMSD for removing residues" }} />
            <CCP4i2TaskElement itemName="SUBMINB" {...props} qualifiers={{ guiLabel: "Subtract the lowest B-value from all B-values" }} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
              PAE File Options (if PAE matrix supplied)
            </Typography>
            <CCP4i2TaskElement itemName="PAEPOWER" {...props} qualifiers={{ guiLabel: "PAE power" }} />
            <CCP4i2TaskElement itemName="PAECUTOFF" {...props} qualifiers={{ guiLabel: "PAE cutoff" }} />
            <CCP4i2TaskElement itemName="PAEGRAPHRES" {...props} qualifiers={{ guiLabel: "PAE graph resolution" }} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Distance Model Options (if suitable model provided)
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="WEIGHTCA" {...props} qualifiers={{ guiLabel: "Weight by CA-CA distance (if distance_model supplied)" }} />
            <CCP4i2TaskElement itemName="DISTPOW" {...props} qualifiers={{ guiLabel: "Distance power (for weighting by CA-CA distance)" }} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;