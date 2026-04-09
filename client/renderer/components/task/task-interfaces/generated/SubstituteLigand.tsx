import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: LIGANDAS } = useTaskItem("LIGANDAS");
  const { value: OBSAS } = useTaskItem("OBSAS");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input Data">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Ligand geometry
          </Typography>
          <CCP4i2TaskElement itemName="LIGANDAS" {...props} qualifiers={{ guiLabel: "Format in which geometry will be specified:" }} />
          {(LIGANDAS === "MOL") && (
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement itemName="MOLIN" {...props} />
            </CCP4i2ContainerElement>
          )}
          {(LIGANDAS === "DICT") && (
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement itemName="DICTIN" {...props} />
            </CCP4i2ContainerElement>
          )}
          {(LIGANDAS === "SMILES") && (
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement itemName="SMILESIN" {...props} />
            </CCP4i2ContainerElement>
          )}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="PIPELINE" {...props} qualifiers={{ guiLabel: "For rigid body refinement use" }} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Starting PDB In the appropriate spacegroup.
          </Typography>
          <CCP4i2TaskElement itemName="XYZIN" {...props} />
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Reflection data
          </Typography>
          <CCP4i2TaskElement itemName="OBSAS" {...props} qualifiers={{ guiLabel: "Format in which reflections will be specified:" }} />
          {(OBSAS === "UNMERGED") && (
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement itemName="UNMERGEDFILES" {...props} />
            </CCP4i2ContainerElement>
          )}
          {(OBSAS === "MERGED") && (
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement itemName="F_SIGF_IN" {...props} />
            </CCP4i2ContainerElement>
          )}
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Free-R flags
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="FREERFLAG_IN" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;