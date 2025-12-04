import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: SEARCH_AFDB } = useTaskItem("SEARCH_AFDB");
  const { value: SEARCH_PDB } = useTaskItem("SEARCH_PDB");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            Sequences from AU content:
          </Typography>
          <CCP4i2TaskElement itemName="ASUIN" {...props} />
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Model databases
          </Typography>
          <CCP4i2TaskElement itemName="SEARCH_PDB" {...props} qualifiers={{ guiLabel: "Search PDB for for possible MR search models" }} />
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            indentNon-redundancy level for homologue search:
          </Typography>
          {(SEARCH_PDB === true) && (
            <CCP4i2TaskElement itemName="REDUNDANCYLEVEL" {...props} qualifiers={{ guiLabel: "indent" }} />
          )}
          <CCP4i2TaskElement itemName="SEARCH_AFDB" {...props} qualifiers={{ guiLabel: "Search EBI-AFDB for possible MR search models" }} />
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            indentEBI-AFDB pLDDT residue score cut-off:
          </Typography>
          {(SEARCH_AFDB === true) && (
            <CCP4i2TaskElement itemName="AFDBLEVEL" {...props} qualifiers={{ guiLabel: "indent" }} />
          )}
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            Maximum no. of search models to create:
          </Typography>
          <CCP4i2TaskElement itemName="MRMAX" {...props} />
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Optional Settings
          </Typography>
          <CCP4i2TaskElement itemName="HHPREDIN" {...props} qualifiers={{ toolTip: "HHPred results" }} />
          <CCP4i2TaskElement itemName="PDBLOCAL" {...props} qualifiers={{ toolTip: "Local PDB mirror" }} />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;