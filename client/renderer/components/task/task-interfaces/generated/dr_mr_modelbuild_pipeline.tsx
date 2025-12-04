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
  const { value: MERGED_OR_UNMERGED } = useTaskItem("MERGED_OR_UNMERGED");
  const { value: SEARCH_AFDB } = useTaskItem("SEARCH_AFDB");
  const { value: SEARCH_PDB } = useTaskItem("SEARCH_PDB");
  const { value: XYZINORMRBUMP } = useTaskItem("XYZINORMRBUMP");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input Data and Protocol">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Input data
          </Typography>
          <CCP4i2TaskElement itemName="MERGED_OR_UNMERGED" {...props} qualifiers={{ guiLabel: "Input data type (run import merged, use merged from 'import merged' job or unmerged:" }} />
          {(MERGED_OR_UNMERGED === "UNMERGED") && (
            <CCP4i2TaskElement itemName="UNMERGEDFILES" {...props} />
          )}
          {(MERGED_OR_UNMERGED === "MERGED") && (
            <CCP4i2TaskElement itemName="F_SIGF_IN" {...props} />
          )}
          {(MERGED_OR_UNMERGED === "MERGED") && (
            <CCP4i2TaskElement itemName="FREER_IN" {...props} />
          )}
          {(MERGED_OR_UNMERGED === "MERGED_F") && (
            <CCP4i2TaskElement itemName="HKLIN" {...props} />
          )}
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Search model
          </Typography>
          <CCP4i2TaskElement itemName="XYZINORMRBUMP" {...props} />
          {(XYZINORMRBUMP === "XYZINPUT") && (
            <CCP4i2TaskElement itemName="XYZIN" {...props} />
          )}
          {(XYZINORMRBUMP === "MRBUMP") && (
            <>
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
              {(XYZINORMRBUMP === "MRBUMP") && (
                <CCP4i2TaskElement itemName="MRMAX" {...props} />
              )}
            </>
          )}
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Sequence of target model
          </Typography>
          <CCP4i2TaskElement itemName="ASUIN" {...props} />
          <CCP4i2TaskElement itemName="NMON" {...props} qualifiers={{ guiLabel: "The number of monomers to search for" }} />
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Options
          </Typography>
          <CCP4i2TaskElement itemName="AUTOCUTOFF" {...props} qualifiers={{ guiLabel: "Run Aimless twice, first to find resolution limit" }} />
          <CCP4i2TaskElement itemName="RUNACORN" {...props} qualifiers={{ guiLabel: "Run phase refinement with acorn before model building" }} />
          <CCP4i2TaskElement itemName="REFMAC_NCYC" {...props} qualifiers={{ guiLabel: "cycles of restrained refinement after MR" }} />
          <CCP4i2TaskElement itemName="BUCC_NCYC" {...props} qualifiers={{ guiLabel: "model building pipeline iterations" }} />
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
        </CCP4i2Tab>
        <CCP4i2Tab key="inputData" label="Advanced options">
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;