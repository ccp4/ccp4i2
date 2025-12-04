import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { useJob } from "../../../../utils";
import { useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: REFERENCE } = useTaskItem("REFERENCE");
  
  if (!container) return <LinearProgress />;
  
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input">
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Reflection objects to manipulate
          </Typography>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="F_SIGF" {...props} />
            <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
          </CCP4i2ContainerElement>
          <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
            Spacegroup and indexing
          </Typography>
          {(REFERENCE === "ANALYSE" || REFERENCE === "EXPAND") && (
            <CCP4i2TaskElement itemName="REFERENCE" {...props} qualifiers={{ guiLabel: "Define new indexing and spacegroup using" }} />
          )}
          {(REFERENCE === "ANALYSE") && (
            <CCP4i2TaskElement itemName="REFERENCE" {...props} qualifiers={{ guiLabel: "Analyse data symmetry" }} />
          )}
          {(REFERENCE === "EXPAND") && (
            <CCP4i2TaskElement itemName="REFERENCE" {...props} qualifiers={{ guiLabel: "Expand to space group P1" }} />
          )}
          {(REFERENCE === "HKLIN_FOBS_REF") && (
            <CCP4i2TaskElement itemName="HKLIN_FOBS_REF" {...props} />
          )}
          {(REFERENCE === "HKLIN_FC_REF") && (
            <CCP4i2TaskElement itemName="HKLIN_FC_REF" {...props} />
          )}
          {(REFERENCE === "HKLIN_FMAP_REF") && (
            <CCP4i2TaskElement itemName="HKLIN_FMAP_REF" {...props} />
          )}
          {(REFERENCE === "XYZIN_REF") && (
            <CCP4i2TaskElement itemName="XYZIN_REF" {...props} />
          )}
          {(REFERENCE === "SPECIFY") && (
            <>
              <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
                Note that reindexing a merged file is only valid within the same point group
              </Typography>
              <CCP4i2TaskElement itemName="CHOOSE_SPACEGROUP" {...props} />
              <CCP4i2TaskElement itemName="REINDEX_OPERATOR" {...props} />
              <CCP4i2TaskElement itemName="USE_REINDEX" {...props} qualifiers={{ guiLabel: "use reindex operator" }} />
              <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                You can give either a spacegroup, or a reindex operator (eg k,l,h), or both
              </Typography>
              <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                If only one of these is given, Pointless is usually able to generate the other automatically, but you should check
              </Typography>
              <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                The validity and consistency will be checked by Pointless: note warnings
              </Typography>
            </>
          )}
          {(REFERENCE === "LATTICE") && (
            <>
              <Typography variant="subtitle1" sx={{ fontWeight: "bold", mt: 2, mb: 1 }}>
                Remove centred lattice absences: BEWARE dangerous
              </Typography>
              <CCP4i2TaskElement itemName="LATTICE_CENTERING" {...props} qualifiers={{ guiLabel: "Desired lattice centering type" }} />
              <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                This option should ONLY be used if you are sure that the wrong cell was used in integration
              </Typography>
              <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                Note that not all centred lattices are consisent with all Bravais lattices, check the result carefully
              </Typography>
              <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                Lattice type "P" is ignored
              </Typography>
            </>
          )}
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;