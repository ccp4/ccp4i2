import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { InlineField } from "../../task-elements/inline-field";
import { useJob } from "../../../../utils";
import { useBoolToggle } from "../../task-elements/shared-hooks";

const COORD_NUMBERS = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 24] as const;

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: MAXIMUM_COORDINATION_NUMBER } = useTaskItem("MAXIMUM_COORDINATION_NUMBER");
  const savePdbmmcif = useBoolToggle(useTaskItem, "SAVE_PDBMMCIF");

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* Single folder: Input Data */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Input Data" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="XYZIN" {...props} qualifiers={{ toolTip: "Atomic model" }} />
        <InlineField label="Monomer code">
          <CCP4i2TaskElement itemName="LIGAND_CODE" {...props} qualifiers={{ guiLabel: " " }} />
        </InlineField>

        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ guiLabel: "Advanced parameters", initiallyOpen: true }}
          containerHint="BlockLevel"
        >
          <InlineField label="Distance threshold: (range 0-1)">
            <CCP4i2TaskElement itemName="DISTANCE_THRESHOLD" {...props} qualifiers={{ guiLabel: " " }} />
          </InlineField>
          <InlineField label="Maximum coordination number:">
            <CCP4i2TaskElement itemName="MAXIMUM_COORDINATION_NUMBER" {...props} qualifiers={{ guiLabel: " " }} />
          </InlineField>
          {COORD_NUMBERS.map((n) => {
            const key = `COORD${String(n).padStart(2, "0")}`;
            return (
              String(MAXIMUM_COORDINATION_NUMBER) === String(n) && (
                <InlineField key={key} label="Coordination class:">
                  <CCP4i2TaskElement itemName={key} {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )
            );
          })}
          <InlineField label="Procrustes distance threshold: (range 0-1)">
            <CCP4i2TaskElement itemName="PROCRUSTES_DISTANCE_THRESHOLD" {...props} qualifiers={{ guiLabel: " " }} />
          </InlineField>
          <InlineField label="Minimum sample size for statistics:">
            <CCP4i2TaskElement itemName="MINIMUM_SAMPLE_SIZE" {...props} qualifiers={{ guiLabel: " " }} />
          </InlineField>
          <CCP4i2TaskElement itemName="USE_PDB" {...props} qualifiers={{ guiLabel: "Use COD structures based on the input PDB/mmCIF coordinates" }} />
          <CCP4i2TaskElement itemName="IDEAL_ANGLES" {...props} qualifiers={{ guiLabel: "Provide only ideal bond angles" }} />
          <CCP4i2TaskElement itemName="SIMPLE" {...props} qualifiers={{ guiLabel: "Simple distance based filtering" }} />
          <CCP4i2TaskElement itemName="SAVE_PDBMMCIF" {...props} qualifiers={{ guiLabel: "Update link records to metal sites in the atomic model" }} />
          {savePdbmmcif.value && (
            <CCP4i2TaskElement itemName="KEEP_LINKS" {...props} qualifiers={{ guiLabel: "Keep existing link records to metal sites" }} />
          )}
        </CCP4i2ContainerElement>
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
