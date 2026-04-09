import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { InlineField } from "../task-elements/inline-field";
import { useBoolToggle } from "../task-elements/shared-hooks";
import { useJob } from "../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const twoDatasets = useBoolToggle(useTaskItem, "TWO_DATASETS");

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Model data" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="XYZIN_1" {...props} />
        <CCP4i2TaskElement itemName="F_SIGF_1" {...props} />
        <InlineField label="Label dataset 1 as">
          <CCP4i2TaskElement
            itemName="NAME_1"
            {...props}
            qualifiers={{ guiLabel: " " }}
          />
        </InlineField>
      </CCP4i2ContainerElement>

      <CCP4i2TaskElement
        itemName="TWO_DATASETS"
        {...props}
        qualifiers={{
          guiLabel:
            "Compare against a second dataset, e.g. before and after model building/refinement",
        }}
      />

      {twoDatasets.value && (
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ initiallyOpen: true }}
          containerHint="BlockLevel"
        >
          <CCP4i2TaskElement itemName="XYZIN_2" {...props} />
          <CCP4i2TaskElement itemName="F_SIGF_2" {...props} />
          <InlineField label="Label dataset 2 as">
            <CCP4i2TaskElement
              itemName="NAME_2"
              {...props}
              qualifiers={{ guiLabel: " " }}
            />
          </InlineField>
        </CCP4i2ContainerElement>
      )}

      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{
          guiLabel: "Which validation and analysis tools should be run?",
        }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="DO_IRIS"
          {...props}
          qualifiers={{
            guiLabel: "Iris: interactive multimetric validation charts",
          }}
        />
        <CCP4i2TaskElement
          itemName="DO_MOLPROBITY"
          {...props}
          qualifiers={{
            guiLabel:
              "MolProbity: checks including C-betas, side-chain flips, omega angles, and clashes",
          }}
        />
        <CCP4i2TaskElement
          itemName="DO_TORTOIZE"
          {...props}
          qualifiers={{
            guiLabel:
              "Tortoize: calculate conformation-dependent Z-scores for backbone geometry",
          }}
        />
        <CCP4i2TaskElement
          itemName="DO_BFACT"
          {...props}
          qualifiers={{
            guiLabel:
              "B-factors: average B-factor graphs and breakdown tables for publications",
          }}
        />
        <CCP4i2TaskElement
          itemName="DO_RAMA"
          {...props}
          qualifiers={{ guiLabel: "Ramachandran plots" }}
        />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
