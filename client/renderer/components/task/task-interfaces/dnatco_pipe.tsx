import { useCallback } from "react";
import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { InlineField } from "../task-elements/inline-field";
import { useJob } from "../../../utils";
import { isTruthy, useBoolToggle } from "../task-elements/shared-hooks";

/**
 * DNATCO — validate one nucleic acid model, or two side by side, and
 * optionally write NtC restraints for the next refinement.
 *
 * Inferred from the Qt task widget (dnatco_pipe_gui.py): one "Input data"
 * folder holding the model, a "compare with another model" toggle revealing
 * the second model, and a "generate restraints" toggle revealing the two
 * restraint parameters.
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const compare = useBoolToggle(useTaskItem, "TOGGLE_XYZIN2");
  const restraints = useBoolToggle(useTaskItem, "GENERATE_RESTRAINTS");
  const { value: XYZIN2, forceUpdate: forceSetXYZIN2 } = useTaskItem("XYZIN2");

  // Turning the comparison off also clears the second model, so a file the
  // user can no longer see is not silently validated by the server.
  const handleCompareToggle = useCallback(
    async (updatedItem: any) => {
      const enabled = isTruthy(updatedItem._value);
      await compare.onChange(updatedItem);
      if (!enabled && XYZIN2?.dbFileId) {
        forceSetXYZIN2({});
      }
    },
    [XYZIN2, forceSetXYZIN2, compare.onChange]
  );

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Input data" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="XYZIN1"
          {...props}
          qualifiers={{
            guiLabel: "Nucleic acid structure model",
            toolTip: "Structure model (PDB or mmCIF) containing DNA and/or RNA",
          }}
        />

        <InlineField label="Compare with another structure model" width="auto">
          <CCP4i2TaskElement
            itemName="TOGGLE_XYZIN2"
            {...props}
            qualifiers={{ guiLabel: " " }}
            onChange={handleCompareToggle}
          />
        </InlineField>
        {compare.value && (
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Second structure model", initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <Typography variant="body2" sx={{ fontStyle: "italic" }}>
              Typically the same structure after refinement or rebuilding. The report sets the
              two models side by side; restraints, when requested, are generated from this one.
            </Typography>
            <CCP4i2TaskElement
              itemName="XYZIN2"
              {...props}
              qualifiers={{ guiLabel: "Second structure model" }}
            />
          </CCP4i2ContainerElement>
        )}

        <InlineField
          label="Generate NtC restraints for refinement"
          hint="external restraints for Refmac/Servalcat"
          width="auto"
        >
          <CCP4i2TaskElement
            itemName="GENERATE_RESTRAINTS"
            {...props}
            qualifiers={{ guiLabel: " " }}
            onChange={restraints.onChange}
          />
        </InlineField>
        {restraints.value && (
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Parameters for restraints generation", initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <InlineField label="Maximum allowed NtC RMSD" hint="Å [default 0.5]">
              <CCP4i2TaskElement itemName="MAX_RMSD" {...props} qualifiers={{ guiLabel: " " }} />
            </InlineField>
            <InlineField label="Restraints sigma factor" hint="[default 1.0]">
              <CCP4i2TaskElement itemName="RESTRAINTS_SIGMA" {...props} qualifiers={{ guiLabel: " " }} />
            </InlineField>
          </CCP4i2ContainerElement>
        )}
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
