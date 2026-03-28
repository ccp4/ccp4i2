/*
 * Copyright (C) 2025-2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
import { Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";
import { useCallback } from "react";
import { useBoolToggle, isTruthy } from "../task-elements/shared-hooks";
import { InlineField } from "../task-elements/inline-field";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);
  const { forceUpdate: forceSetCYCLES } = useTaskItem("CYCLES");

  const useModelPhases = useBoolToggle(useTaskItem, "USE_MODEL_PHASES");

  const handleBASIC = useCallback(
    async (updatedItem: any) => {
      await forceSetCYCLES(isTruthy(updatedItem._value) ? 5 : 25);
    },
    [forceSetCYCLES]
  );

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* Reflection data */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Reflection data" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="F_SIGF" {...props} />
        <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
        <CCP4i2TaskElement
          itemName="USE_MODEL_PHASES"
          {...props}
          qualifiers={{
            guiLabel:
              "Calculate initial phases by refining the starting model (uncheck to specify starting phases, e.g. from experimental phasing)",
          }}
        />
        {!useModelPhases.value && (
          <>
            <CCP4i2TaskElement itemName="PHASES" {...props} />
            <CCP4i2TaskElement
              itemName="UNBIASED"
              {...props}
              qualifiers={{
                guiLabel:
                  "Phases are unbiased and should be used as refinement restraints when the model is poor",
              }}
            />
          </>
        )}
      </CCP4i2ContainerElement>

      {/* Asymmetric unit contents */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Asymmetric unit contents" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="ASUIN" {...props} />
      </CCP4i2ContainerElement>

      {/* Starting model */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Starting model (optional)" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="XYZIN" {...props} />
        <Typography variant="caption" color="text.secondary" sx={{ pl: 1 }}>
          If a starting model is provided it will always be used, regardless of
          the phase source. Clear this field for a fully de novo build.
        </Typography>
      </CCP4i2ContainerElement>

      {/* Options */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Options" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="BASIC"
          {...props}
          qualifiers={{ guiLabel: "Run a quicker basic pipeline" }}
          onChange={handleBASIC}
        />

        <InlineField label="Run for" hint="cycles">
          <CCP4i2TaskElement
            itemName="CYCLES"
            {...props}
            qualifiers={{ guiLabel: " " }}
          />
        </InlineField>

        <InlineField
          width="auto"
          after={
            <InlineField
              label="Stop automatically if R-free does not improve in"
              hint="cycles"
            >
              <CCP4i2TaskElement
                itemName="STOP_CYCLES"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
          }
        >
          <CCP4i2TaskElement
            itemName="AUTO_STOP"
            {...props}
            qualifiers={{ guiLabel: " " }}
            sx={{ width: "auto" }}
          />
        </InlineField>

        <CCP4i2TaskElement
          itemName="SELENOMET"
          {...props}
          qualifiers={{
            guiLabel:
              "Build selenomethionine (MSE) instead of methionine (MET)",
          }}
        />
        <CCP4i2TaskElement
          itemName="TWINNED"
          {...props}
          qualifiers={{ guiLabel: "Use twinned refinement" }}
        />
      </CCP4i2ContainerElement>

      {/* Optional pipeline steps */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Optional pipeline steps" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="SHEETBEND"
          {...props}
          qualifiers={{
            guiLabel:
              "Preliminary low-resolution refinement with Sheetbend",
          }}
        />
        <CCP4i2TaskElement
          itemName="PRUNING"
          {...props}
          qualifiers={{ guiLabel: "Residue and chain pruning" }}
        />
        <CCP4i2TaskElement
          itemName="PARROT"
          {...props}
          qualifiers={{
            guiLabel: "Classical density modification with Parrot",
          }}
        />
        <CCP4i2TaskElement
          itemName="DUMMY_ATOMS"
          {...props}
          qualifiers={{
            guiLabel:
              "Phase improvement through addition and refinement of dummy atoms",
          }}
        />
        <CCP4i2TaskElement
          itemName="WATERS"
          {...props}
          qualifiers={{ guiLabel: "Addition of waters" }}
        />
        <CCP4i2TaskElement
          itemName="SIDE_CHAIN_FIXING"
          {...props}
          qualifiers={{ guiLabel: "Final side-chain fixing" }}
        />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
