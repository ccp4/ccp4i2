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
import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { InlineField } from "../../task-elements/inline-field";
import { useJob } from "../../../../utils";
import { useBoolToggle } from "../../task-elements/shared-hooks";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const provideSequences = useBoolToggle(useTaskItem, "PROVIDESEQUENCES");
  const provideTls = useBoolToggle(useTaskItem, "PROVIDETLS");
  const provideDict = useBoolToggle(useTaskItem, "PROVIDEDICT");

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* Input data - single folder GUI */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Input data" }}
        containerHint="FolderLevel"
      >
        <InlineField label="Refinement in final step used reflection">
          <CCP4i2TaskElement itemName="USINGIORF" {...props} qualifiers={{ guiLabel: " " }} />
        </InlineField>

        <InlineField label="Provide sequences of crystallized species">
          <CCP4i2TaskElement itemName="PROVIDESEQUENCES" {...props} qualifiers={{ guiLabel: " " }} />
        </InlineField>
        {provideSequences.value && (
          <CCP4i2TaskElement itemName="ASUIN" {...props} qualifiers={{ guiLabel: "AU content - sequences" }} />
        )}

        <InlineField label="Provide Refined TLS parameters">
          <CCP4i2TaskElement itemName="PROVIDETLS" {...props} qualifiers={{ guiLabel: " " }} />
        </InlineField>
        {provideTls.value && (
          <CCP4i2TaskElement itemName="TLSIN" {...props} qualifiers={{ guiLabel: "Refined TLS parameter object" }} />
        )}

        <InlineField label="Provide DICT file for ligand in structure">
          <CCP4i2TaskElement itemName="PROVIDEDICT" {...props} qualifiers={{ guiLabel: " " }} />
        </InlineField>
        {provideDict.value && (
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ initiallyOpen: true }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement itemName="DICT" {...props} />
          </CCP4i2ContainerElement>
        )}

        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ guiLabel: "NEW: Allow ANISO use in zero cycle refinement", initiallyOpen: true }}
          containerHint="BlockLevel"
        >
          <InlineField label="Temperature factors">
            <CCP4i2TaskElement itemName="B_REFINEMENT_MODE" {...props} qualifiers={{ guiLabel: " " }} />
          </InlineField>
        </CCP4i2ContainerElement>

        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ initiallyOpen: true }}
          containerHint="BlockLevel"
        >
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            Choose a directory in which to put "Reflections.cif" and "Coordinates.cif"
            which can be uploaded to the deposition service
          </Typography>
          <CCP4i2TaskElement itemName="OUTPUT_DIRECTORY" {...props} />
        </CCP4i2ContainerElement>
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
