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

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: PASTEORREAD } = useTaskItem("PASTEORREAD");

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* Single folder: inputData */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Optional objects from which to start definition" }}
        containerHint="FolderLevel"
      >
        <InlineField label="Paste or read alignment, or extract from HHPred or Blast search">
          <CCP4i2TaskElement itemName="PASTEORREAD" {...props} qualifiers={{ guiLabel: " " }} />
        </InlineField>

        {PASTEORREAD === "ALIGNIN" && (
          <CCP4i2TaskElement itemName="ALIGNIN" {...props} qualifiers={{ toolTip: "Alignment object or file" }} />
        )}
        {PASTEORREAD === "HHPREDIN" && (
          <CCP4i2TaskElement itemName="HHPREDIN" {...props} qualifiers={{ toolTip: "HHPred results" }} />
        )}
        {PASTEORREAD === "BLASTIN" && (
          <CCP4i2TaskElement itemName="BLASTIN" {...props} qualifiers={{ toolTip: "Blast results" }} />
        )}
        {PASTEORREAD === "PASTE" && (
          <CCP4i2TaskElement itemName="SEQUENCETEXT" {...props} />
        )}

        {(PASTEORREAD === "HHPREDIN" || PASTEORREAD === "BLASTIN") && (
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            Choose one alignment from the search file
          </Typography>
        )}

        <Typography variant="body2" color="text.secondary" sx={{ mt: 1, mb: 1 }}>
          Annotation for the alignment
        </Typography>
        <CCP4i2TaskElement itemName="ANNOTATION" {...props} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
