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
import { useJob } from "../../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* Single folder: Xia2 runs */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Xia2 runs" }}
        containerHint="FolderLevel"
      >
        <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
          Xia2 directory to import
        </Typography>
        <CCP4i2TaskElement itemName="XIA2_DIRECTORY" {...props} qualifiers={{ toolTip: "Browse to the directory \"xia2\"" }} />
        <Typography variant="body2" color="text.secondary" sx={{ mt: 1, mb: 1 }}>
          When a valid top-level XIA2 directory is selected above, the list of successful
          data reduction protocols that XIA2 performed will be summarised below.
          Select and delete ("-") the protocols you do not wish to import.
        </Typography>
        <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
          Protocols to keep from this XIA2 directory
        </Typography>
        <CCP4i2TaskElement itemName="runSummaries" {...props} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
