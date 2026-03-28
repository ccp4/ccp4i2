/*
 * Copyright (C) 2026 Newcastle University
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
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";
import { useBoolToggle } from "../task-elements/shared-hooks";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: mode } = useTaskItem("MODE");
  const needUserCellSymm = useBoolToggle(useTaskItem, "NEED_USER_CELLSYMM");
  const freerInHklin = useBoolToggle(useTaskItem, "FREER_IN_HKLIN");

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* Experimental data */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Experimental data" }}
        containerHint="FolderLevel"
      >
        <Typography variant="body2">Enter the file containing the data (supported formats: MTZ, mmCIF)</Typography>
        <CCP4i2TaskElement itemName="HKLIN" {...props} />
        {needUserCellSymm.value && (
          <CCP4i2TaskElement itemName="SPACEGROUPCELL" {...props} />
        )}
        {!freerInHklin.value && (
          <Typography variant="body2" color="warning.main">
            Your data file has no freeR flags, so we will generate them.
          </Typography>
        )}
        <CCP4i2TaskElement itemName="MODE" {...props} qualifiers={{ guiLabel: "How do you expect to solve the structure?" }} />
      </CCP4i2ContainerElement>

      {/* Composition for MR */}
      {mode === "molrep" && (
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ guiLabel: "Composition for MR" }}
          containerHint="FolderLevel"
        >
          <Typography variant="body2">Enter sequence of structure to be solved</Typography>
          <CCP4i2TaskElement itemName="SEQFILEIN" {...props} />
          <CCP4i2TaskElement itemName="CONTAINS_SEMET" {...props} qualifiers={{ guiLabel: "Structure contains selenomethionine" }} />
          <Typography variant="body2">Enter suggested homologous structures</Typography>
          <CCP4i2TaskElement itemName="MR_MODE" {...props} />
          <CCP4i2TaskElement itemName="MR_MODEL_DIR" {...props} />
          <CCP4i2TaskElement itemName="MR_MODEL_LIST" {...props} />
        </CCP4i2ContainerElement>
      )}
    </Paper>
  );
};

export default TaskInterface;
