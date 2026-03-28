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
import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: watMode } = useTaskItem("WAT");

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* Input Data */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Input Data" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="F_SIGF" {...props} />
        <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
        <CCP4i2TaskElement itemName="XYZIN" {...props} />
        <CCP4i2TaskElement itemName="DICT" {...props} />
        <CCP4i2TaskElement itemName="NBCYCLES" {...props} qualifiers={{ guiLabel: "Perform" }} />
        <CCP4i2TaskElement itemName="NSCYCLES" {...props} qualifiers={{ guiLabel: "refinements (\"big cycles\") with a maximum of" }} />
        <CCP4i2TaskElement itemName="WAT" {...props} />
        {watMode === "MAN" && (
          <CCP4i2TaskElement itemName="WATCYC" {...props} qualifiers={{ guiLabel: "From big cycle" }} />
        )}
        <CCP4i2TaskElement itemName="AUTO_NCS" {...props} />
        <CCP4i2TaskElement itemName="RBR" {...props} />
        <CCP4i2TaskElement itemName="TLS" {...props} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
