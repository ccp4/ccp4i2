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
import { useBoolToggle } from "../task-elements/shared-hooks";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const perform = useBoolToggle(useTaskItem, "PERFORM");

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
        <CCP4i2TaskElement itemName="PERFORM" {...props} />
        {!perform.value && (
          <CCP4i2TaskElement itemName="F_SIGF" {...props} />
        )}
        {perform.value && (
          <CCP4i2TaskElement itemName="XYZIN_FIX" {...props} />
        )}
        {!perform.value && (
          <CCP4i2TaskElement itemName="XYZIN" {...props} />
        )}
        {!perform.value && (
          <CCP4i2TaskElement itemName="ASUIN" {...props} />
        )}
        {!perform.value && (
          <CCP4i2TaskElement itemName="NMON" {...props} qualifiers={{ guiLabel: "The number of monomers to search for" }} />
        )}
        {perform.value && (
          <CCP4i2TaskElement itemName="F_PHI_MAP" {...props} />
        )}
        {perform.value && (
          <CCP4i2TaskElement itemName="XYZIN_FIX" {...props} />
        )}
      </CCP4i2ContainerElement>

      {/* Basic Options */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Basic Options" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="SEQ" {...props} />
        <CCP4i2TaskElement itemName="SURF" {...props} />
        <CCP4i2TaskElement itemName="NP" {...props} qualifiers={{ guiLabel: "Number of Rotation Function peaks" }} />
        <CCP4i2TaskElement itemName="NPT" {...props} qualifiers={{ guiLabel: "Number of Translation Function peaks" }} />
      </CCP4i2ContainerElement>

      {/* Advanced Options */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Advanced Options" }}
        containerHint="FolderLevel"
      >
        {perform.value && (
          <CCP4i2TaskElement itemName="PRF" {...props} />
        )}
        {perform.value && (
          <CCP4i2TaskElement itemName="SCORE" {...props} />
        )}
        {perform.value && (
          <CCP4i2TaskElement itemName="NMON_EXP" {...props} qualifiers={{ guiLabel: "Expected number of copies (for contrast calculation only)" }} />
        )}
        {perform.value && (
          <CCP4i2TaskElement itemName="ANISO" {...props} />
        )}
        <CCP4i2TaskElement itemName="HIGH_PATH_VAR" {...props} />
        <CCP4i2TaskElement itemName="LOW_PATH_VAR" {...props} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
