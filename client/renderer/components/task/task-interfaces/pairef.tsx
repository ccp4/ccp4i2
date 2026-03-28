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
  const { value: shType } = useTaskItem("SH_TYPE");
  const usePreref = useBoolToggle(useTaskItem, "USE_PREREF");
  const useShake = useBoolToggle(useTaskItem, "USE_SHAKE");
  const autoWgt = useBoolToggle(useTaskItem, "AUTO_WGT");
  const fixedTls = useBoolToggle(useTaskItem, "FIXED_TLS");

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
        <CCP4i2TaskElement itemName="SH_TYPE" {...props} qualifiers={{ guiLabel: "Run Pairef with" }} />
        <CCP4i2TaskElement itemName="USE_PREREF" {...props} />
        <CCP4i2TaskElement itemName="F_SIGF" {...props} />
        <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
        <CCP4i2TaskElement itemName="XYZIN" {...props} />
        <CCP4i2TaskElement itemName="UNMERGED" {...props} />
        <CCP4i2TaskElement itemName="DICT" {...props} />
        <CCP4i2TaskElement itemName="REFMAC_KEYWORD_FILE" {...props} />
        {shType === "semi" && (
          <>
            <CCP4i2TaskElement itemName="NSHELL" {...props} qualifiers={{ guiLabel: "Add" }} />
            <CCP4i2TaskElement itemName="WSHELL" {...props} qualifiers={{ guiLabel: "resolution shells of width" }} />
          </>
        )}
        {shType === "manual" && (
          <CCP4i2TaskElement itemName="MANSHELL" {...props} qualifiers={{ guiLabel: "Explicitly define shells" }} />
        )}
        <CCP4i2TaskElement itemName="NPRECYCLES" {...props} qualifiers={{ guiLabel: "Number of pre-refinement cycles is" }} />
        {useShake.value && (
          <CCP4i2TaskElement itemName="SHAKE" {...props} />
        )}
        <CCP4i2TaskElement itemName="USE_SHAKE" {...props} />
        <CCP4i2TaskElement itemName="RESETBFAC" {...props} />
        <CCP4i2TaskElement itemName="COMPLETE" {...props} />
        <CCP4i2TaskElement itemName="INIRES" {...props} qualifiers={{ guiLabel: "Manually set the initial resolution (use if required)" }} />
        <CCP4i2TaskElement itemName="NCYCLES" {...props} qualifiers={{ guiLabel: "Number of refinement cycles to perform" }} />
        <CCP4i2TaskElement itemName="AUTO_WGT" {...props} />
        {!autoWgt.value && (
          <CCP4i2TaskElement itemName="WGT_TRM" {...props} />
        )}
        <CCP4i2TaskElement itemName="FIXED_TLS" {...props} />
        {fixedTls.value && (
          <CCP4i2TaskElement itemName="TLSCYC" {...props} />
        )}
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
