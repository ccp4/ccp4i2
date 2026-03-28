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
  const { value: database } = useTaskItem("DATABASE");
  const useApi = useBoolToggle(useTaskItem, "USEAPI");

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
        <CCP4i2TaskElement itemName="SEQIN" {...props} />
        <CCP4i2TaskElement itemName="F_SIGF" {...props} />
        <CCP4i2TaskElement itemName="MAXHITS" {...props} />
        <CCP4i2TaskElement itemName="DATABASE" {...props} />
        {/* PDB settings: shown when database is NOT 'AFDB' (i.e. PDB or All) */}
        {database !== "AFDB" && (
          <>
            <CCP4i2TaskElement itemName="PDBLOCAL" {...props} />
            <CCP4i2TaskElement itemName="PDBSEQDB" {...props} />
          </>
        )}
        {/* AFDB settings: shown when database is NOT 'PDB' (i.e. AFDB or All) */}
        {database !== "PDB" && (
          <>
            <CCP4i2TaskElement itemName="USEAPI" {...props} />
            {!useApi.value && (
              <CCP4i2TaskElement itemName="AFDBSEQDB" {...props} />
            )}
          </>
        )}
        <CCP4i2TaskElement itemName="NPROC" {...props} />
        <CCP4i2TaskElement itemName="DO_CLASSIFY" {...props} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
