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
  const symmetrySource = useBoolToggle(useTaskItem, "SYMMETRY_SOURCE");

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
        <CCP4i2TaskElement itemName="HKLIN" {...props} />
        <CCP4i2TaskElement itemName="HKLIN1" {...props} />
        <CCP4i2TaskElement itemName="HKLIN2" {...props} />
        <CCP4i2TaskElement itemName="N_BINS" {...props} qualifiers={{ guiLabel: "Number of resolution bins" }} />
        <CCP4i2TaskElement itemName="SYMMETRY_SOURCE" {...props} qualifiers={{ guiLabel: "Load symmetry from" }} />
        {symmetrySource.value && (
          <CCP4i2TaskElement itemName="REFERENCEFILE" {...props} />
        )}
        {symmetrySource.value && (
          <CCP4i2TaskElement itemName="CELLFILE" {...props} />
        )}
        {symmetrySource.value && (
          <CCP4i2TaskElement itemName="STREAMFILE" {...props} />
        )}
        <CCP4i2TaskElement itemName="SPACEGROUP" {...props} qualifiers={{ guiLabel: "Space group" }} />
        <CCP4i2TaskElement itemName="CELL" {...props} qualifiers={{ guiLabel: "Unit cell" }} />
        <CCP4i2TaskElement itemName="WAVELENGTH" {...props} qualifiers={{ guiLabel: "Wavelength (A)" }} />
        <CCP4i2TaskElement itemName="D_MAX" {...props} qualifiers={{ guiLabel: "Low resolution cutoff" }} />
        <CCP4i2TaskElement itemName="D_MIN" {...props} qualifiers={{ guiLabel: "High resolution cutoff" }} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
