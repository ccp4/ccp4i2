import React from "react";
import { Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs {...props}>
        {/* ===== Tab 1: Input data ===== */}
        <CCP4i2Tab label="Input data">
          {/* Input reflections */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Input reflections" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="F_SIGF"
              {...props}
              qualifiers={{ guiLabel: "Reflections" }}
            />
          </CCP4i2ContainerElement>

          {/* Search level */}
          <CCP4i2TaskElement
            itemName="SIMBAD_SEARCH_LEVEL"
            {...props}
            qualifiers={{ guiLabel: "Search level:" }}
          />

          {/* Organism */}
          <CCP4i2TaskElement
            itemName="SIMBAD_ORGANISM"
            {...props}
            qualifiers={{ guiLabel: "Organism:" }}
          />

          {/* Number of processors */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              itemName="SIMBAD_NPROC"
              {...props}
              qualifiers={{ guiLabel: "Number of processors" }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ===== Tab 2: Advanced Options ===== */}
        <CCP4i2Tab label="Advanced Options">
          <CCP4i2TaskElement
            itemName="SIMBAD_ROT_PROGRAM"
            {...props}
            qualifiers={{ guiLabel: "Rotation Search Program:" }}
          />
          <CCP4i2TaskElement
            itemName="SIMBAD_MR_PROGRAM"
            {...props}
            qualifiers={{ guiLabel: "Molecular Replacement Program:" }}
          />
          <CCP4i2TaskElement
            itemName="SIMBAD_NMOL"
            {...props}
            qualifiers={{ guiLabel: "Number of Molecules to place" }}
          />
          <CCP4i2TaskElement
            itemName="SIMBAD_PROCESS_ALL"
            {...props}
            qualifiers={{ guiLabel: "Process all possible hits:" }}
          />
          <CCP4i2TaskElement
            itemName="SIMBAD_SGALTERNATIVE"
            {...props}
            qualifiers={{ guiLabel: "Check alternative space groups:" }}
          />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
