import React from "react";
import { Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Atomic model to be moved" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...props}
          itemName="XYZIN_QUERY"
          qualifiers={{ guiLabel: "Atomic model" }}
        />
      </CCP4i2ContainerElement>

      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Fixed (reference) atomic model" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...props}
          itemName="XYZIN_TARGET"
          qualifiers={{ guiLabel: "Atomic model" }}
        />
      </CCP4i2ContainerElement>

      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Options" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...props}
          itemName="ORIGIN_HAND"
          qualifiers={{ guiLabel: "Try all possible origins and hands" }}
        />

        <CCP4i2TaskElement
          {...props}
          itemName="CONNECTIVITY_RADIUS"
          qualifiers={{
            guiLabel:
              "Radius to use in stiching floating fragments to chains",
          }}
        />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
