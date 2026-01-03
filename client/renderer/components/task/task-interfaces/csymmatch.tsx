import React from "react";
import { Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";

/**
 * Task interface component for CSymmatch - Coordinate Symmetry Matching.
 *
 * CSymmatch is used to:
 * - Align and compare protein structures by finding symmetry relationships
 * - Move query coordinates to match reference coordinates
 * - Handle origin shifts and hand changes for optimal alignment
 * - Define connectivity parameters for structural linkage analysis
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab label="Main inputs" key="main">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Input coordinates",
              initiallyOpen: true,
            }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="XYZIN_QUERY"
              qualifiers={{
                guiLabel: "Coordinates to move",
                toolTip:
                  "Query structure that will be moved to align with the reference",
              }}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="XYZIN_TARGET"
              qualifiers={{
                guiLabel: "Reference coordinates",
                toolTip:
                  "Target structure that serves as the reference for alignment",
              }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Alignment options",
              initiallyOpen: true,
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              itemName="ORIGIN_HAND"
              qualifiers={{
                guiLabel: "Try all possible origin shifts, and changes of hand",
                toolTip:
                  "Explore all symmetry operations including origin shifts and enantiomorph changes",
              }}
            />

            <CCP4i2TaskElement
              {...props}
              itemName="CONNECTIVITY_RADIUS"
              qualifiers={{
                guiLabel: "Radius to define linkage",
                toolTip:
                  "Distance threshold for determining atomic connectivity in Angstroms",
              }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
