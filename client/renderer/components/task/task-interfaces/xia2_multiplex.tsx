import { Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab label="Input data" key="input">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Search for DIALS integrated files" }}
            containerHint="FolderLevel"
            initiallyOpen={true}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="SEARCH_ROOT_DIR"
              qualifiers={{ guiLabel: "Root directory" }}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="DIALS_INTEGRATED"
              qualifiers={{ guiLabel: "DIALS .refl" }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Basic parameters" }}
            containerHint="FolderLevel"
            initiallyOpen={true}
          >
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ guiLabel: "Symmetry" }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement
                {...props}
                itemName="symmetry__space_group"
                qualifiers={{ guiLabel: "Space group" }}
              />
            </CCP4i2ContainerElement>

            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ guiLabel: "Resolution" }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement
                {...props}
                itemName="resolution__d_max"
                qualifiers={{ guiLabel: "Low resolution cutoff" }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="resolution__d_min"
                qualifiers={{ guiLabel: "High resolution cutoff" }}
              />
            </CCP4i2ContainerElement>

            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ guiLabel: "Filtering" }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement
                {...props}
                itemName="filtering__method"
                qualifiers={{ guiLabel: "Method" }}
              />
            </CCP4i2ContainerElement>

            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ guiLabel: "Clustering" }}
              containerHint="BlockLevel"
            >
              <CCP4i2TaskElement
                {...props}
                itemName="max_clusters"
                qualifiers={{ guiLabel: "Maximum number of clusters" }}
              />
              <CCP4i2TaskElement
                {...props}
                itemName="cluster_method"
                qualifiers={{ guiLabel: "Metric on which to perform clustering" }}
              />
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="Advanced parameters" key="advanced">
          <CCP4i2TaskElement
            {...props}
            itemName="controlParameters"
            qualifiers={{ guiLabel: "All parameters" }}
          />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
