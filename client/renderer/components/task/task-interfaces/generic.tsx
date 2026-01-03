import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { useApi } from "../../../api";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const api = useApi();
  const { data: container, mutate: mutateContainer } =
    api.get_wrapped_endpoint_json<any>({
      type: "jobs",
      id: props.job.id,
      endpoint: "container",
    });

  if (!container) return <LinearProgress />;

  return (
    <Paper>
      <CCP4i2Tabs>
        <CCP4i2Tab key="Input data" label="Input data">
          <CCP4i2TaskElement
            itemName="inputData"
            {...props}
            qualifiers={{ guiLabel: "Input data" }}
          />
        </CCP4i2Tab>
        <CCP4i2Tab key="Parameters" label="Parameters">
          <CCP4i2TaskElement
            itemName="controlParameters"
            {...props}
            qualifiers={{ guiLabel: "Parameters" }}
          />
        </CCP4i2Tab>
        <CCP4i2Tab key="Keywords" label="Keywords">
          <CCP4i2TaskElement
            itemName="keywords"
            {...props}
            qualifiers={{ guiLabel: "Parameters" }}
          />
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
