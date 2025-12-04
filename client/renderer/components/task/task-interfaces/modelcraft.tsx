import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useApi } from "../../../api";
import { useJob } from "../../../utils";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useCallback } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const api = useApi();
  const { job } = props;
  const { useTaskItem } = useJob(job.id);
  const { value: USE_MODEL_PHASES } = useTaskItem("USE_MODEL_PHASES");
  const { value: XYZIN, update: setXYZIN } = useTaskItem("XYZIN");

  const handleUSE_MODEL_PHASES = useCallback(
    async (new_USE_MODEL_PHASES: any) => {
      //Look at the return value of the USE_MODEL_PHASES item
      //Clear XYZIN if updated USE_MODEL_PHASES is false
      if (!new_USE_MODEL_PHASES._value && XYZIN?.dbFileId) {
        setXYZIN({});
      }
    },
    [XYZIN, setXYZIN]
  );

  return (
    <>
      <CCP4i2Tabs {...props}>
        <CCP4i2Tab label="Main inputs" key="1">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Input data", initiallyOpen: true }}
            key="Input data"
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              key="F_SIGF"
              itemName="F_SIGF"
              qualifiers={{ guiLabel: "Reflections" }}
            />
            <CCP4i2TaskElement
              {...props}
              key="FREERFLAG"
              itemName="FREERFLAG"
              qualifiers={{ guiLabel: "Free R set" }}
            />
            <CCP4i2TaskElement
              {...props}
              key="ASUIN"
              itemName="ASUIN"
              qualifiers={{ guiLabel: "Asymmetric unit contents" }}
            />
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              key="Starting phases"
              qualifiers={{ guiLabel: "Starting phases" }}
            >
              <CCP4i2TaskElement
                {...props}
                key="USE_MODEL_PHASES"
                itemName="USE_MODEL_PHASES"
                qualifiers={{ guiLabel: "Use model phases" }}
                onChange={handleUSE_MODEL_PHASES}
              />
              <CCP4i2TaskElement
                {...props}
                key="PHASES"
                itemName="PHASES"
                qualifiers={{ guiLabel: "Starting phases" }}
                visibility={() => {
                  return !USE_MODEL_PHASES;
                }}
              />
              <CCP4i2TaskElement
                {...props}
                key="XYZIN"
                itemName="XYZIN"
                qualifiers={{ guiLabel: "Coordinates" }}
                visibility={() => {
                  return USE_MODEL_PHASES;
                }}
              />
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Controls" }}
            key="Controls"
          >
            <CCP4i2TaskElement
              {...props}
              key="SELENOMET"
              itemName="SELENOMET"
              qualifiers={{
                guiLabel: "Build methionine (MET) as selenomethionine (MSE)",
              }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </>
  );
};

export default TaskInterface;
