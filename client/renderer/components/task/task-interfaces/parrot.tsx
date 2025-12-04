import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useApi } from "../../../api";
import { useJob } from "../../../utils";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { Button } from "@mui/material";
import { useRunCheck } from "../../../providers/run-check-provider";
import { useCallback, useContext, useEffect } from "react";
import { useRouter } from "next/navigation";
import { Job } from "../../../types/models";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const api = useApi();
  const { job } = props;
  const { useTaskItem, createPeerTask } = useJob(job.id);
  const { value: XYZIN_MODE } = useTaskItem("XYZIN_MODE");
  const { value: ASUIN } = useTaskItem("ASUIN");
  const router = useRouter();

  // 1. Retrieve the function for installing extraDialogActions from the relevant context
  // layer

  const { extraDialogActions, setExtraDialogActions, setRunTaskRequested } =
    useRunCheck();

  // 2. Define a callback to create an ASUIN task if it is not already set

  const createAsuTask = useCallback(async () => {
    await createPeerTask("ProvideAsuContents").then((created_job: Job) => {
      if (created_job) {
        // If the task was created successfully, we can navigate to it
        router.push(`/project/${job.project}/job/${created_job.id}`);
        //Shut down the run check dialog
        setRunTaskRequested(null);
      }
    });
  }, [job, createPeerTask]);

  // 3. Use a useEffect to install the extraDialogActions

  useEffect(() => {
    if (!(ASUIN?.dbFileId?.length > 0)) {
      // If the ASUIN is not set, we add an action to create a ProvideAsuContents task
      // This will be shown in the confirm dialog.  As ever when changing state,
      // we check if the action is already there to avoid unnecessary re-renders.
      if (!extraDialogActions || !extraDialogActions["ASUIN"]) {
        const newExtraDialogActions = {
          ASUIN: (
            <Button variant="contained" onClick={createAsuTask}>
              Create ASU Contents task
            </Button>
          ),
        };
        setExtraDialogActions(newExtraDialogActions);
      }
    }
  }, [ASUIN, extraDialogActions, setExtraDialogActions]);

  //This is a really important cleanup function to avoid memory leaks
  //It ensures that processedErrors and extraDialogActions are cleared when the component unmounts
  useEffect(() => {
    // Cleanup function to reset context values when the component unmounts
    return () => {
      setExtraDialogActions(null);
    };
  }, [setExtraDialogActions]);

  return (
    <CCP4i2Tabs {...props}>
      <CCP4i2Tab label="Main inputs">
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ guiLabel: "Input data" }}
          containerHint="FolderLevel"
          initiallyOpen={true}
        >
          <CCP4i2TaskElement
            {...props}
            itemName="F_SIGF"
            qualifiers={{ guiLabel: "Reflections" }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="FREERFLAG"
            qualifiers={{ guiLabel: "Free R set" }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="ASUIN"
            qualifiers={{ guiLabel: "Asymmetric unit contents" }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="ABCD"
            qualifiers={{ guiLabel: "Starting phases" }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="F_PHI"
            qualifiers={{ guiLabel: "Starting map" }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="XYZIN_MODE"
            qualifiers={{ guiLabel: "Use of NCS", guiMode: "radio" }}
          />
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Model from which to infer NCS" }}
            containerHint="BlockLevel"
            visibility={() => XYZIN_MODE !== "no"}
          >
            <CCP4i2TaskElement
              {...props}
              itemName="XYZIN_HA"
              qualifiers={{ guiLabel: "Heavy atom model" }}
              visibility={() => XYZIN_MODE === "ha"}
            />
            <CCP4i2TaskElement
              {...props}
              itemName="XYZIN_MR"
              qualifiers={{ guiLabel: "Full atom model" }}
              visibility={() => XYZIN_MODE === "mr"}
            />
          </CCP4i2ContainerElement>
        </CCP4i2ContainerElement>

        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ guiLabel: "Controls" }}
          containerHint="FolderLevel"
          initiallyOpen={true}
        >
          <CCP4i2TaskElement
            {...props}
            itemName="CYCLES"
            qualifiers={{ guiLabel: "Number of cycles" }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="ANISOTROPY_CORRECTION"
            qualifiers={{ guiLabel: "Apply anisotropy correction" }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="RESOLUTION"
            qualifiers={{ guiLabel: "Maximum resolution" }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="SOLVENT_CONTENT"
            qualifiers={{ guiLabel: "Estimated solvent content" }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="NCS_MASK_FILTER_RADIUS"
            qualifiers={{ guiLabel: "Filter radius to define NCS mask" }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="VERBOSE"
            qualifiers={{ guiLabel: "Verbosity of log file" }}
          />
        </CCP4i2ContainerElement>
      </CCP4i2Tab>
      <CCP4i2Tab label="Reference structures" key="2">
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ guiLabel: "Reference density and atmomic models" }}
          containerHint="FolderLevel"
          initiallyOpen={true}
        >
          <span>
            <b>
              <em>
                You should normally let Parrot choose reference structures
              </em>
            </b>
          </span>
          <CCP4i2TaskElement
            {...props}
            itemName="F_SIGF_REF"
            qualifiers={{ guiLabel: "Reference density reflections" }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="ABCD_REF"
            qualifiers={{ guiLabel: "Reference density phases" }}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="XYZIN_REF"
            qualifiers={{
              guiLabel: "Model to define region of space to build",
            }}
          />
        </CCP4i2ContainerElement>
      </CCP4i2Tab>
    </CCP4i2Tabs>
  );
};

export default TaskInterface;
