import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* Single folder: Xia2 runs */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Xia2 runs" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="XIA2_DIRECTORY"
          {...props}
          qualifiers={{
            guiLabel: "xia2 run directory, or a directory of runs",
            toolTip: "The directory of one xia2 run (it holds DataFiles), or one holding several runs",
          }}
        />
        <Typography variant="body2" color="text.secondary" sx={{ mt: 1 }}>
          The runs are found in it: the directory itself if it is a xia2 run,
          otherwise each sub-directory that is. Each run&apos;s merged data, free
          set and unmerged integrated reflections are imported.
        </Typography>
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
