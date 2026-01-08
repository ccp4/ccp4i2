import { LinearProgress, Paper, Stack, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useJob } from "../../../utils";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useCallback } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem, fetchDigest } = useJob(job.id);

  const { forceUpdate: forceSetSEQUENCETEXT } = useTaskItem("SEQUENCETEXT");
  const { item: SEQINItem } = useTaskItem("SEQIN");

  // Handle SEQIN file change - explicitly fetch digest and populate SEQUENCETEXT
  const handleSeqInChange = useCallback(async () => {
    if (!SEQINItem?._objectPath) return;

    const digestData = await fetchDigest(SEQINItem._objectPath);
    if (!digestData?.sequence) return;

    // Build the formatted sequence text
    const newSequence = `>${digestData.identifier || ""}\n${digestData.sequence}`.replace("*", "");
    await forceSetSEQUENCETEXT(newSequence);
  }, [SEQINItem?._objectPath, fetchDigest, forceSetSEQUENCETEXT]);

  return (
    <CCP4i2Tabs {...props}>
      <CCP4i2Tab label="Main inputs">
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ guiLabel: "Key files" }}
          containerHint="FolderLevel"
          initiallyOpen={true}
        >
          <CCP4i2TaskElement
            {...props}
            itemName="SEQUENCETEXT"
            qualifiers={{ guiLabel: "Sequence", guiMode: "multiLine" }}
            sx={{ minWidth: "100%", minHeight: "10rem" }}
          />

          <CCP4i2TaskElement
            {...props}
            itemName="SEQIN"
            qualifiers={{ guiLabel: "File from which to extract sequence" }}
            onChange={handleSeqInChange}
          />

          <CCP4i2TaskElement
            {...props}
            itemName="XYZIN"
            qualifiers={{ guiLabel: "MTZFile (for Matthews volumne calc)" }}
          />
        </CCP4i2ContainerElement>
      </CCP4i2Tab>
    </CCP4i2Tabs>
  );
};
export default TaskInterface;
