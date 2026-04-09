import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";
import { useCallback, useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container, fetchDigest, useFileDigest } = useJob(
    props.job.id
  );
  const { value: alignMode } = useTaskItem("ALIGNMENTORSEQUENCEIN");
  const { item: alignItem } = useTaskItem("ALIGNIN");
  const { forceUpdate: setTargetIndex } = useTaskItem("TARGETINDEX");

  // Fetch alignment digest to get sequence identifiers
  const { data: alignDigest } = useFileDigest(alignItem?._objectPath || "");

  // Build enumerators and labels from the alignment digest identifiers
  const { enumerators, menuText } = useMemo(() => {
    const ids: string[] = alignDigest?.identifiers || [];
    if (ids.length === 0) return { enumerators: [], menuText: [] };
    return {
      enumerators: ids.map((_: string, i: number) => i),
      menuText: ids,
    };
  }, [alignDigest?.identifiers]);

  // When ALIGNIN changes, fetch its digest to update identifiers
  const handleAlignInChange = useCallback(async () => {
    if (!alignItem?._objectPath) return;
    const digest = await fetchDigest(alignItem._objectPath);
    if (digest?.identifiers?.length > 0) {
      await setTargetIndex(0);
    }
  }, [alignItem?._objectPath, fetchDigest, setTargetIndex]);

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
        <CCP4i2TaskElement itemName="XYZIN" {...props} />
        {alignMode === "SEQUENCE" && (
          <CCP4i2TaskElement
            itemName="CHAINIDS"
            {...props}
            qualifiers={{ guiLabel: "Chain ID in template to manipulate" }}
          />
        )}
        <CCP4i2TaskElement
          itemName="ALIGNMENTORSEQUENCEIN"
          {...props}
          qualifiers={{ guiLabel: "Provide target sequence as:" }}
        />
        {alignMode === "ALIGNMENT" && (
          <CCP4i2TaskElement
            itemName="ALIGNIN"
            {...props}
            onChange={handleAlignInChange}
          />
        )}
        {alignMode === "ALIGNMENT" && enumerators.length > 0 && (
          <CCP4i2TaskElement
            itemName="TARGETINDEX"
            {...props}
            qualifiers={{
              guiLabel:
                "Identifier of the target sequence in this alignment",
              enumerators,
              menuText,
              onlyEnumerators: true,
            }}
          />
        )}
        {alignMode === "SEQUENCE" && (
          <CCP4i2TaskElement itemName="SEQUENCEIN" {...props} />
        )}
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
