import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob, usePrevious } from "../../../utils";
import { useEffect, useRef } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { container, useFileDigest, useTaskItem } = useJob(job.id);

  const { value: HKLINValue } = useTaskItem("HKLIN");
  const oldHKLINValue = usePrevious(HKLINValue);
  const hasUploadedFile = Boolean(HKLINValue?.dbFileId);
  const digestObjectPath = hasUploadedFile ? "cif2mtz.inputData.HKLIN" : "";
  const { data: HKLINDigest, mutate: mutateDigest } =
    useFileDigest(digestObjectPath);

  const { forceUpdate: forceUpdateSPACEGROUPCELL } =
    useTaskItem("SPACEGROUPCELL");

  // Track which digest we've already processed to avoid re-running
  const processedDigestRef = useRef<string | null>(null);

  // When the file changes, invalidate the cached digest so SWR re-fetches
  useEffect(() => {
    if (!HKLINValue?.dbFileId || !oldHKLINValue?.dbFileId) return;
    if (HKLINValue.dbFileId !== oldHKLINValue.dbFileId) {
      processedDigestRef.current = null;
      mutateDigest();
    }
  }, [HKLINValue?.dbFileId, oldHKLINValue?.dbFileId, mutateDigest]);

  // Auto-populate spacegroup and cell from file digest
  useEffect(() => {
    if (!HKLINDigest || job?.status !== 1) return;

    const digestKey = JSON.stringify({
      spaceGroup: HKLINDigest.spaceGroup,
      cell: HKLINDigest.cell,
    });

    if (digestKey === processedDigestRef.current) return;
    processedDigestRef.current = digestKey;

    const processDigest = async () => {
      if (!forceUpdateSPACEGROUPCELL) return;

      const update: Record<string, any> = {};

      if (HKLINDigest.spaceGroup) {
        update.spaceGroup = String(HKLINDigest.spaceGroup).replace(
          /\s+/g,
          ""
        );
      }

      if (HKLINDigest.cell) {
        update.cell = HKLINDigest.cell;
      }

      if (Object.keys(update).length > 0) {
        try {
          await forceUpdateSPACEGROUPCELL(update);
        } catch (e) {
          console.error("[cif2mtz] Failed to update SPACEGROUPCELL:", e);
        }
      }
    };

    processDigest();
  }, [HKLINDigest, job?.status]);

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Input Data" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="HKLIN" {...props} />
        <CCP4i2TaskElement itemName="SPACEGROUPCELL" {...props} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
