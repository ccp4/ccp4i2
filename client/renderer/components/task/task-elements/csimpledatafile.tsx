import { CDataFileElement, IconMenuItem } from "./cdatafile";
import { CCP4i2TaskElementProps } from "./task-element";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { readFilePromise, useJob, useProject } from "../../../utils";

interface CSimpleDataFileElementProps extends CCP4i2TaskElementProps {
  hasValidationError?: boolean;
  forceExpanded?: boolean;
  iconMenuItems?: IconMenuItem[];
}

export const CSimpleDataFileElement: React.FC<CSimpleDataFileElementProps> = (
  props
) => {
  const { job, itemName, onChange, visibility } = props;
  const { useTaskItem, useFileDigest, uploadFileParam } = useJob(job.id);
  const { mutateFiles, mutateJobs } = useProject(job.project);
  const { item } = useTaskItem(itemName);
  const [selectedFiles, setSelectedFiles] = useState<FileList | null>(null);
  const { data: fileDigest, mutate: mutateDigest } = useFileDigest(
    item?._objectPath
  );
  const previousSelectedFiles = useRef<FileList | null>(null);

  const processFirstFile = useCallback(async () => {
    if (!selectedFiles || selectedFiles.length == 0 || !item) return;
    if (selectedFiles === previousSelectedFiles.current) return;
    previousSelectedFiles.current = selectedFiles;
    const fileBuffer = await readFilePromise(selectedFiles[0], "ArrayBuffer");

    // Use centralized uploadFileParam with local cache patching
    const uploadResult = await uploadFileParam({
      objectPath: item._objectPath,
      file: new Blob([fileBuffer as ArrayBuffer], { type: item._qualifiers.mimeTypeName }),
      fileName: selectedFiles[0].name,
    });

    // Handle response
    if (uploadResult?.success && uploadResult.data?.updated_item) {
      onChange?.(uploadResult.data.updated_item);
    }
    setSelectedFiles(null);

    // Execute additional mutations not handled by uploadFileParam
    await Promise.all([
      mutateJobs(),
      mutateFiles(),
      mutateDigest(),
    ]);
  }, [
    item,
    selectedFiles,
    onChange,
    uploadFileParam,
    mutateJobs,
    mutateFiles,
    mutateDigest,
  ]);

  // Auto-process files when selected.
  //
  // processFirstFile is async; previously it was fired here un-awaited and with
  // no error handler, so ANY failure (file read, upload POST, digest fetch) was
  // silently swallowed — the file picker appeared to do nothing at all, with no
  // console message. Attach a catch so the failure is at least visible/diagnosable
  // and the component is left in a clean state to retry.
  useEffect(() => {
    if (selectedFiles && processFirstFile) {
      processFirstFile().catch((err) => {
        console.error(
          `File upload failed for ${itemName} ` +
            `(${selectedFiles?.[0]?.name ?? "unknown file"}):`,
          err
        );
        setSelectedFiles(null);
      });
    }
  }, [selectedFiles, processFirstFile, itemName]);

  const isVisible = useMemo(
    () =>
      !visibility ||
      (typeof visibility === "function" ? visibility() : visibility),
    [visibility]
  );

  if (!isVisible) return null;

  return <CDataFileElement {...props} setFiles={setSelectedFiles} />;
};
