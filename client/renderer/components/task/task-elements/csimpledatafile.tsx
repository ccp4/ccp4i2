import { CDataFileElement, IconMenuItem } from "./cdatafile";
import { CCP4i2TaskElementProps } from "./task-element";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useJob, useProject } from "../../../utils";
import { useImportProvenance } from "../../../providers/import-provenance-provider";

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
  const { requestImportProvenance } = useImportProvenance();
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

    // Ask for a provenance note first (a no-op unless the user has turned the
    // preference on), so it rides along in the same upload POST. null means
    // "don't attach"; "" means the user chose Skip.
    const provenance = await requestImportProvenance(selectedFiles[0].name);

    // Hand over the picked File itself, not a re-read copy: a File is a Blob,
    // so the upload needs nothing more, and only a real File lets the desktop
    // import it by path (no bytes through the browser) and lets a served
    // deployment stage it in chunks. Re-wrapping in a Blob defeated both and
    // read the whole file into memory first.
    const uploadResult = await uploadFileParam({
      objectPath: item._objectPath,
      file: selectedFiles[0],
      fileName: selectedFiles[0].name,
      description: provenance ?? undefined,
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
    requestImportProvenance,
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
