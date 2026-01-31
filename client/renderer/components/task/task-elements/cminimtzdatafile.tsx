import { Stack } from "@mui/material";
import { CDataFileElement } from "./cdatafile";
import { CCP4i2TaskElementProps } from "./task-element";
import { useCallback, useMemo } from "react";
import { BaseSpacegroupCellElement } from "./base-spacegroup-cell-element";
import { readFilePromise, useJob, useProject } from "../../../utils";
import { useCCP4i2Window } from "../../../app-context";
import { selectMtzColumnsEnhanced, SiblingInput } from "./mtz-column-dialog";

/** MTZ-related class names that are siblings of interest */
const MTZ_SIBLING_CLASSES = [
  "CFreeRDataFile",
  "CObsDataFile",
  "CMapCoeffsDataFile",
  "CMiniMtzDataFile",
  "CMtzDataFile",
];

export const CMiniMtzDataFileElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => {
  const { job, itemName, onChange, visibility } = props;
  const { useTaskItem, useFileDigest, uploadFileParam, container } = useJob(job.id);
  const { mutateJobs, mutateFiles } = useProject(job.project);
  const { item, value } = useTaskItem(itemName);
  const { cootModule } = useCCP4i2Window();

  // Only fetch digest when a file has been uploaded (has dbFileId)
  const hasFile = Boolean(value?.dbFileId);
  const digestPath = hasFile && item?._objectPath ? item._objectPath : "";
  const { data: fileDigest, mutate: mutateDigest } = useFileDigest(digestPath);

  const infoContent = useMemo(
    () => <BaseSpacegroupCellElement data={fileDigest} />,
    [fileDigest]
  );

  /**
   * Get sibling inputs from the same parent container.
   * Siblings are items under the same inputData container (e.g., task.inputData.*)
   * that are MTZ-related data files.
   */
  const getSiblingInputs = useCallback((): SiblingInput[] => {
    if (!item?._objectPath || !container?.lookup) {
      return [];
    }

    // Extract the parent path (e.g., "task.inputData" from "task.inputData.F_SIGF")
    const parts = item._objectPath.split(".");
    if (parts.length < 2) return [];

    // Remove the last part to get the parent path
    const parentPath = parts.slice(0, -1).join(".");

    const siblings: SiblingInput[] = [];

    // Search through container lookup for sibling items
    for (const [objectPath, lookupItem] of Object.entries(container.lookup)) {
      // Check if this item is a sibling (same parent, different item)
      if (
        objectPath.startsWith(parentPath + ".") &&
        objectPath !== item._objectPath &&
        lookupItem &&
        typeof lookupItem === "object" &&
        "_class" in lookupItem &&
        MTZ_SIBLING_CLASSES.includes(lookupItem._class as string)
      ) {
        siblings.push({
          objectPath,
          className: lookupItem._class as string,
        });
      }
    }

    return siblings;
  }, [item?._objectPath, container?.lookup]);

  /**
   * Handle file selection from the file picker.
   * This is the single entry point for file uploads - it handles:
   * 1. Showing the column selection dialog (for MTZ files) with sibling awareness
   * 2. Uploading the file with the selected columns (with local cache patching)
   * 3. Updating the UI state
   */
  const handleFileSelection = useCallback(
    async (files: FileList | null) => {
      if (!files || files.length === 0 || !item || !cootModule) {
        return;
      }

      const file = files[0];

      try {
        // Get sibling inputs for context-aware dialog
        const siblingInputs = getSiblingInputs();

        // Show enhanced column selection dialog with sibling awareness
        const result = await selectMtzColumnsEnhanced({
          file,
          item,
          cootModule,
          siblingInputs,
          multiSelectMode: false, // Single-select for now, can be enabled per-task
        });

        // User cancelled the dialog
        if (result === null) {
          return;
        }

        // Read file and upload using centralized uploadFileParam (with local cache patching)
        const fileBuffer = await readFilePromise(file, "ArrayBuffer");
        const fileBlob = new Blob([fileBuffer as ArrayBuffer], { type: "application/CCP4-mtz-file" });

        // Use enhanced columnSelectors if available, otherwise fall back to single columnSelector
        const uploadResult = await uploadFileParam({
          objectPath: item._objectPath,
          file: fileBlob,
          fileName: file.name,
          // Send both for backward compatibility
          columnSelector: result.columnSelector || undefined,
          columnSelectors: result.reflectionSelections,
        });

        // Handle response
        if (uploadResult?.success && uploadResult.data?.updated_item) {
          onChange?.(uploadResult.data.updated_item);
        }

        // If FreeR was selected, upload to the sibling CFreeRDataFile input
        if (result.freeRSelection) {
          const freeRSibling = siblingInputs.find(
            (sibling) => sibling.className === "CFreeRDataFile"
          );

          if (freeRSibling) {
            console.log("Uploading FreeR to sibling:", freeRSibling.objectPath, result.freeRSelection);
            await uploadFileParam({
              objectPath: freeRSibling.objectPath,
              file: fileBlob,
              fileName: file.name,
              columnSelector: result.freeRSelection.columnSelector,
            });
          }
        }

        // Trigger additional mutations not handled by uploadFileParam
        await Promise.all([
          mutateJobs(),
          mutateFiles(),
          mutateDigest(),
        ]);
      } catch (error) {
        console.error("File upload failed:", error);
        // Could show an error toast/snackbar here
      }
    },
    [item, cootModule, getSiblingInputs, onChange, uploadFileParam, mutateJobs, mutateFiles, mutateDigest]
  );

  const isVisible = useMemo(
    () =>
      !visibility ||
      (typeof visibility === "function" ? visibility() : visibility),
    [visibility]
  );

  if (!isVisible) return null;

  // Wrap async handler for the sync setFiles prop
  const setFiles = useCallback(
    (files: FileList | null) => {
      handleFileSelection(files);
    },
    [handleFileSelection]
  );

  return (
    <Stack direction="column">
      <CDataFileElement
        {...props}
        infoContent={infoContent}
        setFiles={setFiles}
      />
    </Stack>
  );
};
