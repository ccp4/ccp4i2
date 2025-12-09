import {
  Checkbox,
  FormControlLabel,
  FormGroup,
  Stack,
  Typography,
  Chip,
} from "@mui/material";
import { CSimpleDataFileElement } from "./csimpledatafile";
import { CCP4i2TaskElementProps } from "./task-element";
import { useCallback, useEffect, useMemo, useState } from "react";
import { useJob } from "../../../utils";

/**
 * Sequence entry from the CAsuDataFile digest
 */
interface SequenceEntry {
  index: number;
  name: string;
  polymerType: string;
  nCopies: number;
  sequenceLength: number;
  sequencePreview: string;
  description: string;
  selected: boolean;
}

/**
 * Digest structure returned by the backend for CAsuDataFile
 */
interface CAsuDataFileDigest {
  sequences?: SequenceEntry[];
  sequenceCount?: number;
  status?: string;
  reason?: string;
}

/**
 * CAsuDataFileElement - Renders a CAsuDataFile selector with sequence selection checkboxes.
 *
 * When a file is selected, the digest is fetched and the sequences are displayed
 * as checkboxes. Users can select/deselect individual sequences, which updates
 * the selection CDict on the file object.
 */
export const CAsuDataFileElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => {
  const { job, itemName, qualifiers } = props;
  const { useTaskItem, useFileDigest, setParameter, mutateContainer } = useJob(job.id);
  const { item, value } = useTaskItem(itemName);

  // Only fetch digest when a file has been uploaded (has dbFileId)
  const hasFile = Boolean(value?.dbFileId);
  const digestPath = hasFile && item?._objectPath ? item._objectPath : "";
  const { data: fileDigest, mutate: mutateDigest } = useFileDigest(digestPath) as {
    data: CAsuDataFileDigest | undefined;
    mutate: () => void;
  };

  // Local state for checkbox values (optimistic updates)
  const [localSelections, setLocalSelections] = useState<Record<string, boolean>>({});

  // Track if we have selections to show
  const hasSequences = useMemo(() => {
    return fileDigest?.sequences && fileDigest.sequences.length > 0;
  }, [fileDigest]);

  // Initialize local selections from digest
  useEffect(() => {
    if (fileDigest?.sequences) {
      const initialSelections: Record<string, boolean> = {};
      fileDigest.sequences.forEach((seq) => {
        initialSelections[seq.name] = seq.selected;
      });
      setLocalSelections(initialSelections);
    }
  }, [fileDigest?.sequences]);

  // Handle checkbox change
  const handleSelectionChange = useCallback(
    async (sequenceName: string, checked: boolean) => {
      if (!item?._objectPath || job.status !== 1) return;

      // Optimistic update
      setLocalSelections((prev) => ({
        ...prev,
        [sequenceName]: checked,
      }));

      // Update the selection CDict on the backend
      // The selection is stored at item._objectPath.selection[sequenceName]
      const selectionPath = `${item._objectPath}.selection`;

      try {
        // We need to update the CDict - use a special format that the backend understands
        // The CDict stores key-value pairs, so we send the full updated selection
        const updatedSelections = {
          ...localSelections,
          [sequenceName]: checked,
        };

        await setParameter({
          object_path: selectionPath,
          value: updatedSelections,
        });

        // Refresh data
        await Promise.all([mutateContainer(), mutateDigest()]);
      } catch (error) {
        console.error("Error updating sequence selection:", error);
        // Revert optimistic update on error
        setLocalSelections((prev) => ({
          ...prev,
          [sequenceName]: !checked,
        }));
      }
    },
    [item?._objectPath, job.status, localSelections, setParameter, mutateContainer, mutateDigest]
  );

  // Determine if we should force the panel expanded
  const forceExpanded = useMemo(() => {
    if (!hasSequences) return false;
    // Expand if any sequence is deselected (non-default state)
    return Object.values(localSelections).some((selected) => !selected);
  }, [hasSequences, localSelections]);

  // Override qualifiers to enable selectionMode display
  const overriddenQualifiers = useMemo(() => {
    return { ...item?._qualifiers, ...qualifiers };
  }, [item, qualifiers]);

  // Check visibility
  const inferredVisibility = useMemo(() => {
    if (!props.visibility) return true;
    if (typeof props.visibility === "function") {
      return props.visibility();
    }
    return props.visibility;
  }, [props.visibility]);

  if (!inferredVisibility) return null;

  return (
    <CSimpleDataFileElement {...props} forceExpanded={forceExpanded}>
      {hasSequences && overriddenQualifiers.selectionMode !== undefined && (
        <Stack spacing={1} sx={{ mt: 1 }}>
          <Typography variant="subtitle2" color="text.secondary">
            Sequence Selection
          </Typography>
          <FormGroup>
            {fileDigest?.sequences?.map((seq) => (
              <FormControlLabel
                key={seq.index}
                control={
                  <Checkbox
                    checked={localSelections[seq.name] ?? seq.selected}
                    onChange={(e) =>
                      handleSelectionChange(seq.name, e.target.checked)
                    }
                    disabled={job.status !== 1}
                    size="small"
                  />
                }
                label={
                  <Stack direction="row" spacing={1} alignItems="center">
                    <Typography variant="body2" fontWeight="medium">
                      {seq.name}
                    </Typography>
                    <Chip
                      label={seq.polymerType}
                      size="small"
                      variant="outlined"
                      sx={{ height: 20, fontSize: "0.7rem" }}
                    />
                    <Typography variant="caption" color="text.secondary">
                      {seq.nCopies} {seq.nCopies === 1 ? "copy" : "copies"} |{" "}
                      {seq.sequenceLength} residues
                    </Typography>
                    {seq.description && (
                      <Typography
                        variant="caption"
                        color="text.secondary"
                        sx={{
                          maxWidth: 200,
                          overflow: "hidden",
                          textOverflow: "ellipsis",
                          whiteSpace: "nowrap",
                        }}
                      >
                        - {seq.description}
                      </Typography>
                    )}
                  </Stack>
                }
                sx={{ ml: 0 }}
              />
            ))}
          </FormGroup>
        </Stack>
      )}
    </CSimpleDataFileElement>
  );
};
