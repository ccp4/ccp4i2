import {
  Box,
  Button,
  Checkbox,
  FormControlLabel,
  FormGroup,
  Radio,
  RadioGroup,
  Stack,
  Tooltip,
  Typography,
  Chip,
} from "@mui/material";
import { Add as AddIcon } from "@mui/icons-material";
import { CSimpleDataFileElement } from "./csimpledatafile";
import { CCP4i2TaskElementProps } from "./task-element";
import { useCallback, useEffect, useMemo, useState } from "react";
import { useApi } from "../../../api";
import { useJob, useProject, useProjectFiles } from "../../../utils";
import { File as DjangoFile, Project } from "../../../types/models";
import { InlineTaskModal } from "./inline-task-modal";
import { useContainerField } from "./hooks/useContainerField";

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

/** One sequence's line in the selector: name, polymer type, copies, length. */
const SequenceLabel: React.FC<{ seq: SequenceEntry }> = ({ seq }) => (
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
);

/**
 * CAsuDataFileElement - Renders a CAsuDataFile selector with sequence selection.
 *
 * The task's selectionMode qualifier (as in Qt i2) decides the selector:
 *   0 - none; the task uses every sequence
 *   2 - checkboxes, any non-empty subset; collapsed unless something is deselected
 *   1 - radio buttons, exactly one; always drawn open, so the user sees which
 *       sequence the task uses. With several sequences the task cannot run until
 *       one is chosen (the server's CAsuDataFile.validity() reports it as an error)
 * A choice updates the selection CDict on the file object.
 *
 * If no file is selected, shows a "Create ASU Content" button that opens an
 * inline modal to configure and run ProvideAsuContents, then auto-selects the output.
 */
export const CAsuDataFileElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => {
  const { job, itemName, qualifiers } = props;
  const api = useApi();
  const { useFileDigest, fileItemToParameterArg, mutateContainer } = useJob(job.id);
  const { item, unwrappedValue: value, isVisible, commit } = useContainerField<any>({
    job,
    itemName,
    visibility: props.visibility,
    disabled: props.disabled,
    onChange: props.onChange,
  });
  const { files: projectFiles } = useProjectFiles(job.project);
  const { jobs: projectJobs } = useProject(job.project);
  const { data: projects } = api.get<Project[]>("projects");

  // State for inline task modal
  const [modalOpen, setModalOpen] = useState(false);

  // Check if there are any existing CAsuDataFiles in the project
  const existingAsuFiles = useMemo(() => {
    if (!projectFiles) return [];
    return projectFiles.filter((file) => file.type === "application/CCP4-asu-content");
  }, [projectFiles]);

  // Only fetch digest when a file has been uploaded (has dbFileId)
  const hasFile = Boolean(value?.dbFileId);
  const digestPath = hasFile && item?._objectPath ? item._objectPath : "";
  const { data: fileDigest, mutate: mutateDigest } = useFileDigest(digestPath) as {
    data: CAsuDataFileDigest | undefined;
    mutate: () => void;
  };

  const overriddenQualifiers = useMemo(() => {
    return { ...item?._qualifiers, ...qualifiers };
  }, [item, qualifiers]);
  const selectionMode = Number(overriddenQualifiers.selectionMode ?? 0);

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

  // Commit a whole selection, optimistically
  const commitSelections = useCallback(
    async (updatedSelections: Record<string, boolean>) => {
      if (!item?._objectPath || job.status !== 1) return;

      const previous = localSelections;
      setLocalSelections(updatedSelections);

      const result = await commit(updatedSelections, { subPath: ".selection" });
      if (result && !result.success) {
        setLocalSelections(previous);
        return;
      }

      await Promise.all([mutateContainer(), mutateDigest()]);
    },
    [item?._objectPath, job.status, localSelections, commit, mutateContainer, mutateDigest]
  );

  // Checkbox (mode 2): toggle one sequence
  const handleSelectionChange = useCallback(
    (sequenceName: string, checked: boolean) =>
      commitSelections({ ...localSelections, [sequenceName]: checked }),
    [commitSelections, localSelections]
  );

  // Radio (mode 1): this sequence and no other
  const handleSingleSelection = useCallback(
    (sequenceName: string) =>
      commitSelections(
        Object.fromEntries(
          (fileDigest?.sequences ?? []).map((seq) => [seq.name, seq.name === sequenceName])
        )
      ),
    [commitSelections, fileDigest?.sequences]
  );

  // Mode 1's choice, or "" while none (or several) are selected
  const singleSelection = useMemo(() => {
    const chosen = Object.entries(localSelections).filter(([, sel]) => sel);
    return chosen.length === 1 ? chosen[0][0] : "";
  }, [localSelections]);

  // Determine if we should force the panel expanded
  const forceExpanded = useMemo(() => {
    if (!hasSequences) return false;
    // Mode 1: always open, so the user sees which sequence the task uses
    if (selectionMode === 1) return true;
    // Mode 2: expand if any sequence is deselected (non-default state)
    if (selectionMode === 2)
      return Object.values(localSelections).some((selected) => !selected);
    return false;
  }, [hasSequences, selectionMode, localSelections]);

  /**
   * Handle the output file from the inline ProvideAsuContents task.
   * Auto-selects it in this widget via setParameter.
   */
  const handleOutputFileReady = useCallback(
    async (outputFile: DjangoFile) => {
      if (!item?._objectPath || !projectJobs) return;

      const paramArg = fileItemToParameterArg(
        outputFile,
        item._objectPath,
        projectJobs,
        projects || []
      );

      await commit(paramArg.value);
      await Promise.all([mutateContainer(), mutateDigest()]);
    },
    [item?._objectPath, projectJobs, projects, fileItemToParameterArg, commit, mutateContainer, mutateDigest]
  );

  if (!isVisible) return null;

  // Determine if we should show the expanded panel (for create button when no files exist)
  const shouldForceExpand = forceExpanded || (!hasFile && existingAsuFiles.length === 0);

  return (
    <>
      <CSimpleDataFileElement {...props} forceExpanded={shouldForceExpand}>
        {/* Create ASU Content action - shown when no file selected */}
        {!hasFile && job.status === 1 && (
          <Box sx={{ mb: hasSequences ? 2 : 0 }}>
            <Stack direction="row" spacing={1} alignItems="center">
              <Tooltip title="Create a new ASU content file using the ProvideAsuContents task">
                <Button
                  variant="outlined"
                  size="small"
                  startIcon={<AddIcon />}
                  onClick={() => setModalOpen(true)}
                  sx={{ textTransform: "none" }}
                >
                  Create ASU Content
                </Button>
              </Tooltip>
              {existingAsuFiles.length === 0 && (
                <Typography variant="caption" color="text.secondary">
                  No ASU content files in project
                </Typography>
              )}
            </Stack>
          </Box>
        )}

        {/* Sequence selection - its form follows the task's selectionMode */}
        {hasSequences && (selectionMode === 1 || selectionMode === 2) && (
          <Stack spacing={1} sx={{ mt: 1 }}>
            <Typography variant="subtitle2" color="text.secondary">
              {selectionMode === 1 ? "Select one sequence" : "Select one or more sequences"}
            </Typography>
            {selectionMode === 1 ? (
              <RadioGroup
                value={singleSelection}
                onChange={(e) => handleSingleSelection(e.target.value)}
              >
                {fileDigest?.sequences?.map((seq) => (
                  <FormControlLabel
                    key={seq.index}
                    value={seq.name}
                    control={<Radio size="small" disabled={job.status !== 1} />}
                    label={<SequenceLabel seq={seq} />}
                    sx={{ ml: 0 }}
                  />
                ))}
              </RadioGroup>
            ) : (
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
                    label={<SequenceLabel seq={seq} />}
                    sx={{ ml: 0 }}
                  />
                ))}
              </FormGroup>
            )}
          </Stack>
        )}
      </CSimpleDataFileElement>

      {/* Inline modal for creating ProvideAsuContents */}
      <InlineTaskModal
        open={modalOpen}
        onClose={() => setModalOpen(false)}
        taskName="ProvideAsuContents"
        parentJob={job}
        onOutputFileReady={handleOutputFileReady}
        title="Create ASU Contents"
      />
    </>
  );
};
