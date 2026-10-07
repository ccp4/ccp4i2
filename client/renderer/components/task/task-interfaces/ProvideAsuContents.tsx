import {
  Alert,
  Box,
  Grid2,
  LinearProgress,
  Stack,
  Typography,
} from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useJob, valueOfItem } from "../../../utils";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useCallback, useEffect, useRef } from "react";
import { useApi } from "../../../api";
import { BaseSpacegroupCellElement } from "../task-elements/base-spacegroup-cell-element";
import { usePopcorn } from "../../../providers/popcorn-provider";
import {
  chainsFromDigest,
  deduplicateChains,
  mergeAsuEntries,
  type AsuSequenceEntry,
} from "../task-elements/mmcif-sequence-parser";

/**
 * ProvideAsuContents: say what is in the asymmetric unit.
 *
 * The page is built around the table of sequences (ASU_CONTENT), which is what
 * the job writes out. Below it, three file slots fill the table: an existing
 * ASU file, a model, or a sequence file. Each is an ordinary file widget, so
 * each file can be browsed for on disk, picked from the project, or fetched
 * from the web (PDBe / RCSB / AlphaFold for a model, UniProt or a PDB entry
 * for a sequence). Loading a file APPENDS to the table; a sequence already
 * there gains copies rather than a duplicate row; Clear on the table starts
 * over. Last, a reflection file gives the cell for a Matthews estimate of how
 * many copies the crystal holds.
 */

/** The ASU-file digest already lists rows; the other two digests list chains. */
function entriesFromDigest(digest: any, label: string): { entries: AsuSequenceEntry[]; chains: number } {
  if (Array.isArray(digest?.sequences) && !digest?.composition) {
    const entries = digest.sequences.map((seq: any) => ({
      name: seq.name,
      sequence: seq.sequence,
      polymerType: seq.polymerType,
      description: seq.description,
      nCopies: seq.nCopies ?? 1,
    }));
    return { entries, chains: entries.length };
  }
  const chains = chainsFromDigest(digest, label);
  return { entries: deduplicateChains(chains), chains: chains.length };
}

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const api = useApi();
  const { setMessage } = usePopcorn();
  const { useTaskItem, useFileDigest, fetchDigest, getErrors, mutateValidation, mutateContainer } =
    useJob(job.id);
  const { forceUpdate: forceSetAsuContent, item: asuContentItem } = useTaskItem("ASU_CONTENT");
  const { item: asuContentInItem, value: asuContentInValue } = useTaskItem("ASUCONTENTIN");
  const { value: HKLINValue } = useTaskItem("HKLIN");

  // File digest for HKLIN (used for Matthews calculation)
  // Only fetch when a file has been uploaded (has dbFileId) - otherwise digest endpoint fails
  const hasHKLINFile = Boolean(HKLINValue?.dbFileId);
  const hklinDigestPath = hasHKLINFile ? "ProvideAsuContents.inputData.HKLIN" : "";
  const { data: HKLINDigest, mutate: mutateHKLINDigest } = useFileDigest(hklinDigestPath);

  // ASU content is valid when there are no validation errors for it
  const asuContentErrors = getErrors(asuContentItem);
  const rowCount = asuContentItem?._value?.length ?? 0;
  const isAsuContentValid = asuContentErrors.length === 0 && rowCount > 0;

  /** Molecular weight of the table's contents; only asked for when they are valid. */
  const { data: molWeight, mutate: mutateMolWeight } = api.objectMethod<any>(
    job.id,
    "ProvideAsuContents.inputData.ASU_CONTENT",
    "molecularWeight",
    {},
    [],
    isAsuContentValid
  );

  // The hook keeps its last result while switched off, so after Clear it
  // would still report the weight of what was just removed. Only a valid
  // table has a weight.
  const molWeightValue: number | undefined = isAsuContentValid ? molWeight?.data?.result : undefined;

  /** Matthews analysis: needs the molecular weight and the reflection file's cell. */
  const { data: matthewsAnalysis, mutate: mutateMatthews } = api.objectMethod<any>(
    job.id,
    "ProvideAsuContents.inputData.HKLIN.fileContent",
    "matthewsCoeff",
    { molWt: molWeightValue },
    [molWeightValue, HKLINDigest],
    !!(molWeightValue && HKLINDigest)
  );

  const refreshDerived = useCallback(async () => {
    mutateValidation();
    await mutateMolWeight();
    mutateMatthews();
  }, [mutateValidation, mutateMolWeight, mutateMatthews]);

  /**
   * A source file was chosen: read its digest and append its sequences to the
   * table. `updatedItem` comes from onChange and carries the fresh object path;
   * when a pulldown choice has not yet reached the container, that is the only
   * path that names the new file.
   */
  const fillFrom = useCallback(
    async (updatedItem: any, fallbackItem: any) => {
      const objectPath = updatedItem?._objectPath || fallbackItem?._objectPath;
      if (!objectPath) return;
      const chosen = updatedItem ? valueOfItem(updatedItem) : valueOfItem(fallbackItem);
      if (!chosen?.dbFileId && !chosen?.baseName) return; // cleared, not chosen
      const label = chosen?.annotation || chosen?.baseName || "file";

      const digest = await fetchDigest(objectPath);
      if (!digest || digest.status === "Failed") {
        setMessage(`Could not read ${label}${digest?.reason ? `: ${digest.reason}` : ""}`);
        return;
      }
      const { entries, chains } = entriesFromDigest(digest, label);
      if (entries.length === 0) {
        setMessage(`No polymer sequences found in ${label}`);
        return;
      }
      const existing: AsuSequenceEntry[] = Array.isArray(valueOfItem(asuContentItem))
        ? valueOfItem(asuContentItem)
        : [];
      const merged = mergeAsuEntries(existing, entries);
      await forceSetAsuContent(merged);
      // Force a full container refetch so the table re-renders (a function-updater
      // patch with revalidate:false may not reach every subscriber).
      await mutateContainer();
      await refreshDerived();

      const added = merged.length - existing.length;
      const folded = entries.length - added;
      const parts = [`${chains} chain${chains === 1 ? "" : "s"} from ${label}`];
      if (added) parts.push(`${added} new row${added === 1 ? "" : "s"}`);
      if (folded) parts.push(`${folded} added as copies of rows already present`);
      setMessage(parts.join(": "));
    },
    [fetchDigest, asuContentItem, forceSetAsuContent, mutateContainer, refreshDerived, setMessage]
  );

  // When ASUCONTENTIN is already set (e.g. imperatively by another task) and the
  // table is empty, fill the table from it once on mount.
  const hasAutoPopulated = useRef(false);
  useEffect(() => {
    if (hasAutoPopulated.current) return;
    const hasFile = Boolean(asuContentInValue?.dbFileId);
    if (hasFile && rowCount === 0 && asuContentInItem?._objectPath) {
      hasAutoPopulated.current = true;
      fillFrom(asuContentInItem, asuContentInItem);
    }
  }, [asuContentInValue?.dbFileId, rowCount, asuContentInItem, fillFrom]);

  // Why the Matthews panel is empty, when it is.
  const matthewsResults =
    molWeightValue && matthewsAnalysis?.success ? matthewsAnalysis?.data?.result?.results : null;
  let matthewsReason: string | null = null;
  if (!hasHKLINFile) matthewsReason = "Choose a reflection file: its cell sets how many copies fit.";
  else if (rowCount === 0) matthewsReason = "Add at least one sequence to the table first.";
  else if (!isAsuContentValid) matthewsReason = "Fix the sequences flagged in the table first.";
  else if (!molWeightValue || !HKLINDigest) matthewsReason = "Calculating…";
  else if (!matthewsResults) matthewsReason = "No Matthews estimate for this cell and contents.";

  return (
    <CCP4i2Tabs {...props}>
      <CCP4i2Tab label="Main inputs">
        {/* 1. The table: what the job writes out */}
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ guiLabel: "Contents of the asymmetric unit" }}
          containerHint="FolderLevel"
          initiallyOpen={true}
        >
          <CCP4i2TaskElement
            {...props}
            itemName="ASU_CONTENT"
            qualifiers={{ guiLabel: "ASU contents" }}
            onChange={refreshDerived}
          />
          <Stack direction="row" spacing={2} alignItems="center" sx={{ mt: 1 }}>
            <Typography variant="body2">
              Molecular weight:{" "}
              {molWeightValue
                ? `${molWeightValue.toFixed(0)} Da`
                : isAsuContentValid
                  ? "calculating…"
                  : rowCount === 0
                    ? "–"
                    : "not available until the table is valid"}
            </Typography>
            {rowCount > 0 && asuContentErrors.length > 0 && (
              <Typography variant="body2" color="error">
                {asuContentErrors.length} problem{asuContentErrors.length === 1 ? "" : "s"} in the table
              </Typography>
            )}
          </Stack>
        </CCP4i2ContainerElement>

        {/* 2. Sources: each slot appends to the table */}
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ guiLabel: "Fill the table from a file" }}
          containerHint="FolderLevel"
          initiallyOpen={true}
        >
          <Typography variant="caption" color="text.secondary" sx={{ display: "block", mb: 1 }}>
            Each file adds its sequences to the table; one already there gains copies
            instead of a duplicate row. Files can be browsed for, picked from the project,
            or fetched from the web with the buttons beside each slot. For a whole PDB entry use the model slot.
          </Typography>
          <CCP4i2TaskElement
            {...props}
            itemName="ASUCONTENTIN"
            qualifiers={{ guiLabel: "An existing ASU content file", allowCreate: false }}
            onChange={(updated: any) => fillFrom(updated, asuContentInItem)}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="XYZIN"
            qualifiers={{ guiLabel: "A model: every polymer chain is added (fetch a whole PDB entry here, from PDBe, RCSB or AlphaFold)" }}
            onChange={(updated: any) => fillFrom(updated, null)}
          />
          <CCP4i2TaskElement
            {...props}
            itemName="SEQIN"
            qualifiers={{ guiLabel: "One sequence: a file, or fetched from UniProt (a PDB entry gives one chain; use the model slot for all of it)" }}
            onChange={(updated: any) => fillFrom(updated, null)}
          />
        </CCP4i2ContainerElement>

        {/* 3. Cell → how many copies fit */}
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ guiLabel: "Solvent analysis" }}
          containerHint="FolderLevel"
          initiallyOpen={true}
        >
          <Grid2 container spacing={2}>
            <Grid2 size={{ xs: 12, sm: 8 }}>
              <CCP4i2TaskElement
                {...props}
                itemName="HKLIN"
                qualifiers={{ guiLabel: "Reflections (for the Matthews analysis)" }}
                onChange={() => mutateHKLINDigest()}
              />
              {HKLINDigest && (
                <Stack spacing={1} sx={{ mt: 1 }}>
                  <BaseSpacegroupCellElement data={HKLINDigest} />
                  {HKLINDigest.cell && (!HKLINDigest.wavelengths || HKLINDigest.wavelengths.length === 0) && (
                    <Alert severity="info" sx={{ py: 0 }}>
                      No wavelength information in this MTZ file
                    </Alert>
                  )}
                </Stack>
              )}
            </Grid2>
            <Grid2 size={{ xs: 12, sm: 4 }}>
              {matthewsResults ? (
                <Box
                  sx={{
                    p: 1.5,
                    borderRadius: 1,
                    bgcolor: "action.hover",
                    border: 1,
                    borderColor: "divider",
                  }}
                >
                  <Typography variant="caption" color="text.secondary" sx={{ mb: 1, display: "block" }}>
                    Matthews analysis: copies of the table's contents per ASU
                  </Typography>
                  <Stack spacing={1}>
                    {matthewsResults.map(
                      (result: { nmol_in_asu: number; percent_solvent: number; prob_matth: number }) => {
                        const probability = result.prob_matth;
                        const isLikely = probability > 0.5;
                        return (
                          <Box
                            key={result.nmol_in_asu}
                            sx={{
                              p: 1,
                              borderRadius: 0.5,
                              bgcolor: isLikely ? "success.main" : "background.paper",
                              color: isLikely ? "success.contrastText" : "text.primary",
                              border: 1,
                              borderColor: isLikely ? "success.main" : "divider",
                            }}
                          >
                            <Stack direction="row" justifyContent="space-between" alignItems="center">
                              <Typography variant="body2" fontWeight="medium">
                                {result.nmol_in_asu} &times; these contents
                              </Typography>
                              <Typography variant="body2" fontWeight="bold">
                                {(probability * 100).toFixed(0)}%
                              </Typography>
                            </Stack>
                            <Stack direction="row" spacing={2} sx={{ mt: 0.5 }}>
                              <Typography variant="caption" sx={{ opacity: isLikely ? 0.9 : 0.7 }}>
                                {result.percent_solvent.toFixed(1)}% solvent
                              </Typography>
                              <LinearProgress
                                variant="determinate"
                                value={probability * 100}
                                sx={{
                                  flex: 1,
                                  alignSelf: "center",
                                  height: 4,
                                  borderRadius: 2,
                                  bgcolor: isLikely ? "success.light" : "action.disabledBackground",
                                  "& .MuiLinearProgress-bar": {
                                    bgcolor: isLikely ? "success.contrastText" : "primary.main",
                                  },
                                }}
                              />
                            </Stack>
                          </Box>
                        );
                      }
                    )}
                  </Stack>
                </Box>
              ) : (
                <Typography variant="body2" color="text.secondary">
                  {matthewsReason}
                </Typography>
              )}
            </Grid2>
          </Grid2>
        </CCP4i2ContainerElement>
      </CCP4i2Tab>
    </CCP4i2Tabs>
  );
};
export default TaskInterface;
