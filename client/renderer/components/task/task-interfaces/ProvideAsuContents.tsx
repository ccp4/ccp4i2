import {
  Alert,
  Box,
  Chip,
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
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useApi } from "../../../api";
import { apiPost } from "../../../api-fetch";
import { BaseSpacegroupCellElement } from "../task-elements/base-spacegroup-cell-element";
import { usePopcorn } from "../../../providers/popcorn-provider";
import {
  MatthewsRow,
  bestMatthewsRow,
  commonFactor,
  copiesOfRows,
  describeCopies,
  namesOfRows,
} from "../../../lib/asu-copies";

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
 * over. The server does the appending when the file is set
 * (ProvideAsuContents.afterParameterSet), so a job set up any other way (i2run,
 * an agent) fills the same way; this page shows the table it then holds.
 *
 * Last, a reflection file gives the cell for a Matthews estimate of how many
 * copies the crystal holds. The estimate is of multiples of the table's
 * smallest whole set (the copies divided by their common factor), so the
 * rows stay put while the copies change (#673, #687). Clicking a row sets
 * every entry's copies to that multiple (the server's scaleCopies, so an
 * agent can do the same). When a sequence file fills the table the most
 * probable row is applied by default, until the user picks one: a sequence
 * says nothing about copies. A model does (its chains), and an AU-contents
 * file is the user's own answer, so both are taken as they are (Martin,
 * 2026-10-31: 1H1S's 2 + 2 from the model would have become 3 + 3).
 */

const MatthewsPicker: React.FC<{
  results: MatthewsRow[];
  /** The table's entry names and copies, so each row can say what it would set. */
  names: string[];
  copies: number[];
  factor: number;
  current: number;
  best: MatthewsRow | null;
  editable: boolean;
  onPick: (multiple: number) => void;
}> = ({ results, names, copies, factor, current, best, editable, onPick }) => (
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
      Matthews analysis: what the asymmetric unit could hold
      {editable ? ". Click a row to set the copies to it." : "."}
    </Typography>
    <Stack spacing={1}>
      {results.map((result) => {
        const probability = result.prob_matth;
        const isCurrent = result.nmol_in_asu === current;
        const isBest = best?.nmol_in_asu === result.nmol_in_asu;
        const clickable = editable && !isCurrent;
        return (
          <Box
            key={result.nmol_in_asu}
            role={clickable ? "button" : undefined}
            tabIndex={clickable ? 0 : undefined}
            aria-label={clickable ? `Set the copies to ${result.nmol_in_asu} times` : undefined}
            onClick={clickable ? () => onPick(result.nmol_in_asu) : undefined}
            onKeyDown={
              clickable
                ? (ev) => {
                    if (ev.key === "Enter" || ev.key === " ") {
                      ev.preventDefault();
                      onPick(result.nmol_in_asu);
                    }
                  }
                : undefined
            }
            sx={{
              p: 1,
              borderRadius: 0.5,
              bgcolor: isCurrent ? "success.main" : "background.paper",
              color: isCurrent ? "success.contrastText" : "text.primary",
              border: isBest && !isCurrent ? 2 : 1,
              borderColor: isCurrent ? "success.main" : isBest ? "success.main" : "divider",
              cursor: clickable ? "pointer" : "default",
              "&:hover": clickable ? { bgcolor: "action.selected" } : undefined,
              "&:focus-visible": clickable
                ? { outline: "2px solid", outlineColor: "primary.main" }
                : undefined,
            }}
          >
            {/* The heading is what the table would hold, spelled out per
                entry: a bare multiplier ("2 x") read as "twice the current
                table" when the chip beside it said "current". */}
            <Stack direction="row" justifyContent="space-between" alignItems="flex-start" spacing={1}>
              <Stack direction="row" alignItems="center" spacing={1} flexWrap="wrap" useFlexGap>
                <Typography variant="body2" fontWeight="medium" sx={{ wordBreak: "break-word" }}>
                  {describeCopies(names, copies, result.nmol_in_asu)}
                </Typography>
                {isCurrent && (
                  <Chip label="current" size="small" color="success" variant="outlined"
                    sx={{ color: "inherit", borderColor: "currentColor" }} />
                )}
                {isBest && !isCurrent && (
                  <Chip label="best guess" size="small" color="success" variant="outlined" />
                )}
              </Stack>
              <Typography variant="body2" fontWeight="bold" sx={{ flexShrink: 0 }}>
                {(probability * 100).toFixed(0)}%
              </Typography>
            </Stack>
            <Stack direction="row" spacing={2} sx={{ mt: 0.5 }}>
              <Typography variant="caption" sx={{ opacity: isCurrent ? 0.9 : 0.7, whiteSpace: "nowrap" }}>
                {result.nmol_in_asu} &times; the set &middot; {result.percent_solvent.toFixed(1)}% solvent
              </Typography>
              <LinearProgress
                variant="determinate"
                value={probability * 100}
                sx={{
                  flex: 1,
                  alignSelf: "center",
                  height: 4,
                  borderRadius: 2,
                  bgcolor: isCurrent ? "success.light" : "action.disabledBackground",
                  "& .MuiLinearProgress-bar": {
                    bgcolor: isCurrent ? "success.contrastText" : "primary.main",
                  },
                }}
              />
            </Stack>
          </Box>
        );
      })}
    </Stack>
  </Box>
);

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const api = useApi();
  const { setMessage } = usePopcorn();
  const { useTaskItem, useFileDigest, getErrors, mutateValidation, mutateContainer } =
    useJob(job.id);
  const { item: asuContentItem } = useTaskItem("ASU_CONTENT");
  const { item: asuContentInItem, value: asuContentInValue } = useTaskItem("ASUCONTENTIN");
  const { value: HKLINValue } = useTaskItem("HKLIN");
  const editable = job.status === 1;

  // File digest for HKLIN (used for Matthews calculation)
  // Only fetch when a file has been uploaded (has dbFileId) - otherwise digest endpoint fails
  const hasHKLINFile = Boolean(HKLINValue?.dbFileId);
  const hklinDigestPath = hasHKLINFile ? "ProvideAsuContents.inputData.HKLIN" : "";
  const { data: HKLINDigest, mutate: mutateHKLINDigest } = useFileDigest(hklinDigestPath);

  // ASU content is valid when there are no validation errors for it
  const asuContentErrors = getErrors(asuContentItem);
  const rows: any[] = asuContentItem?._value ?? [];
  const rowCount = rows.length;
  const isAsuContentValid = asuContentErrors.length === 0 && rowCount > 0;

  // The copies as they stand, and their common factor: Matthews is run on
  // the contents divided by it, so the table's rows do not move when the
  // copies change, and the row equal to the factor is the current contents.
  const copies = useMemo(() => copiesOfRows(rows), [rows]);
  const names = useMemo(() => namesOfRows(rows), [rows]);
  const factor = useMemo(() => commonFactor(copies), [copies]);

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
  const reducedMolWt: number | undefined = molWeightValue ? molWeightValue / factor : undefined;

  /** Matthews analysis: needs the molecular weight and the reflection file's cell. */
  const { data: matthewsAnalysis, mutate: mutateMatthews } = api.objectMethod<any>(
    job.id,
    "ProvideAsuContents.inputData.HKLIN.fileContent",
    "matthewsCoeff",
    { molWt: reducedMolWt },
    [reducedMolWt, HKLINDigest],
    !!(reducedMolWt && HKLINDigest)
  );

  const refreshDerived = useCallback(async () => {
    mutateValidation();
    await mutateMolWeight();
    mutateMatthews();
  }, [mutateValidation, mutateMolWeight, mutateMatthews]);

  // Why the Matthews panel is empty, when it is.
  const matthewsResults: MatthewsRow[] | null =
    reducedMolWt && matthewsAnalysis?.success ? matthewsAnalysis?.data?.result?.results : null;
  const best = useMemo(() => bestMatthewsRow(matthewsResults), [matthewsResults]);
  let matthewsReason: string | null = null;
  if (!hasHKLINFile) matthewsReason = "Choose a reflection file: its cell sets how many copies fit.";
  else if (rowCount === 0) matthewsReason = "Add at least one sequence to the table first.";
  else if (!isAsuContentValid) matthewsReason = "Fix the sequences flagged in the table first.";
  else if (!reducedMolWt || !HKLINDigest) matthewsReason = "Calculating…";
  else if (!matthewsResults) matthewsReason = "No Matthews estimate for this cell and contents.";

  /** Set every entry's copies to `multiple` of the smallest whole set. */
  const applyingRef = useRef(false);
  const pickMultiple = useCallback(
    async (multiple: number, how: "chosen" | "default") => {
      if (applyingRef.current) return;
      applyingRef.current = true;
      try {
        const response: any = await apiPost(`jobs/${job.id}/object_method/`, {
          object_path: "ProvideAsuContents",
          method_name: "scaleCopies",
          args: [multiple],
          kwargs: {},
        });
        const result = response?.data?.result;
        if (!result?.success) {
          setMessage(`Could not set the copies: ${result?.error || response?.error || "unknown error"}`);
          return;
        }
        await mutateContainer();
        await refreshDerived();
        const row = matthewsResults?.find((r) => r.nmol_in_asu === multiple);
        const detail = row ? ` (${(row.prob_matth * 100).toFixed(0)}%, ${row.percent_solvent.toFixed(0)}% solvent)` : "";
        const what = describeCopies(names, copies, multiple);
        setMessage(
          how === "default"
            ? `Copies set to the Matthews best guess: ${what}${detail}. Click another row to change them.`
            : `Copies set to ${what}${detail}`
        );
      } catch (err: any) {
        setMessage(`Could not set the copies: ${err?.message || err}`);
      } finally {
        applyingRef.current = false;
      }
    },
    [job.id, mutateContainer, refreshDerived, setMessage, matthewsResults, names, copies]
  );

  // The default is the best guess, after a sequence file fills the table
  // only: the most probable multiple is applied once the analysis is in,
  // unless the user has picked a row. A fresh token per event, so stale
  // results from before it cannot satisfy it.
  const userPickedRef = useRef(false);
  const [autoApplyToken, setAutoApplyToken] = useState(0);
  const handledTokenRef = useRef(0);
  const requestDefault = useCallback(() => {
    if (editable && !userPickedRef.current) setAutoApplyToken((t) => t + 1);
  }, [editable]);
  useEffect(() => {
    if (autoApplyToken === 0 || handledTokenRef.current === autoApplyToken) return;
    if (!matthewsResults || !best) return;
    handledTokenRef.current = autoApplyToken;
    if (best.nmol_in_asu !== factor) pickMultiple(best.nmol_in_asu, "default");
  }, [autoApplyToken, matthewsResults, best, factor, pickMultiple]);

  const onPick = useCallback(
    (multiple: number) => {
      userPickedRef.current = true;
      pickMultiple(multiple, "chosen");
    },
    [pickMultiple]
  );

  /**
   * A source file was chosen. The server added its sequences to the table when
   * the file was set; show the table it now holds, and what was added.
   */
  const fillFrom = useCallback(
    async (updatedItem: any, fallbackItem: any, note?: string | null, thenDefault = false) => {
      const chosen = updatedItem ? valueOfItem(updatedItem) : valueOfItem(fallbackItem);
      if (!chosen?.dbFileId && !chosen?.baseName) return; // cleared, not chosen
      const label = chosen?.annotation || chosen?.baseName || "file";
      await mutateContainer();
      await refreshDerived();
      setMessage(note || `Sequences from ${label} added to the table`);
      if (thenDefault) requestDefault();
    },
    [mutateContainer, refreshDerived, setMessage, requestDefault]
  );

  // When ASUCONTENTIN is already set (e.g. imperatively by another task) and the
  // table is empty, fill the table from it once on mount.
  const hasAutoPopulated = useRef(false);
  useEffect(() => {
    if (hasAutoPopulated.current) return;
    const hasFile = Boolean(asuContentInValue?.dbFileId);
    if (hasFile && rowCount === 0 && asuContentInItem?._objectPath) {
      hasAutoPopulated.current = true;
      // The file reached this job without being set here (carried over from
      // an earlier job), so the server has not filled the table from it yet
      apiPost(`jobs/${job.id}/object_method/`, {
        object_path: "ProvideAsuContents",
        method_name: "fillFromFile",
        args: ["ASUCONTENTIN"],
        kwargs: {},
      })
        .then((response: any) =>
          fillFrom(asuContentInItem, asuContentInItem, response?.data?.result?.note))
        .catch((err: any) => setMessage(`Could not fill the table: ${err?.message || err}`));
    }
  }, [asuContentInValue?.dbFileId, rowCount, asuContentInItem, fillFrom, job.id, setMessage]);

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
              <Stack spacing={0}>
                {asuContentErrors.map((e) => {
                  // "...ASU_CONTENT[1].sequence" -> the entry's name
                  const index = Number(/\[(\d+)\]/.exec(e.path)?.[1]);
                  const name = Number.isFinite(index)
                    ? asuContentItem?._value?.[index]?._value?.name?._value || `Entry ${index + 1}`
                    : null;
                  const text = (e.error?.messages ?? []).join("; ") || "problem";
                  return (
                    <Typography key={e.path} variant="body2" color="error">
                      {name ? `${name}: ${text}` : text}
                    </Typography>
                  );
                })}
              </Stack>
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
            onChange={(updated: any) => fillFrom(updated, null, undefined, true)}
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
                <MatthewsPicker
                  results={matthewsResults}
                  names={names}
                  copies={copies}
                  factor={factor}
                  current={factor}
                  best={best}
                  editable={editable}
                  onPick={onPick}
                />
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
