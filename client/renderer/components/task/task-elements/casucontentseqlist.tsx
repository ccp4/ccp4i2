import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import {
  Autocomplete,
  Box,
  Button,
  Card,
  CardActionArea,
  CardContent,
  Chip,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  IconButton,
  Stack,
  TextField,
  Tooltip,
  Typography,
} from "@mui/material";
import { useApi } from "../../../api";
import { useJob, usePrevious, valueOfItem } from "../../../utils";
import { useCallback, useEffect, useRef, useState } from "react";
import { Add, Delete, Download, Science } from "@mui/icons-material";
import { apiFetch, apiText } from "../../../api-fetch";
import { usePopcorn } from "../../../providers/popcorn-provider";
import {
  parseMultiChainFasta,
  deduplicateChains,
  ChainSequenceInfo,
} from "./mmcif-sequence-parser";

/** Get abbreviated polymer type label */
const getPolymerTypeLabel = (type: string): string => {
  switch (type?.toUpperCase()) {
    case "PROTEIN":
      return "Protein";
    case "DNA":
      return "DNA";
    case "RNA":
      return "RNA";
    default:
      return type || "Unknown";
  }
};

/** Get color for polymer type */
const getPolymerTypeColor = (
  type: string
): "primary" | "secondary" | "success" | "warning" | "info" => {
  switch (type?.toUpperCase()) {
    case "PROTEIN":
      return "primary";
    case "DNA":
      return "info";
    case "RNA":
      return "success";
    default:
      return "warning";
  }
};

/** Format sequence for preview (first N residues with ellipsis) */
const formatSequencePreview = (sequence: string, maxLength = 60): string => {
  if (!sequence) return "No sequence";
  const seq = sequence.replace(/\s/g, "");
  if (seq.length <= maxLength) return seq;
  return seq.substring(0, maxLength) + "…";
};

/** PDB import sources */
const PDB_SOURCES = [
  { value: "ebi", label: "PDBe (EBI)" },
  { value: "rcsb", label: "RCSB PDB" },
] as const;

type PdbSource = (typeof PDB_SOURCES)[number]["value"];

/**
 * Parse PDBe molecules API response into ChainSequenceInfo[].
 */
function parsePdbeMolecules(data: any, pdbId: string): ChainSequenceInfo[] {
  const entities = data[pdbId.toLowerCase()];
  if (!Array.isArray(entities)) return [];
  const chains: ChainSequenceInfo[] = [];
  for (const entity of entities) {
    if (!entity.sequence || entity.molecule_type === "water" || entity.molecule_type === "bound") continue;
    const polymerType = entity.molecule_type?.includes("polypeptide") ? "PROTEIN"
      : entity.molecule_type?.includes("polyribonucleotide") && !entity.molecule_type?.includes("deoxy") ? "RNA"
      : entity.molecule_type?.includes("polydeoxyribonucleotide") ? "DNA"
      : "OTHER" as const;
    const chainIds: string[] = entity.in_chains || [];
    const moleculeName = Array.isArray(entity.molecule_name) ? entity.molecule_name[0] : entity.molecule_name || "";
    for (const chainId of chainIds) {
      chains.push({
        chainId,
        sequence: entity.sequence,
        polymerType,
        length: entity.sequence.length,
        description: moleculeName,
      });
    }
  }
  return chains;
}

export const CAsuContentSeqListElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => {
  const api = useApi();
  const { itemName, job } = props;
  const { setMessage } = usePopcorn();

  // Use separate state for dialog open vs which item is selected
  // This prevents re-renders from closing the dialog
  const [isDialogOpen, setIsDialogOpen] = useState(false);

  // Import from PDB state
  const [importDialogOpen, setImportDialogOpen] = useState(false);
  const [importSource, setImportSource] = useState<PdbSource>("ebi");
  const [importPdbId, setImportPdbId] = useState("");
  const [importInFlight, setImportInFlight] = useState(false);

  // Store the object path when we open, so it doesn't change during editing
  const selectedObjectPathRef = useRef<string | null>(null);

  const {
    useTaskItem,
    setParameter,
    mutateContainer,
    getValidationColor,
  } = useJob(job.id);

  const { item, update: updateList, value: itemValue } = useTaskItem(itemName);
  const previousItemValue = usePrevious(itemValue);
  const previousListLength = usePrevious(item?._value?.length);

  const handleOpenDialog = useCallback(
    (index: number) => {
      if (!item?._value?.[index]) return;

      // Store the object path at the time of opening
      selectedObjectPathRef.current = item._value[index]._objectPath;
      setIsDialogOpen(true);
    },
    [item?._value]
  );

  const handleCloseDialog = useCallback(() => {
    setIsDialogOpen(false);
    selectedObjectPathRef.current = null;
    mutateContainer();
    // Notify parent that content may have changed (triggers validity/molWeight recalc)
    props.onChange?.(item);
  }, [mutateContainer, props.onChange, item]);

  const extendListItem = useCallback(async () => {
    if (!updateList) return;
    var taskElement = JSON.parse(JSON.stringify(item._subItem));
    taskElement._objectPath = taskElement._objectPath.replace(
      "[?]",
      "[" + item._value.length + "]"
    );
    for (var valueElementKey in taskElement._value) {
      var valueElement = taskElement._value[valueElementKey];
      valueElement._objectPath = valueElement._objectPath.replace(
        "[?]",
        "[" + item._value.length + "]"
      );
    }
    const listValue = Array.isArray(valueOfItem(item)) ? valueOfItem(item) : [];
    let newItemValue = valueOfItem(taskElement);
    listValue.push(newItemValue);
    await updateList(listValue);
  }, [item, updateList]);

  // When the list grows by exactly one item (user clicked "Add"), open the new item.
  // Skip when growing by more than one (bulk load from file) to avoid unwanted dialog.
  useEffect(() => {
    const currentLength = item?._value?.length ?? 0;
    if (
      previousListLength !== undefined &&
      currentLength === previousListLength + 1
    ) {
      const newIndex = currentLength - 1;
      handleOpenDialog(newIndex);
    }
  }, [item?._value?.length, previousListLength, handleOpenDialog]);

  const deleteItem = useCallback(
    async (index: number) => {
      const array = item._value;
      if (index > -1 && index < array.length) {
        array.splice(index, 1);
        const setParameterArg = {
          object_path: item._objectPath,
          value: valueOfItem(item),
        };
        const updateResult: any = await setParameter(setParameterArg);
        if (props.onChange) {
          await props.onChange(updateResult.updated_item);
        }
      }
    },
    [item, setParameter, props.onChange]
  );

  /**
   * Import all polymer chains from a PDB entry as ASU content sequences.
   * Fetches chain data from PDBe or RCSB, then bulk-populates the sequence list.
   */
  const handleImportFromPdb = useCallback(async () => {
    if (!importPdbId.trim() || !updateList) return;
    const id = importPdbId.trim().toLowerCase();
    setImportInFlight(true);
    try {
      setMessage(`Fetching composition for ${id.toUpperCase()}...`);
      let chains: ChainSequenceInfo[];

      if (importSource === "ebi") {
        const url = `https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/${id}`;
        const result = await apiFetch(url);
        const data = await result.json();
        chains = parsePdbeMolecules(data, id);
      } else {
        const fastaUrl = `https://www.rcsb.org/fasta/entry/${id.toUpperCase()}`;
        const fastaText = await apiText(fastaUrl);
        chains = parseMultiChainFasta(fastaText);
      }

      if (chains.length === 0) {
        setMessage(`No polymer chains found in ${id.toUpperCase()}`);
        return;
      }

      // Deduplicate identical chains (e.g. CDK2 chains A,C → nCopies=2)
      const seqList = deduplicateChains(chains);

      await updateList(seqList);
      await mutateContainer();
      if (props.onChange) {
        props.onChange(item);
      }
      const uniqueCount = seqList.length;
      const totalChains = chains.length;
      const msg = uniqueCount === totalChains
        ? `Imported ${totalChains} chain(s) from ${id.toUpperCase()}`
        : `Imported ${totalChains} chains as ${uniqueCount} unique sequence(s) from ${id.toUpperCase()}`;
      setMessage(msg);
      setImportDialogOpen(false);
      setImportPdbId("");
    } catch (err: any) {
      setMessage(err.message || `Failed to fetch ${id.toUpperCase()}`);
    } finally {
      setImportInFlight(false);
    }
  }, [importPdbId, importSource, updateList, mutateContainer, item, props.onChange, setMessage]);

  useEffect(() => {
    console.debug("CAsuContentSeqListElement mounted");
    return () => {
      console.debug("CAsuContentSeqListElement unmounted");
    };
  }, []);

  const hasSequences = (item?._value?.length ?? 0) > 0;

  return (
    item && (
      <>
        {/* Header with add/import buttons - only show when there are sequences */}
        {hasSequences && (
          <Stack
            direction="row"
            justifyContent="space-between"
            alignItems="center"
            sx={{ mb: 2 }}
          >
            <Typography variant="body2" color="text.secondary">
              Click a card to edit sequence details
            </Typography>
            <Stack direction="row" spacing={1}>
              {job.status === 1 && (
                <Button
                  variant="outlined"
                  size="small"
                  startIcon={<Download />}
                  onClick={() => setImportDialogOpen(true)}
                >
                  Import from PDB
                </Button>
              )}
              <Button
                variant="contained"
                size="small"
                startIcon={<Add />}
                onClick={async () => {
                  await extendListItem();
                }}
              >
                Add Sequence
              </Button>
            </Stack>
          </Stack>
        )}

        {/* Sequence cards */}
        <Stack spacing={1.5}>
          {item?._value?.map((contentElement: any, iElement: number) => {
            const rowValidationColor = getValidationColor(contentElement);
            const name = contentElement._value?.name?._value || `Sequence ${iElement + 1}`;
            const polymerType = contentElement._value?.polymerType?._value || "";
            const description = contentElement._value?.description?._value || "";
            const nCopies = contentElement._value?.nCopies?._value ?? 1;
            const sequence = contentElement._value?.sequence?._value || "";
            const seqLength = sequence.replace(/\s/g, "").length;

            return (
              <Card
                key={`${iElement}`}
                variant="outlined"
                sx={{
                  borderLeft:
                    rowValidationColor !== "inherit"
                      ? `4px solid ${rowValidationColor}`
                      : undefined,
                  transition: "all 0.2s",
                  "&:hover": {
                    boxShadow: 2,
                    borderColor: "primary.main",
                  },
                }}
              >
                <CardActionArea onClick={() => handleOpenDialog(iElement)}>
                  <CardContent sx={{ py: 1.5, "&:last-child": { pb: 1.5 } }}>
                    <Stack spacing={1}>
                      {/* Top row: name, type, copies, delete */}
                      <Stack
                        direction="row"
                        alignItems="center"
                        justifyContent="space-between"
                      >
                        <Stack direction="row" alignItems="center" spacing={1.5}>
                          <Science fontSize="small" color="action" />
                          <Typography variant="subtitle2" fontWeight="bold">
                            {name}
                          </Typography>
                          <Chip
                            label={getPolymerTypeLabel(polymerType)}
                            size="small"
                            color={getPolymerTypeColor(polymerType)}
                            variant="outlined"
                          />
                          <Chip
                            label={`×${nCopies} ${nCopies === 1 ? "copy" : "copies"}`}
                            size="small"
                            variant="filled"
                            color="default"
                          />
                        </Stack>
                        <Stack direction="row" alignItems="center" spacing={1}>
                          <Typography variant="caption" color="text.secondary">
                            {seqLength} residues
                          </Typography>
                          <Tooltip title="Delete sequence">
                            <IconButton
                              size="small"
                              color="error"
                              onClick={(ev) => {
                                ev.stopPropagation();
                                ev.preventDefault();
                                deleteItem(iElement);
                              }}
                            >
                              <Delete fontSize="small" />
                            </IconButton>
                          </Tooltip>
                        </Stack>
                      </Stack>

                      {/* Description if present */}
                      {description && (
                        <Typography
                          variant="body2"
                          color="text.secondary"
                          sx={{
                            overflow: "hidden",
                            textOverflow: "ellipsis",
                            whiteSpace: "nowrap",
                          }}
                        >
                          {description}
                        </Typography>
                      )}

                      {/* Sequence preview */}
                      <Box
                        sx={{
                          p: 1,
                          bgcolor: "action.hover",
                          borderRadius: 0.5,
                          fontFamily: "monospace",
                          fontSize: "0.75rem",
                          letterSpacing: "0.05em",
                          color: "text.secondary",
                          overflow: "hidden",
                          textOverflow: "ellipsis",
                          whiteSpace: "nowrap",
                        }}
                      >
                        {formatSequencePreview(sequence)}
                      </Box>
                    </Stack>
                  </CardContent>
                </CardActionArea>
              </Card>
            );
          })}
        </Stack>

        {/* Empty state */}
        {!hasSequences && (
          <Box
            sx={{
              mt: 2,
              p: 3,
              textAlign: "center",
              borderRadius: 2,
              bgcolor: "action.hover",
              border: 1,
              borderColor: "divider",
              borderStyle: "dashed",
            }}
          >
            <Science sx={{ fontSize: 40, color: "action.disabled", mb: 1 }} />
            <Typography variant="body2" color="text.secondary" gutterBottom>
              No sequences defined yet
            </Typography>
            <Stack direction="row" spacing={1} justifyContent="center" sx={{ mt: 1 }}>
              <Button
                variant="outlined"
                size="small"
                startIcon={<Download />}
                onClick={() => setImportDialogOpen(true)}
                disabled={job.status !== 1}
              >
                Import from PDB
              </Button>
              <Button
                variant="outlined"
                size="small"
                startIcon={<Add />}
                onClick={async () => {
                  await extendListItem();
                }}
              >
                Add Manually
              </Button>
            </Stack>
          </Box>
        )}

        {/* Edit dialog */}
        <Dialog
          open={isDialogOpen}
          onClose={handleCloseDialog}
          fullWidth
          maxWidth={false}
          slotProps={{
            paper: {
              style: {
                margin: "1rem",
                width: "calc(100% - 2rem)",
              },
            },
          }}
        >
          <DialogContent>
            {/* Use the stored object path ref so it doesn't change during editing */}
            {isDialogOpen && selectedObjectPathRef.current && (
              <CCP4i2TaskElement
                {...props}
                itemName={selectedObjectPathRef.current}
              />
            )}
          </DialogContent>
          <DialogActions>
            <Button onClick={handleCloseDialog} variant="contained">
              OK
            </Button>
          </DialogActions>
        </Dialog>

        {/* Import from PDB dialog */}
        <Dialog
          open={importDialogOpen}
          onClose={() => setImportDialogOpen(false)}
          maxWidth="sm"
          fullWidth
        >
          <DialogTitle>Import sequences from PDB entry</DialogTitle>
          <DialogContent>
            <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
              Fetch all polymer chains from a PDB entry and add them as ASU content sequences.
              This will replace any existing sequences.
            </Typography>
            <Stack spacing={2} sx={{ mt: 1 }}>
              <Autocomplete
                value={PDB_SOURCES.find((s) => s.value === importSource) || PDB_SOURCES[0]}
                onChange={(_, newValue) => {
                  if (newValue) setImportSource(newValue.value);
                }}
                options={[...PDB_SOURCES]}
                getOptionLabel={(option) => option.label}
                renderInput={(params) => (
                  <TextField {...params} label="Source" size="small" />
                )}
                disableClearable
              />
              <TextField
                label="PDB ID"
                placeholder="e.g. 1jst"
                value={importPdbId}
                onChange={(e) => setImportPdbId(e.target.value)}
                size="small"
                onKeyDown={(e) => {
                  if (e.key === "Enter" && importPdbId.trim() && !importInFlight) {
                    handleImportFromPdb();
                  }
                }}
              />
            </Stack>
          </DialogContent>
          <DialogActions>
            <Button onClick={() => setImportDialogOpen(false)} disabled={importInFlight}>
              Cancel
            </Button>
            <Button
              onClick={handleImportFromPdb}
              disabled={!importPdbId.trim() || importInFlight}
              variant="contained"
            >
              {importInFlight ? "Importing..." : "Import"}
            </Button>
          </DialogActions>
        </Dialog>
      </>
    )
  );
};
