import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import {
  Box,
  Button,
  Chip,
  Dialog,
  DialogActions,
  DialogContent,
  IconButton,
  Stack,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  Tooltip,
  Typography,
} from "@mui/material";
import { useJob, usePrevious, valueOfItem } from "../../../utils";
import { useCallback, useEffect, useRef, useState } from "react";
import { Add, Delete, DeleteSweep, Science } from "@mui/icons-material";

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
const formatSequencePreview = (sequence: string, maxLength = 40): string => {
  if (!sequence) return "";
  const seq = sequence.replace(/\s/g, "");
  if (seq.length <= maxLength) return seq;
  return seq.substring(0, maxLength) + "…";
};

/**
 * The asymmetric unit's contents as a table: one row per distinct sequence,
 * with its polymer type, length and number of copies. Rows are added and
 * removed here and edited in a detail dialog (CAsuContentSeqElement); the
 * table is filled in bulk by the page that owns it, from an ASU, coordinate
 * or sequence file, so it offers no file or web fetching of its own.
 */
export const CAsuContentSeqListElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => {
  const { itemName, job } = props;
  const editable = job.status === 1;

  // Open state and the row's object path are kept apart: the path lives in a
  // ref so a re-render while editing does not close the dialog.
  const [isDialogOpen, setIsDialogOpen] = useState(false);
  const selectedObjectPathRef = useRef<string | null>(null);

  const { useTaskItem, setParameter, mutateContainer, getValidationColor } =
    useJob(job.id);

  const { item, update: updateList } = useTaskItem(itemName);
  const previousListLength = usePrevious(item?._value?.length);

  const openRow = useCallback(
    (index: number) => {
      if (!item?._value?.[index]) return;
      selectedObjectPathRef.current = item._value[index]._objectPath;
      setIsDialogOpen(true);
    },
    [item?._value]
  );

  const closeRow = useCallback(() => {
    setIsDialogOpen(false);
    selectedObjectPathRef.current = null;
    mutateContainer();
    // Notify parent that content may have changed (triggers validity/molWeight recalc)
    props.onChange?.(item);
  }, [mutateContainer, props.onChange, item]);

  const addRow = useCallback(async () => {
    if (!updateList) return;
    const taskElement = JSON.parse(JSON.stringify(item._subItem));
    taskElement._objectPath = taskElement._objectPath.replace(
      "[?]",
      "[" + item._value.length + "]"
    );
    for (const key in taskElement._value) {
      const valueElement = taskElement._value[key];
      valueElement._objectPath = valueElement._objectPath.replace(
        "[?]",
        "[" + item._value.length + "]"
      );
    }
    const listValue = Array.isArray(valueOfItem(item)) ? valueOfItem(item) : [];
    listValue.push(valueOfItem(taskElement));
    await updateList(listValue);
  }, [item, updateList]);

  // When the list grows by exactly one row (Add), open it for editing. A bulk
  // fill from a file grows it by more, and opening one row then is a surprise.
  useEffect(() => {
    const currentLength = item?._value?.length ?? 0;
    if (previousListLength !== undefined && currentLength === previousListLength + 1) {
      openRow(currentLength - 1);
    }
  }, [item?._value?.length, previousListLength, openRow]);

  const deleteRow = useCallback(
    async (index: number) => {
      const array = item._value;
      if (index > -1 && index < array.length) {
        array.splice(index, 1);
        const updateResult: any = await setParameter({
          object_path: item._objectPath,
          value: valueOfItem(item),
        });
        if (props.onChange) {
          await props.onChange(updateResult.updated_item);
        }
      }
    },
    [item, setParameter, props.onChange]
  );

  const clearRows = useCallback(async () => {
    if (!updateList) return;
    await updateList([]);
    await mutateContainer();
    props.onChange?.(item);
  }, [updateList, mutateContainer, props.onChange, item]);

  if (!item) return null;
  const rows: any[] = item._value ?? [];
  const hasRows = rows.length > 0;
  const dialogOpen = isDialogOpen && Boolean(selectedObjectPathRef.current);

  return (
    <>
      <Stack direction="row" justifyContent="space-between" alignItems="center" sx={{ mb: 1 }}>
        <Typography variant="body2" color="text.secondary">
          {hasRows
            ? "Click a row to edit it."
            : "Add a sequence by hand, or fill the table from a file below."}
        </Typography>
        <Stack direction="row" spacing={1}>
          {hasRows && (
            <Tooltip title="Remove every row">
              <span>
                <Button
                  variant="outlined"
                  size="small"
                  color="inherit"
                  startIcon={<DeleteSweep />}
                  onClick={clearRows}
                  disabled={!editable}
                >
                  Clear
                </Button>
              </span>
            </Tooltip>
          )}
          <Button
            variant="contained"
            size="small"
            startIcon={<Add />}
            onClick={addRow}
            disabled={!editable}
          >
            Add sequence
          </Button>
        </Stack>
      </Stack>

      {hasRows ? (
        <Table size="small" sx={{ "& td, & th": { whiteSpace: "nowrap" } }}>
          <TableHead>
            <TableRow>
              <TableCell>Name</TableCell>
              <TableCell>Type</TableCell>
              <TableCell align="right">Residues</TableCell>
              <TableCell align="right">Copies</TableCell>
              <TableCell sx={{ width: "100%" }}>Sequence</TableCell>
              <TableCell padding="checkbox" />
            </TableRow>
          </TableHead>
          <TableBody>
            {rows.map((row: any, index: number) => {
              const color = getValidationColor(row);
              const name = row._value?.name?._value || `Sequence ${index + 1}`;
              const polymerType = row._value?.polymerType?._value || "";
              const description = row._value?.description?._value || "";
              const nCopies = row._value?.nCopies?._value ?? 1;
              const sequence = row._value?.sequence?._value || "";
              const seqLength = sequence.replace(/\s/g, "").length;
              return (
                <TableRow
                  key={index}
                  hover
                  onClick={() => openRow(index)}
                  sx={{
                    cursor: "pointer",
                    "& td:first-of-type": {
                      borderLeft: color !== "inherit" ? `4px solid ${color}` : undefined,
                    },
                  }}
                >
                  <TableCell>
                    <Tooltip title={description || ""} placement="top-start">
                      <Typography variant="body2" fontWeight="medium">
                        {name}
                      </Typography>
                    </Tooltip>
                  </TableCell>
                  <TableCell>
                    <Chip
                      label={getPolymerTypeLabel(polymerType)}
                      size="small"
                      color={getPolymerTypeColor(polymerType)}
                      variant="outlined"
                    />
                  </TableCell>
                  <TableCell align="right">{seqLength || "–"}</TableCell>
                  <TableCell align="right">{nCopies}</TableCell>
                  <TableCell
                    sx={{
                      fontFamily: "monospace",
                      fontSize: "0.75rem",
                      color: "text.secondary",
                      maxWidth: 0,
                      overflow: "hidden",
                      textOverflow: "ellipsis",
                    }}
                  >
                    {formatSequencePreview(sequence) || "(no sequence)"}
                  </TableCell>
                  <TableCell padding="checkbox">
                    <Tooltip title="Remove this sequence">
                      <span>
                        <IconButton
                          size="small"
                          color="error"
                          disabled={!editable}
                          onClick={(ev) => {
                            ev.stopPropagation();
                            deleteRow(index);
                          }}
                        >
                          <Delete fontSize="small" />
                        </IconButton>
                      </span>
                    </Tooltip>
                  </TableCell>
                </TableRow>
              );
            })}
          </TableBody>
        </Table>
      ) : (
        <Box
          sx={{
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
          <Typography variant="body2" color="text.secondary">
            Nothing in the asymmetric unit yet
          </Typography>
        </Box>
      )}

      {/* Edit dialog */}
      <Dialog
        open={dialogOpen}
        onClose={closeRow}
        fullWidth
        maxWidth={false}
        slotProps={{
          paper: { style: { margin: "1rem", width: "calc(100% - 2rem)" } },
        }}
      >
        <DialogContent>
          {/* Use the stored object path ref so it doesn't change during editing */}
          {dialogOpen && selectedObjectPathRef.current && (
            <CCP4i2TaskElement {...props} itemName={selectedObjectPathRef.current} />
          )}
        </DialogContent>
        <DialogActions>
          <Button onClick={closeRow} variant="contained">
            OK
          </Button>
        </DialogActions>
      </Dialog>
    </>
  );
};
