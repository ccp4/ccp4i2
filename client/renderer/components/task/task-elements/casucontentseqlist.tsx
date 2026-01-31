import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import {
  Box,
  Button,
  Card,
  CardActionArea,
  CardContent,
  Chip,
  Dialog,
  DialogActions,
  DialogContent,
  IconButton,
  Stack,
  Tooltip,
  Typography,
} from "@mui/material";
import { useApi } from "../../../api";
import { useJob, usePrevious, valueOfItem } from "../../../utils";
import { useCallback, useEffect, useRef, useState } from "react";
import { Add, Delete, Science } from "@mui/icons-material";

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

export const CAsuContentSeqListElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => {
  const api = useApi();
  const { itemName, job } = props;

  // Use separate state for dialog open vs which item is selected
  // This prevents re-renders from closing the dialog
  const [isDialogOpen, setIsDialogOpen] = useState(false);

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

  // When the list grows, open the newly added item
  useEffect(() => {
    const currentLength = item?._value?.length ?? 0;
    if (
      previousListLength !== undefined &&
      currentLength > previousListLength
    ) {
      // Open the newly added item (last in list)
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
        {/* Header with add button - only show when there are sequences */}
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
            <Button
              variant="outlined"
              size="small"
              startIcon={<Add />}
              onClick={async () => {
                await extendListItem();
              }}
              sx={{ mt: 1 }}
            >
              Add First Sequence
            </Button>
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
      </>
    )
  );
};
