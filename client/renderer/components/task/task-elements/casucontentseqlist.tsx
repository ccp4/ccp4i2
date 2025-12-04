import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  Toolbar,
  Typography,
} from "@mui/material";
import { useApi } from "../../../api";
import { useJob, usePrevious, valueOfItem } from "../../../utils";
import { useCallback, useEffect, useRef, useState } from "react";
import { Add, Delete } from "@mui/icons-material";
import { useParameterChangeIntent } from "../../../providers/parameter-change-intent-provider";

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

  const { useTaskItem, setParameter, container, mutateContainer, getValidationColor } = useJob(
    job.id
  );
  const { intent, setIntent, clearIntent } = useParameterChangeIntent();

  const { item, update: updateList, value: itemValue } = useTaskItem(itemName);
  const previousItemValue = usePrevious(itemValue);

  const handleOpenDialog = useCallback((index: number) => {
    if (!item?._value?.[index]) return;

    // Store the object path at the time of opening
    selectedObjectPathRef.current = item._value[index]._objectPath;
    setIsDialogOpen(true);
  }, [item?._value]);

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
  }, [item, job, updateList, setIntent, mutateContainer]);

  // After reload, select the new item if intent matches
  useEffect(() => {
    if (
      item?._value &&
      intent &&
      typeof intent === "object" &&
      "jobId" in intent &&
      job.id === (intent as any).jobId &&
      "reason" in intent &&
      intent.reason === "UserEdit" &&
      "parameterPath" in intent &&
      intent.parameterPath === item._objectPath &&
      "previousValue" in intent &&
      intent.previousValue !== undefined &&
      JSON.stringify(intent.previousValue) !== JSON.stringify(itemValue)
    ) {
      if (intent.previousValue.length < itemValue.length) {
        // Open the newly added item (last in list)
        const newIndex = item._value.length - 1;
        handleOpenDialog(newIndex);
      }
      clearIntent();
    }
  }, [item?._value, intent, job.id, clearIntent, itemValue, handleOpenDialog]);

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

  return (
    item && (
      <>
        <Toolbar sx={{ m: 0, p: 0, mt: 0, pt: 0 }}>
          <Typography variant="body1" sx={{ mt: 0, pt: 0, flexGrow: 1 }}>
            Click a table row to edit the constituents of the ASU
          </Typography>
          <Button
            variant="contained"
            startIcon={<Add />}
            onClick={async () => {
              await extendListItem();
            }}
          />
        </Toolbar>

        <Table>
          <TableHead>
            <TableRow>
              <TableCell style={{ maxWidth: "5rem" }}>Name</TableCell>
              <TableCell style={{ maxWidth: "5rem" }}>Type</TableCell>
              <TableCell>Description</TableCell>
              <TableCell style={{ maxWidth: "5rem" }}>Number in AU</TableCell>
              <TableCell>Sequence</TableCell>
              <TableCell style={{ maxWidth: "5rem" }}>Delete</TableCell>
            </TableRow>
          </TableHead>
          <TableBody>
            {item?._value?.map((contentElement: any, iElement: number) => {
              // Get validation color for the whole row (CAsuContentSeq element)
              const rowValidationColor = getValidationColor(contentElement);
              return (
                <TableRow
                  key={`${iElement}`}
                  onClick={() => handleOpenDialog(iElement)}
                  sx={{
                    transition: "box-shadow 0.2s",
                    cursor: "pointer",
                    // Apply left border based on row validation status
                    borderLeft: rowValidationColor !== "inherit" ? `4px solid ${rowValidationColor}` : undefined,
                    "&:hover": {
                      boxShadow: 3,
                      backgroundColor: "rgba(0, 0, 0, 0.04)",
                    },
                  }}
                >
                  {[
                    "name",
                    "polymerType",
                    "description",
                    "nCopies",
                    "sequence",
                  ].map((property) => {
                    // Get validation color for each individual field
                    const cellValidationColor = getValidationColor(contentElement._value[property]);
                    return (
                      <TableCell
                        key={property}
                        sx={{
                          maxWidth: ["name", "polymerType", "nCopies"].includes(
                            property
                          )
                            ? "5rem"
                            : property === "description"
                              ? "10rem"
                              : property === "sequence"
                                ? "20rem"
                                : undefined,
                          // Apply bottom border if field has validation error
                          borderBottom: cellValidationColor !== "inherit"
                            ? `3px solid ${cellValidationColor}`
                            : undefined,
                        }}
                      >
                        <div
                          style={{
                            maxHeight: "12rem",
                            overflowY: "auto",
                            wordWrap: "break-word",
                            whiteSpace: "pre-wrap",
                          }}
                        >
                          {contentElement._value[property]?._value ?? ""}
                        </div>
                      </TableCell>
                    );
                  })}
                  <TableCell style={{ maxWidth: "5rem" }}>
                    <Button
                      startIcon={<Delete />}
                      size="small"
                      onClick={(ev: any) => {
                        ev.stopPropagation();
                        ev.preventDefault();
                        deleteItem(iElement);
                      }}
                    />
                  </TableCell>
                </TableRow>
              );
            })}
          </TableBody>
        </Table>
        {(item?._value?.length ?? 0) === 0 && (
          <div
            style={{
              margin: "1rem auto",
              maxWidth: 340,
              background: "rgba(255,255,200,0.95)",
              borderRadius: 10,
              boxShadow: "0 1px 4px rgba(0,0,0,0.06)",
              padding: "0.75rem 1rem",
              textAlign: "center",
              color: "#666",
              fontSize: 15,
              lineHeight: 1.3,
              position: "relative",
            }}
          >
            <span
              role="img"
              aria-label="hint"
              style={{ fontSize: 20, marginBottom: 4 }}
            >
              💡
            </span>
            <div>
              No content yet. Use the <b>+</b> button above to add a constituent
              to the ASU.
            </div>
          </div>
        )}
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
