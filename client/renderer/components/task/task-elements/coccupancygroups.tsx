import React, { useCallback, useEffect, useMemo, useState } from "react";
import {
  Box,
  Button,
  IconButton,
  Stack,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  TextField,
  Tooltip,
  Typography,
} from "@mui/material";
import { Add, Delete } from "@mui/icons-material";
import { CCP4i2TaskElementProps } from "./task-element";
import { useJob, valueOfItem } from "../../../utils";

// ============================================================
// Shared utilities
// ============================================================

/** Validate that a list operation can proceed */
const validateListOperation = (item: any, job: any): boolean => {
  if (!item) return false;
  if (job.status !== 1) return false;
  return true;
};

/** Replace [?] placeholder with actual index in all nested object paths */
const updateObjectPath = (element: any, newIndex: number): any => {
  if (!element) return element;
  const updated = { ...element };
  updated._objectPath = element._objectPath?.replace(
    "[?]",
    `[${newIndex}]`
  );
  if (typeof updated._value === "object" && updated._value) {
    updated._value = Object.keys(updated._value).reduce(
      (acc: any, key: string) => {
        const val = updated._value[key];
        if (val?._objectPath) {
          acc[key] = {
            ...val,
            _objectPath: val._objectPath.replace("[?]", `[${newIndex}]`),
          };
        } else {
          acc[key] = val;
        }
        return acc;
      },
      {}
    );
  }
  return updated;
};

// ============================================================
// Inline cell editor
// ============================================================

interface InlineCellProps {
  objectPath: string;
  serverValue: any;
  type?: "number" | "text";
  placeholder?: string;
  disabled?: boolean;
  setParameter: (arg: { object_path: string; value: any }) => Promise<any>;
  sx?: any;
}

const InlineCell: React.FC<InlineCellProps> = ({
  objectPath,
  serverValue,
  type = "text",
  placeholder,
  disabled,
  setParameter,
  sx,
}) => {
  const displayValue =
    serverValue === null || serverValue === undefined ? "" : String(serverValue);
  const [localValue, setLocalValue] = useState(displayValue);

  // Sync from server when it changes
  useEffect(() => {
    const newDisplay =
      serverValue === null || serverValue === undefined
        ? ""
        : String(serverValue);
    setLocalValue(newDisplay);
  }, [serverValue]);

  const handleBlur = useCallback(async () => {
    const parsed = type === "number" ? Number(localValue) || 0 : localValue;
    const current =
      type === "number" ? Number(displayValue) || 0 : displayValue;
    if (parsed !== current) {
      await setParameter({ object_path: objectPath, value: parsed });
    }
  }, [localValue, displayValue, type, objectPath, setParameter]);

  const handleKeyDown = useCallback(
    (e: React.KeyboardEvent) => {
      if (e.key === "Enter") {
        (e.target as HTMLInputElement).blur();
      }
    },
    []
  );

  return (
    <TextField
      size="small"
      variant="standard"
      type={type}
      value={localValue}
      onChange={(e) => setLocalValue(e.target.value)}
      onBlur={handleBlur}
      onKeyDown={handleKeyDown}
      placeholder={placeholder}
      disabled={disabled}
      slotProps={{
        input: { sx: { fontSize: "0.875rem", ...sx } },
      }}
      fullWidth
    />
  );
};

// ============================================================
// COccRefmacSelectionListElement
//
// Renders COccRefmacSelectionList as a table:
// | Group ID | Chain(s) | From | To | Atom | Alt | [delete] |
// ============================================================

export const COccRefmacSelectionListElement: React.FC<
  CCP4i2TaskElementProps
> = (props) => {
  const { itemName, job, onChange } = props;

  const { useTaskItem, setParameter } = useJob(job.id);
  const { item } = useTaskItem(itemName);

  const isEditable = job.status === 1;
  const rows: any[] = item?._value ?? [];

  const handleAddRow = useCallback(async () => {
    if (!validateListOperation(item, job)) return;
    try {
      const taskElement = JSON.parse(JSON.stringify(item._subItem));
      const newIndex = rows.length;
      const updated = updateObjectPath(taskElement, newIndex);
      const currentListValue = Array.isArray(valueOfItem(item))
        ? [...valueOfItem(item)]
        : [];
      currentListValue.push(valueOfItem(updated));
      const result: any = await setParameter({
        object_path: item._objectPath,
        value: currentListValue,
      });
      if (result?.success && result.data?.updated_item && onChange) {
        await onChange(result.data.updated_item);
      }
    } catch (error) {
      console.error("Error adding occupancy group row:", error);
    }
  }, [item, job, rows.length, setParameter, onChange]);

  const handleDeleteRow = useCallback(
    async (index: number) => {
      if (!validateListOperation(item, job)) return;
      try {
        const currentListValue = Array.isArray(valueOfItem(item))
          ? [...valueOfItem(item)]
          : [];
        currentListValue.splice(index, 1);
        const result: any = await setParameter({
          object_path: item._objectPath,
          value: currentListValue,
        });
        if (result?.success && result.data?.updated_item && onChange) {
          await onChange(result.data.updated_item);
        }
      } catch (error) {
        console.error("Error deleting occupancy group row:", error);
      }
    },
    [item, job, setParameter, onChange]
  );

  if (!item) return null;

  return (
    <Box sx={{ mt: 1 }}>
      <TableContainer>
        <Table size="small">
          <TableHead>
            <TableRow>
              <TableCell sx={{ fontWeight: "bold", width: "5rem" }}>
                Group ID
              </TableCell>
              <TableCell sx={{ fontWeight: "bold", width: "6rem" }}>
                Chain(s)
              </TableCell>
              <TableCell sx={{ fontWeight: "bold", width: "5rem" }}>
                From
              </TableCell>
              <TableCell sx={{ fontWeight: "bold", width: "5rem" }}>
                To
              </TableCell>
              <TableCell sx={{ fontWeight: "bold", width: "5rem" }}>
                Atom
              </TableCell>
              <TableCell sx={{ fontWeight: "bold", width: "4rem" }}>
                Alt
              </TableCell>
              <TableCell sx={{ width: "3rem" }} />
            </TableRow>
          </TableHead>
          <TableBody>
            {rows.map((row: any, index: number) => {
              const v = row._value || {};
              return (
                <TableRow key={row._objectPath || index}>
                  <TableCell sx={{ py: 0.5 }}>
                    <InlineCell
                      objectPath={v.groupId?._objectPath}
                      serverValue={v.groupId?._value}
                      type="number"
                      placeholder="1"
                      disabled={!isEditable}
                      setParameter={setParameter}
                    />
                  </TableCell>
                  <TableCell sx={{ py: 0.5 }}>
                    <InlineCell
                      objectPath={v.chainIds?._objectPath}
                      serverValue={v.chainIds?._value}
                      placeholder="A"
                      disabled={!isEditable}
                      setParameter={setParameter}
                    />
                  </TableCell>
                  <TableCell sx={{ py: 0.5 }}>
                    <InlineCell
                      objectPath={v.firstRes?._objectPath}
                      serverValue={v.firstRes?._value}
                      type="number"
                      placeholder="1"
                      disabled={!isEditable}
                      setParameter={setParameter}
                    />
                  </TableCell>
                  <TableCell sx={{ py: 0.5 }}>
                    <InlineCell
                      objectPath={v.lastRes?._objectPath}
                      serverValue={v.lastRes?._value}
                      type="number"
                      placeholder="100"
                      disabled={!isEditable}
                      setParameter={setParameter}
                    />
                  </TableCell>
                  <TableCell sx={{ py: 0.5 }}>
                    <InlineCell
                      objectPath={v.atoms?._objectPath}
                      serverValue={v.atoms?._value}
                      placeholder="All"
                      disabled={!isEditable}
                      setParameter={setParameter}
                    />
                  </TableCell>
                  <TableCell sx={{ py: 0.5 }}>
                    <InlineCell
                      objectPath={v.alt?._objectPath}
                      serverValue={v.alt?._value}
                      placeholder="-"
                      disabled={!isEditable}
                      setParameter={setParameter}
                    />
                  </TableCell>
                  <TableCell sx={{ py: 0.5 }}>
                    <Tooltip title="Delete row">
                      <span>
                        <IconButton
                          size="small"
                          color="error"
                          disabled={!isEditable}
                          onClick={() => handleDeleteRow(index)}
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
      </TableContainer>

      {rows.length === 0 && (
        <Typography
          variant="body2"
          color="text.secondary"
          sx={{ fontStyle: "italic", mt: 1, ml: 1 }}
        >
          No occupancy groups defined.
        </Typography>
      )}

      <Box sx={{ mt: 1 }}>
        <Button
          variant="outlined"
          size="small"
          startIcon={<Add />}
          onClick={handleAddRow}
          disabled={!isEditable}
        >
          Add Group
        </Button>
      </Box>
    </Box>
  );
};

// ============================================================
// COccRelationRefmacListElement
//
// Renders COccRelationRefmacList as a simple table:
// | Group IDs (space separated) | [delete] |
// ============================================================

export const COccRelationRefmacListElement: React.FC<
  CCP4i2TaskElementProps
> = (props) => {
  const { itemName, job, onChange } = props;

  const { useTaskItem, setParameter } = useJob(job.id);
  const { item } = useTaskItem(itemName);

  const isEditable = job.status === 1;
  const rows: any[] = item?._value ?? [];
  const guiLabel =
    item?._objectPath?.split(".").at(-1) ||
    "Group relations";

  const handleAddRow = useCallback(async () => {
    if (!validateListOperation(item, job)) return;
    try {
      const taskElement = JSON.parse(JSON.stringify(item._subItem));
      const newIndex = rows.length;
      const updated = updateObjectPath(taskElement, newIndex);
      const currentListValue = Array.isArray(valueOfItem(item))
        ? [...valueOfItem(item)]
        : [];
      currentListValue.push(valueOfItem(updated));
      const result: any = await setParameter({
        object_path: item._objectPath,
        value: currentListValue,
      });
      if (result?.success && result.data?.updated_item && onChange) {
        await onChange(result.data.updated_item);
      }
    } catch (error) {
      console.error("Error adding group relation row:", error);
    }
  }, [item, job, rows.length, setParameter, onChange]);

  const handleDeleteRow = useCallback(
    async (index: number) => {
      if (!validateListOperation(item, job)) return;
      try {
        const currentListValue = Array.isArray(valueOfItem(item))
          ? [...valueOfItem(item)]
          : [];
        currentListValue.splice(index, 1);
        const result: any = await setParameter({
          object_path: item._objectPath,
          value: currentListValue,
        });
        if (result?.success && result.data?.updated_item && onChange) {
          await onChange(result.data.updated_item);
        }
      } catch (error) {
        console.error("Error deleting group relation row:", error);
      }
    },
    [item, job, setParameter, onChange]
  );

  if (!item) return null;

  return (
    <Box sx={{ mt: 1 }}>
      <TableContainer>
        <Table size="small">
          <TableHead>
            <TableRow>
              <TableCell sx={{ fontWeight: "bold" }}>
                Group IDs (space separated)
              </TableCell>
              <TableCell sx={{ width: "3rem" }} />
            </TableRow>
          </TableHead>
          <TableBody>
            {rows.map((row: any, index: number) => {
              const v = row._value || {};
              return (
                <TableRow key={row._objectPath || index}>
                  <TableCell sx={{ py: 0.5 }}>
                    <InlineCell
                      objectPath={v.groupIds?._objectPath}
                      serverValue={v.groupIds?._value}
                      placeholder="e.g. 1 2 3"
                      disabled={!isEditable}
                      setParameter={setParameter}
                    />
                  </TableCell>
                  <TableCell sx={{ py: 0.5 }}>
                    <Tooltip title="Delete row">
                      <span>
                        <IconButton
                          size="small"
                          color="error"
                          disabled={!isEditable}
                          onClick={() => handleDeleteRow(index)}
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
      </TableContainer>

      {rows.length === 0 && (
        <Typography
          variant="body2"
          color="text.secondary"
          sx={{ fontStyle: "italic", mt: 1, ml: 1 }}
        >
          No group relations defined.
        </Typography>
      )}

      <Box sx={{ mt: 1 }}>
        <Button
          variant="outlined"
          size="small"
          startIcon={<Add />}
          onClick={handleAddRow}
          disabled={!isEditable}
        >
          Add Relation
        </Button>
      </Box>
    </Box>
  );
};
