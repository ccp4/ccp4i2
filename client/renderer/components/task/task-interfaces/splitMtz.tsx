import { useCallback, useMemo, useState } from "react";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useJob, valueOfItem } from "../../../utils";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import {
  Alert,
  Box,
  Button,
  Chip,
  CircularProgress,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  Divider,
  FormControl,
  IconButton,
  InputLabel,
  List,
  ListItemButton,
  ListItemIcon,
  ListItemText,
  MenuItem,
  Paper,
  Select,
  Stack,
  Typography,
} from "@mui/material";
import ChevronLeftIcon from "@mui/icons-material/ChevronLeft";
import ChevronRightIcon from "@mui/icons-material/ChevronRight";
import KeyboardDoubleArrowLeftIcon from "@mui/icons-material/KeyboardDoubleArrowLeft";
import KeyboardDoubleArrowRightIcon from "@mui/icons-material/KeyboardDoubleArrowRight";
import TableChartIcon from "@mui/icons-material/TableChart";
import AddIcon from "@mui/icons-material/Add";
import DeleteIcon from "@mui/icons-material/Delete";
import BuildIcon from "@mui/icons-material/Build";

// Type definitions
interface MtzColumn {
  columnLabel: string;
  columnType: string;
  dataset: string;
  groupIndex: number;
}

interface ColumnGroup {
  columnGroupType: string;
  contentFlag: number;
  dataset: string;
  columnList: MtzColumn[];
}

/**
 * Generate a unique key for a column group for comparison purposes.
 * Uses dataset + columnGroupType + column labels to identify a group.
 */
function getGroupKey(group: ColumnGroup): string {
  const colLabels = group.columnList.map((c) => c.columnLabel).join(",");
  return `${group.dataset}:${group.columnGroupType}:${colLabels}`;
}

/**
 * Pattern definitions for column type signatures.
 * Format: [typeSignature, groupType, contentFlag]
 *
 * contentFlag values must match backend CONTENT_SIGNATURE_LIST order in CObsDataFileStub:
 * - Obs: 1=Anom I (KMKM), 2=Anom SF (GLGL), 3=Mean I (JQ), 4=Mean SF (FQ)
 * - Phs: 1=HL (AAAA), 2=Phi-FOM (PW)
 * - MapCoeffs: 1=F,Phi (FP), 2=F,sigF,Phi (FQP)
 * - FreeR: 1=Integer (I)
 */
const COLUMN_PATTERNS: [string, string, number][] = [
  // Observations - longer patterns first (anomalous before mean)
  ["KMKM", "Obs", 1], // I+, sigI+, I-, sigI- (Anomalous Intensities)
  ["GLGL", "Obs", 2], // F+, sigF+, F-, sigF- (Anomalous SFs)
  ["JQ", "Obs", 3], // I, sigI (Mean Intensities)
  ["FQ", "Obs", 4], // F, sigF (Mean SFs)
  // Phases
  ["AAAA", "Phs", 1], // HLA, HLB, HLC, HLD (Hendrickson-Lattman)
  ["PW", "Phs", 2], // Phi, FOM
  // Map coefficients
  ["FQP", "MapCoeffs", 2], // F, sigF, Phi
  ["FP", "MapCoeffs", 1], // F, Phi
  // FreeR flag
  ["I", "FreeR", 1], // Integer flag
];

/**
 * Group columns using pattern matching on column type signatures.
 */
function groupColumnsByPattern(columns: MtzColumn[]): ColumnGroup[] {
  if (!columns || columns.length === 0) return [];

  const groups: ColumnGroup[] = [];
  const used = new Array(columns.length).fill(false);
  const typeString = columns.map((c) => c.columnType).join("");

  let i = 0;
  while (i < columns.length) {
    if (used[i]) {
      i++;
      continue;
    }

    let matched = false;

    for (const [pattern, groupType, contentFlag] of COLUMN_PATTERNS) {
      const patternLen = pattern.length;
      if (i + patternLen > columns.length) continue;

      const candidate = typeString.slice(i, i + patternLen);
      if (candidate !== pattern) continue;

      // Special check for FreeR: label must contain 'free'
      if (groupType === "FreeR") {
        const label = columns[i].columnLabel;
        if (!label.toLowerCase().includes("free")) continue;
      }

      // For Obs data, columns must come from the same dataset
      if (groupType === "Obs") {
        const firstDataset = columns[i].dataset;
        if (firstDataset) {
          const sameDataset = Array.from({ length: patternLen }).every(
            (_, j) => columns[i + j].dataset === firstDataset
          );
          if (!sameDataset) continue;
        }
      }

      groups.push({
        columnGroupType: groupType,
        contentFlag,
        dataset: columns[i].dataset || "",
        columnList: columns.slice(i, i + patternLen),
      });

      for (let j = 0; j < patternLen; j++) {
        used[i + j] = true;
      }

      i += patternLen;
      matched = true;
      break;
    }

    if (!matched) {
      i++;
    }
  }

  return groups;
}

interface MtzDigest {
  columnGroups?: ColumnGroup[];
  listOfColumns?: MtzColumn[];
  format?: string;
  merged?: boolean;
  cell?: {
    a?: number;
    b?: number;
    c?: number;
    alpha?: number;
    beta?: number;
    gamma?: number;
  };
  spaceGroup?: string;
  status?: string;
  reason?: string;
}

// Color mapping for column group types
const TYPE_COLORS: Record<string, "primary" | "secondary" | "success" | "warning" | "info" | "default"> = {
  Obs: "primary",
  Phs: "secondary",
  MapCoeffs: "success",
  FreeR: "warning",
};

const TYPE_LABELS: Record<string, string> = {
  Obs: "Observations",
  Phs: "Phases",
  MapCoeffs: "Map Coefficients",
  FreeR: "Free R Flag",
};

interface ColumnRequirement {
  types: string[];
  labels: string[];
  description: string;
}

/**
 * Column type requirements for each file type + contentFlag combination.
 * Used by the custom column group builder to validate user selections.
 *
 * contentFlag values must match backend CONTENT_SIGNATURE_LIST order:
 * - Obs: 1=Anom I, 2=Anom SF, 3=Mean I, 4=Mean SF
 * - Phs: 1=HL, 2=Phi-FOM
 * - MapCoeffs: 1=F,Phi, 2=F,sigF,Phi
 * - FreeR: 1=Integer
 */
const COLUMN_REQUIREMENTS: Record<string, Record<number, ColumnRequirement>> = {
  Obs: {
    1: { types: ["K", "M", "K", "M"], labels: ["I+", "sigI+", "I-", "sigI-"], description: "Anomalous Is" },
    2: { types: ["G", "L", "G", "L"], labels: ["F+", "sigF+", "F-", "sigF-"], description: "Anomalous SFs" },
    3: { types: ["J", "Q"], labels: ["I", "SIGI"], description: "Mean Is (I, sigI)" },
    4: { types: ["F", "Q"], labels: ["F", "SIGF"], description: "Mean SFs (F, sigF)" },
  },
  MapCoeffs: {
    1: { types: ["F", "P"], labels: ["F", "PHI"], description: "F, Phi" },
    2: { types: ["F", "Q", "P"], labels: ["F", "SIGF", "PHI"], description: "F, sigF, Phi" },
  },
  Phs: {
    1: { types: ["A", "A", "A", "A"], labels: ["HLA", "HLB", "HLC", "HLD"], description: "Hendrickson-Lattman" },
    2: { types: ["P", "W"], labels: ["PHI", "FOM"], description: "Phi-FOM" },
  },
  FreeR: {
    1: { types: ["I"], labels: ["FREE"], description: "Integer flag" },
  },
};

/**
 * Get available content flags for a given file type.
 */
function getContentFlagsForType(fileType: string): { value: number; label: string }[] {
  const requirements = COLUMN_REQUIREMENTS[fileType];
  if (!requirements) return [];
  return Object.entries(requirements).map(([flag, req]) => ({
    value: parseInt(flag),
    label: req.description,
  }));
}

/**
 * Helper to extract value, handling _value wrapper from JSON encoder.
 */
function extractValue(val: any): any {
  if (val && typeof val === "object" && "_value" in val) {
    return val._value;
  }
  return val;
}

/**
 * Convert container item value to ColumnGroup array.
 * Handles both direct values and _value wrappers from JSON encoder.
 */
function containerValueToColumnGroups(containerValue: any[]): ColumnGroup[] {
  if (!containerValue || !Array.isArray(containerValue)) return [];

  return containerValue.map((item) => {
    const columnListRaw = extractValue(item.columnList) || [];
    const columnList = Array.isArray(columnListRaw) ? columnListRaw : [];

    return {
      columnGroupType: extractValue(item.columnGroupType) || "",
      contentFlag: extractValue(item.contentFlag) || 0,
      dataset: extractValue(item.dataset) || "",
      columnList: columnList.map((col: any) => ({
        columnLabel: extractValue(col.columnLabel) || "",
        columnType: extractValue(col.columnType) || "",
        dataset: extractValue(col.dataset) || "",
        groupIndex: extractValue(col.groupIndex) || 0,
      })),
    };
  });
}

/**
 * Convert ColumnGroup array to the format expected by set_parameter.
 * Items in COLUMNGROUPLIST are marked as selected=true for the backend to process.
 */
function columnGroupsToContainerValue(groups: ColumnGroup[]): any[] {
  return groups.map((group) => ({
    columnGroupType: group.columnGroupType,
    contentFlag: group.contentFlag,
    dataset: group.dataset,
    selected: true, // Mark as selected for backend processing
    columnList: group.columnList.map((col) => ({
      columnLabel: col.columnLabel,
      columnType: col.columnType,
      dataset: col.dataset,
      groupIndex: col.groupIndex,
    })),
  }));
}

/**
 * ColumnGroupItem - Single item in the available/selected lists
 */
const ColumnGroupItem: React.FC<{
  group: ColumnGroup;
  selected: boolean;
  onClick: () => void;
  disabled: boolean;
}> = ({ group, selected, onClick, disabled }) => {
  const chipColor = TYPE_COLORS[group.columnGroupType] || "default";
  const typeLabel = TYPE_LABELS[group.columnGroupType] || group.columnGroupType;
  const columnsStr = group.columnList.map((c) => c.columnLabel).join(", ");

  return (
    <ListItemButton
      selected={selected}
      onClick={onClick}
      disabled={disabled}
      sx={{ py: 0.5 }}
    >
      <ListItemIcon sx={{ minWidth: 36 }}>
        <TableChartIcon fontSize="small" color={chipColor as any} />
      </ListItemIcon>
      <ListItemText
        primary={
          <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
            <Chip
              label={typeLabel}
              color={chipColor}
              size="small"
              sx={{ height: 20, fontSize: "0.7rem" }}
            />
            <Typography variant="body2" color="text.secondary">
              {group.dataset || "-"}
            </Typography>
          </Box>
        }
        secondary={
          <Typography
            variant="caption"
            sx={{ fontFamily: "monospace", color: "text.secondary" }}
          >
            {columnsStr}
          </Typography>
        }
      />
    </ListItemButton>
  );
};

/**
 * TwoColumnSelector - Displays available and selected column groups side by side
 */
const TwoColumnSelector: React.FC<{
  availableGroups: ColumnGroup[];
  selectedGroups: ColumnGroup[];
  onSelect: (group: ColumnGroup) => void;
  onDeselect: (group: ColumnGroup) => void;
  onSelectAll: () => void;
  onDeselectAll: () => void;
  disabled: boolean;
}> = ({ availableGroups, selectedGroups, onSelect, onDeselect, onSelectAll, onDeselectAll, disabled }) => {
  // Track which items are highlighted for transfer
  const [highlightedAvailable, setHighlightedAvailable] = useState<Set<string>>(new Set());
  const [highlightedSelected, setHighlightedSelected] = useState<Set<string>>(new Set());

  const toggleAvailableHighlight = (key: string) => {
    setHighlightedAvailable((prev) => {
      const next = new Set(prev);
      if (next.has(key)) {
        next.delete(key);
      } else {
        next.add(key);
      }
      return next;
    });
  };

  const toggleSelectedHighlight = (key: string) => {
    setHighlightedSelected((prev) => {
      const next = new Set(prev);
      if (next.has(key)) {
        next.delete(key);
      } else {
        next.add(key);
      }
      return next;
    });
  };

  const handleMoveRight = () => {
    // Move highlighted available items to selected
    availableGroups
      .filter((g) => highlightedAvailable.has(getGroupKey(g)))
      .forEach((g) => onSelect(g));
    setHighlightedAvailable(new Set());
  };

  const handleMoveLeft = () => {
    // Move highlighted selected items to available
    selectedGroups
      .filter((g) => highlightedSelected.has(getGroupKey(g)))
      .forEach((g) => onDeselect(g));
    setHighlightedSelected(new Set());
  };

  const canMoveRight = highlightedAvailable.size > 0;
  const canMoveLeft = highlightedSelected.size > 0;

  return (
    <Stack direction="row" spacing={1} sx={{ mt: 2 }}>
      {/* Available column groups */}
      <Paper variant="outlined" sx={{ flex: 1, minHeight: 200 }}>
        <Typography
          variant="subtitle2"
          sx={{ p: 1, bgcolor: "action.hover", borderBottom: 1, borderColor: "divider" }}
        >
          Available ({availableGroups.length})
        </Typography>
        <List dense sx={{ maxHeight: 300, overflow: "auto" }}>
          {availableGroups.length === 0 ? (
            <Typography variant="body2" color="text.secondary" sx={{ p: 2, textAlign: "center" }}>
              No available groups
            </Typography>
          ) : (
            availableGroups.map((group) => {
              const key = getGroupKey(group);
              return (
                <ColumnGroupItem
                  key={key}
                  group={group}
                  selected={highlightedAvailable.has(key)}
                  onClick={() => toggleAvailableHighlight(key)}
                  disabled={disabled}
                />
              );
            })
          )}
        </List>
      </Paper>

      {/* Transfer buttons */}
      <Stack justifyContent="center" spacing={1}>
        <Button
          variant="outlined"
          size="small"
          onClick={onSelectAll}
          disabled={disabled || availableGroups.length === 0}
          sx={{ minWidth: 40, px: 1 }}
          title="Select all"
        >
          <KeyboardDoubleArrowRightIcon />
        </Button>
        <Button
          variant="outlined"
          size="small"
          onClick={handleMoveRight}
          disabled={disabled || !canMoveRight}
          sx={{ minWidth: 40, px: 1 }}
          title="Select highlighted"
        >
          <ChevronRightIcon />
        </Button>
        <Button
          variant="outlined"
          size="small"
          onClick={handleMoveLeft}
          disabled={disabled || !canMoveLeft}
          sx={{ minWidth: 40, px: 1 }}
          title="Deselect highlighted"
        >
          <ChevronLeftIcon />
        </Button>
        <Button
          variant="outlined"
          size="small"
          onClick={onDeselectAll}
          disabled={disabled || selectedGroups.length === 0}
          sx={{ minWidth: 40, px: 1 }}
          title="Deselect all"
        >
          <KeyboardDoubleArrowLeftIcon />
        </Button>
      </Stack>

      {/* Selected column groups */}
      <Paper variant="outlined" sx={{ flex: 1, minHeight: 200 }}>
        <Typography
          variant="subtitle2"
          sx={{ p: 1, bgcolor: "success.light", color: "success.contrastText", borderBottom: 1, borderColor: "divider" }}
        >
          Selected for Export ({selectedGroups.length})
        </Typography>
        <List dense sx={{ maxHeight: 300, overflow: "auto" }}>
          {selectedGroups.length === 0 ? (
            <Typography variant="body2" color="text.secondary" sx={{ p: 2, textAlign: "center" }}>
              No groups selected
            </Typography>
          ) : (
            selectedGroups.map((group) => {
              const key = getGroupKey(group);
              return (
                <ColumnGroupItem
                  key={key}
                  group={group}
                  selected={highlightedSelected.has(key)}
                  onClick={() => toggleSelectedHighlight(key)}
                  disabled={disabled}
                />
              );
            })
          )}
        </List>
      </Paper>
    </Stack>
  );
};

/**
 * MtzFileInfo - Display summary info about the loaded MTZ file
 */
const MtzFileInfo: React.FC<{ digest: MtzDigest }> = ({ digest }) => {
  if (!digest || digest.status === "Failed") {
    return null;
  }

  const cell = digest.cell;
  const cellString = cell
    ? `${cell.a?.toFixed(2)} ${cell.b?.toFixed(2)} ${cell.c?.toFixed(2)} ${cell.alpha?.toFixed(1)} ${cell.beta?.toFixed(1)} ${cell.gamma?.toFixed(1)}`
    : "N/A";

  return (
    <Box sx={{ mt: 1, mb: 2 }}>
      <Typography variant="body2" color="text.secondary">
        <strong>Format:</strong> {digest.format || "MTZ"} |{" "}
        <strong>Merged:</strong> {digest.merged ? "Yes" : "No"} |{" "}
        <strong>Space Group:</strong> {digest.spaceGroup || "N/A"} |{" "}
        <strong>Cell:</strong> {cellString}
      </Typography>
    </Box>
  );
};

/**
 * CustomColumnGroupDialog - Modal for creating custom column groups
 */
const CustomColumnGroupDialog: React.FC<{
  open: boolean;
  onClose: () => void;
  onAdd: (group: ColumnGroup) => void;
  onDelete: (group: ColumnGroup) => void;
  availableColumns: MtzColumn[];
  customGroups: ColumnGroup[];
  disabled: boolean;
}> = ({ open, onClose, onAdd, onDelete, availableColumns, customGroups, disabled }) => {
  const [fileType, setFileType] = useState<string>("");
  const [contentFlag, setContentFlag] = useState<number>(0);
  const [selectedColumns, setSelectedColumns] = useState<(MtzColumn | null)[]>([]);
  const [dataset, setDataset] = useState<string>("");

  // Get available datasets from columns
  const datasets = useMemo(() => {
    const unique = new Set(availableColumns.map((c) => c.dataset).filter(Boolean));
    return Array.from(unique);
  }, [availableColumns]);

  // Get requirements for current selection
  const requirements = useMemo(() => {
    if (!fileType || !contentFlag) return null;
    return COLUMN_REQUIREMENTS[fileType]?.[contentFlag] || null;
  }, [fileType, contentFlag]);

  // Reset column selections when requirements change
  const handleFileTypeChange = (newType: string) => {
    setFileType(newType);
    setContentFlag(0);
    setSelectedColumns([]);
  };

  const handleContentFlagChange = (newFlag: number) => {
    setContentFlag(newFlag);
    const req = COLUMN_REQUIREMENTS[fileType]?.[newFlag];
    if (req) {
      setSelectedColumns(new Array(req.types.length).fill(null));
    } else {
      setSelectedColumns([]);
    }
  };

  // Get columns matching a specific type
  const getColumnsForType = (typeCode: string): MtzColumn[] => {
    return availableColumns.filter((col) => col.columnType === typeCode);
  };

  // Check if the current selection is valid
  const isValid = useMemo(() => {
    if (!requirements) return false;
    if (selectedColumns.length !== requirements.types.length) return false;
    return selectedColumns.every((col, idx) => {
      if (!col) return false;
      return col.columnType === requirements.types[idx];
    });
  }, [requirements, selectedColumns]);

  // Handle adding the custom group
  const handleAdd = () => {
    if (!isValid || !requirements) return;

    const newGroup: ColumnGroup = {
      columnGroupType: fileType,
      contentFlag: contentFlag,
      dataset: dataset || selectedColumns[0]?.dataset || "",
      columnList: selectedColumns.filter((c): c is MtzColumn => c !== null),
    };

    onAdd(newGroup);

    // Reset form
    setFileType("");
    setContentFlag(0);
    setSelectedColumns([]);
    setDataset("");
  };

  const fileTypes = Object.keys(COLUMN_REQUIREMENTS);
  const contentFlags = getContentFlagsForType(fileType);

  return (
    <Dialog open={open} onClose={onClose} maxWidth="md" fullWidth>
      <DialogTitle>
        <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
          <BuildIcon />
          Custom Column Groups
        </Box>
      </DialogTitle>
      <DialogContent>
        {/* Existing custom groups */}
        {customGroups.length > 0 && (
          <Box sx={{ mb: 3 }}>
            <Typography variant="subtitle2" sx={{ mb: 1 }}>
              Existing Custom Groups
            </Typography>
            <Paper variant="outlined">
              <List dense>
                {customGroups.map((group, idx) => {
                  const key = `custom-${idx}-${getGroupKey(group)}`;
                  const chipColor = TYPE_COLORS[group.columnGroupType] || "default";
                  const typeLabel = TYPE_LABELS[group.columnGroupType] || group.columnGroupType;
                  const columnsStr = group.columnList.map((c) => c.columnLabel).join(", ");
                  return (
                    <ListItemButton key={key} sx={{ py: 0.5 }}>
                      <ListItemIcon sx={{ minWidth: 36 }}>
                        <TableChartIcon fontSize="small" color={chipColor as any} />
                      </ListItemIcon>
                      <ListItemText
                        primary={
                          <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
                            <Chip label={typeLabel} color={chipColor} size="small" sx={{ height: 20, fontSize: "0.7rem" }} />
                            <Chip label="Custom" size="small" variant="outlined" sx={{ height: 20, fontSize: "0.7rem" }} />
                            <Typography variant="body2" color="text.secondary">
                              {group.dataset || "-"}
                            </Typography>
                          </Box>
                        }
                        secondary={
                          <Typography variant="caption" sx={{ fontFamily: "monospace", color: "text.secondary" }}>
                            {columnsStr}
                          </Typography>
                        }
                      />
                      <IconButton
                        size="small"
                        onClick={() => onDelete(group)}
                        disabled={disabled}
                        title="Remove custom group"
                      >
                        <DeleteIcon fontSize="small" />
                      </IconButton>
                    </ListItemButton>
                  );
                })}
              </List>
            </Paper>
          </Box>
        )}

        {/* Create new custom group */}
        <Typography variant="subtitle2" sx={{ mb: 2 }}>
          Create New Custom Group
        </Typography>

        <Stack spacing={2}>
          {/* File type selector */}
          <FormControl fullWidth size="small">
            <InputLabel>File Type</InputLabel>
            <Select
              value={fileType}
              label="File Type"
              onChange={(e) => handleFileTypeChange(e.target.value)}
              disabled={disabled}
            >
              {fileTypes.map((type) => (
                <MenuItem key={type} value={type}>
                  {TYPE_LABELS[type] || type}
                </MenuItem>
              ))}
            </Select>
          </FormControl>

          {/* Content flag selector */}
          {fileType && (
            <FormControl fullWidth size="small">
              <InputLabel>Data Type</InputLabel>
              <Select
                value={contentFlag || ""}
                label="Data Type"
                onChange={(e) => handleContentFlagChange(Number(e.target.value))}
                disabled={disabled}
              >
                {contentFlags.map((cf) => (
                  <MenuItem key={cf.value} value={cf.value}>
                    {cf.label}
                  </MenuItem>
                ))}
              </Select>
            </FormControl>
          )}

          {/* Dataset selector (optional) */}
          {requirements && datasets.length > 0 && (
            <FormControl fullWidth size="small">
              <InputLabel>Dataset (optional)</InputLabel>
              <Select
                value={dataset}
                label="Dataset (optional)"
                onChange={(e) => setDataset(e.target.value)}
                disabled={disabled}
              >
                <MenuItem value="">
                  <em>Auto-detect</em>
                </MenuItem>
                {datasets.map((ds) => (
                  <MenuItem key={ds} value={ds}>
                    {ds}
                  </MenuItem>
                ))}
              </Select>
            </FormControl>
          )}

          {/* Column selectors */}
          {requirements && (
            <Box>
              <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                Select columns for each slot:
              </Typography>
              <Stack spacing={1}>
                {requirements.labels.map((label, idx) => {
                  const requiredType = requirements.types[idx];
                  const matchingColumns = getColumnsForType(requiredType);
                  const selectedCol = selectedColumns[idx];

                  return (
                    <FormControl key={idx} fullWidth size="small">
                      <InputLabel>
                        {label} (type: {requiredType})
                      </InputLabel>
                      <Select
                        value={selectedCol?.columnLabel || ""}
                        label={`${label} (type: ${requiredType})`}
                        onChange={(e) => {
                          const col = matchingColumns.find((c) => c.columnLabel === e.target.value);
                          const newSelected = [...selectedColumns];
                          newSelected[idx] = col || null;
                          setSelectedColumns(newSelected);
                        }}
                        disabled={disabled}
                        error={matchingColumns.length === 0}
                      >
                        {matchingColumns.length === 0 ? (
                          <MenuItem disabled>
                            <em>No columns of type {requiredType} available</em>
                          </MenuItem>
                        ) : (
                          matchingColumns.map((col) => (
                            <MenuItem key={col.columnLabel} value={col.columnLabel}>
                              {col.columnLabel} ({col.dataset || "no dataset"})
                            </MenuItem>
                          ))
                        )}
                      </Select>
                    </FormControl>
                  );
                })}
              </Stack>
            </Box>
          )}
        </Stack>
      </DialogContent>
      <DialogActions>
        <Button onClick={onClose}>Close</Button>
        <Button
          variant="contained"
          startIcon={<AddIcon />}
          onClick={handleAdd}
          disabled={disabled || !isValid}
        >
          Add to Selection
        </Button>
      </DialogActions>
    </Dialog>
  );
};

/**
 * SplitMtz Task Interface
 *
 * Architecture:
 * - SUBSTRATE (Declarative): Column groups from digest - read-only reference data
 * - SELECTION (Imperative): COLUMNGROUPLIST contains only selected groups
 * - UI: Two-column selector with Available | Selected panels
 *
 * Transfer buttons trigger imperative set_param calls to update COLUMNGROUPLIST.
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem, useFileDigest } = useJob(job.id);

  // Task item hooks
  const { item: hklinItem } = useTaskItem("HKLIN");
  const { value: containerColumnGroups, update: updateColumnGroupList } =
    useTaskItem("COLUMNGROUPLIST");

  const isEditable = job.status === 1;

  // Custom group builder dialog state
  const [customDialogOpen, setCustomDialogOpen] = useState(false);
  const [customGroups, setCustomGroups] = useState<ColumnGroup[]>([]);

  // Extract file UUID for display purposes and digest guard
  const hklinValue = useMemo(() => {
    if (!hklinItem) return null;
    return valueOfItem(hklinItem);
  }, [hklinItem]);

  const fileUuid = hklinValue?.dbFileId || null;

  // Only fetch digest when a file has been uploaded (has dbFileId)
  const hasFile = Boolean(fileUuid);
  const digestPath = hasFile ? `splitMtz.inputData.HKLIN` : "";
  const {
    data: digest,
    error: digestError,
    isLoading: loading,
    mutate: mutateDigest,
  } = useFileDigest(digestPath);

  // Get raw columns from digest for custom group builder
  const availableColumns = useMemo(() => {
    if (!digest?.listOfColumns) return [];
    return digest.listOfColumns;
  }, [digest]);

  // Compute substrate groups from digest
  const substrateGroups = useMemo(() => {
    if (!digest) return [];

    if (digest.listOfColumns && digest.listOfColumns.length > 0) {
      return groupColumnsByPattern(digest.listOfColumns);
    } else if (digest.columnGroups) {
      return digest.columnGroups.map((grp: ColumnGroup) => ({ ...grp }));
    }
    return [];
  }, [digest]);

  // Combined substrate = auto-detected groups + custom groups
  const combinedSubstrate = useMemo(() => {
    return [...substrateGroups, ...customGroups];
  }, [substrateGroups, customGroups]);

  // Selected groups come from container (COLUMNGROUPLIST)
  const selectedGroups = useMemo(() => {
    return containerValueToColumnGroups(containerColumnGroups || []);
  }, [containerColumnGroups]);

  // Available groups = combined substrate minus selected
  const availableGroups = useMemo(() => {
    const selectedKeys = new Set(selectedGroups.map(getGroupKey));
    return combinedSubstrate.filter((g) => !selectedKeys.has(getGroupKey(g)));
  }, [combinedSubstrate, selectedGroups]);

  /**
   * Handle HKLIN file change - clear selection and trigger digest refresh.
   * Selection must be cleared because column groups from the old file are no longer valid.
   */
  const handleHklinChange = useCallback(async () => {
    // Clear the selection - old column groups are no longer valid
    await updateColumnGroupList([]);
    // Clear custom groups too - they reference columns from the old file
    setCustomGroups([]);
    // Trigger SWR to revalidate/refetch the digest
    mutateDigest();
  }, [mutateDigest, updateColumnGroupList]);

  /**
   * Handle selecting a group (move from available to selected).
   * Imperatively updates COLUMNGROUPLIST.
   */
  const handleSelectGroup = useCallback(
    async (group: ColumnGroup) => {
      const newSelected = [...selectedGroups, group];
      await updateColumnGroupList(columnGroupsToContainerValue(newSelected));
    },
    [selectedGroups, updateColumnGroupList]
  );

  /**
   * Handle deselecting a group (move from selected to available).
   * Imperatively updates COLUMNGROUPLIST.
   */
  const handleDeselectGroup = useCallback(
    async (group: ColumnGroup) => {
      const keyToRemove = getGroupKey(group);
      const newSelected = selectedGroups.filter((g) => getGroupKey(g) !== keyToRemove);
      await updateColumnGroupList(columnGroupsToContainerValue(newSelected));
    },
    [selectedGroups, updateColumnGroupList]
  );

  /**
   * Handle selecting all available groups.
   */
  const handleSelectAll = useCallback(async () => {
    const allGroups = [...selectedGroups, ...availableGroups];
    await updateColumnGroupList(columnGroupsToContainerValue(allGroups));
  }, [selectedGroups, availableGroups, updateColumnGroupList]);

  /**
   * Handle deselecting all groups (clear selection).
   */
  const handleDeselectAll = useCallback(async () => {
    await updateColumnGroupList([]);
  }, [updateColumnGroupList]);

  /**
   * Handle adding a custom group - adds to both custom groups list and selected groups.
   */
  const handleAddCustomGroup = useCallback(
    async (group: ColumnGroup) => {
      // Add to custom groups list
      setCustomGroups((prev) => [...prev, group]);
      // Also add to selected groups
      const newSelected = [...selectedGroups, group];
      await updateColumnGroupList(columnGroupsToContainerValue(newSelected));
    },
    [selectedGroups, updateColumnGroupList]
  );

  /**
   * Handle deleting a custom group - removes from both custom groups list and selected groups.
   */
  const handleDeleteCustomGroup = useCallback(
    async (group: ColumnGroup) => {
      const keyToRemove = getGroupKey(group);
      // Remove from custom groups list
      setCustomGroups((prev) => prev.filter((g) => getGroupKey(g) !== keyToRemove));
      // Also remove from selected groups if present
      const newSelected = selectedGroups.filter((g) => getGroupKey(g) !== keyToRemove);
      await updateColumnGroupList(columnGroupsToContainerValue(newSelected));
    },
    [selectedGroups, updateColumnGroupList]
  );

  return (
    <CCP4i2Tabs {...props}>
      <CCP4i2Tab label="Input" key="input">
        <CCP4i2ContainerElement
          {...props}
          key="FileInput"
          itemName=""
          containerHint="BlockLevel"
          qualifiers={{ initiallyOpen: true, guiLabel: "MTZ File Input" }}
        >
          <CCP4i2TaskElement
            {...props}
            key="HKLIN"
            itemName="HKLIN"
            qualifiers={{
              guiLabel: "Input MTZ file",
              toolTip: "Select an MTZ file to split into separate column groups",
            }}
            onChange={handleHklinChange}
          />
        </CCP4i2ContainerElement>

        <CCP4i2ContainerElement
          {...props}
          key="ColumnGroups"
          itemName=""
          containerHint="BlockLevel"
          qualifiers={{
            initiallyOpen: true,
            guiLabel: `Column Groups${selectedGroups.length > 0 ? ` (${selectedGroups.length} selected)` : ""}`,
          }}
        >
          {loading && (
            <Box sx={{ display: "flex", alignItems: "center", gap: 2, py: 2 }}>
              <CircularProgress size={20} />
              <Typography>Loading MTZ file information...</Typography>
            </Box>
          )}

          {digestError && (
            <Alert severity="error" sx={{ mt: 1 }}>
              {digestError.message || "Failed to load MTZ file information"}
            </Alert>
          )}

          {!loading && !digestError && !fileUuid && (
            <Alert severity="info" sx={{ mt: 1 }}>
              Select an MTZ file above to view its column groups.
            </Alert>
          )}

          {!loading && !digestError && fileUuid && digest && (
            <>
              <MtzFileInfo digest={digest} />
              <Divider sx={{ my: 1 }} />
              <Box sx={{ display: "flex", justifyContent: "space-between", alignItems: "center", mb: 1 }}>
                <Typography variant="body2" color="text.secondary">
                  Click items to highlight, then use arrows to transfer between Available and Selected.
                </Typography>
                <Button
                  variant="outlined"
                  size="small"
                  startIcon={<BuildIcon />}
                  onClick={() => setCustomDialogOpen(true)}
                  disabled={!isEditable || availableColumns.length === 0}
                >
                  Custom Group{customGroups.length > 0 ? ` (${customGroups.length})` : ""}
                </Button>
              </Box>
              <TwoColumnSelector
                availableGroups={availableGroups}
                selectedGroups={selectedGroups}
                onSelect={handleSelectGroup}
                onDeselect={handleDeselectGroup}
                onSelectAll={handleSelectAll}
                onDeselectAll={handleDeselectAll}
                disabled={!isEditable}
              />
              {selectedGroups.length > 0 && (
                <Typography variant="body2" color="text.secondary" sx={{ mt: 2 }}>
                  {selectedGroups.length} column group{selectedGroups.length !== 1 ? "s" : ""} will be
                  exported as separate files.
                </Typography>
              )}

              {/* Custom Column Group Builder Dialog */}
              <CustomColumnGroupDialog
                open={customDialogOpen}
                onClose={() => setCustomDialogOpen(false)}
                onAdd={handleAddCustomGroup}
                onDelete={handleDeleteCustomGroup}
                availableColumns={availableColumns}
                customGroups={customGroups}
                disabled={!isEditable}
              />
            </>
          )}
        </CCP4i2ContainerElement>
      </CCP4i2Tab>
    </CCP4i2Tabs>
  );
};

export default TaskInterface;
