"use client";

import React, {
  SyntheticEvent,
  useCallback,
  useEffect,
  useMemo,
  useReducer,
  useState,
} from "react";
import {
  Autocomplete,
  AutocompleteChangeReason,
  Button,
  Checkbox,
  DialogActions,
  DialogContent,
  DialogTitle,
  Divider,
  FormControlLabel,
  Radio,
  RadioGroup,
  Stack,
  TextField,
  Typography,
} from "@mui/material";
import SimpleDialog from "@mui/material/Dialog";
import { createRoot } from "react-dom/client";
import { v4 as uuid4 } from "uuid";

// Constants
const SIGNATURE_MAP: Record<string, string> = {
  KMKM: "Intensity Friedel pairs (I+/I-)",
  GLGL: "Structure factor Friedel pairs (F+/F-)",
  JQ: "Mean intensities (I/σI)",
  FQ: "Mean structure factors (F/σF)",
  I: "FreeR flags",
} as const;

// Map signatures to content flags for CObsDataFile
const SIGNATURE_TO_CONTENT_FLAG: Record<string, number> = {
  KMKM: 1, // CONTENT_FLAG_IPAIR
  GLGL: 2, // CONTENT_FLAG_FPAIR
  JQ: 3,   // CONTENT_FLAG_IMEAN
  FQ: 4,   // CONTENT_FLAG_FMEAN
} as const;

// Map signatures to file suffixes for multi-representation output
const SIGNATURE_TO_SUFFIX: Record<string, string> = {
  KMKM: "_asIPAIR",
  GLGL: "_asFPAIR",
  JQ: "_asIMEAN",
  FQ: "_asFMEAN",
  I: "_asFREER",
} as const;

const GENERIC_SIGNATURES = ["KMKM", "GLGL", "JQ", "FQ"] as const;
const FREER_SIGNATURES = ["I"] as const;

// Types
interface ColumnOptions {
  [signature: string]: string[];
}

interface ColumnNames {
  [columnLabel: string]: string;
}

interface ColumnCounters {
  [columnType: string]: number;
}

interface ValuesState {
  [signature: string]: string;
}

interface ValuesAction {
  type: "SET_VALUE";
  signature: string;
  value: string;
}

interface ItemQualifiers {
  correctColumns?: string[];
}

interface MtzItem {
  _class?: string;
  _objectPath?: string;
  _qualifiers?: ItemQualifiers;
}

interface MtzColumnDialogProps {
  columnNames: ColumnNames;
  item: MtzItem;
  onAccept: (columnSelector: string) => void;
  onCancel: () => void;
}

// Reducer
const valuesReducer = (
  state: ValuesState,
  action: ValuesAction
): ValuesState => {
  switch (action.type) {
    case "SET_VALUE":
      return {
        ...state,
        [action.signature]: action.value,
      };
    default:
      return state;
  }
};

// Helper functions
const createSignatureOptions = (
  signature: string,
  sortedColumnNames: Record<string, string[]>
): string[] => {
  const signatureOptions: string[][] = [];
  let optionIndex = 0;

  const columnCounters: ColumnCounters = {};
  Object.keys(sortedColumnNames).forEach((key) => {
    columnCounters[key] = 0;
  });

  let shouldContinue = true;

  while (shouldContinue) {
    signatureOptions[optionIndex] = [];

    for (let i = 0; i < signature.length; i++) {
      const columnType = signature.charAt(i);

      if (!Object.keys(sortedColumnNames).includes(columnType)) {
        shouldContinue = false;
        break;
      }

      if (columnCounters[columnType] >= sortedColumnNames[columnType].length) {
        shouldContinue = false;
        break;
      }

      signatureOptions[optionIndex].push(
        sortedColumnNames[columnType][columnCounters[columnType]]
      );
      columnCounters[columnType] += 1;
    }

    optionIndex += 1;
  }

  return signatureOptions
    .filter((option) => option.length > 0)
    .map((option) => `/*/*/[${option.join(",")}]`);
};

const buildColumnOptions = (
  allColumnNames: ColumnNames,
  item: MtzItem
): ColumnOptions => {
  const sortedColumnNames: Record<string, string[]> = {};

  // Group columns by type
  Object.entries(allColumnNames).forEach(([label, columnType]) => {
    if (!sortedColumnNames[columnType]) {
      sortedColumnNames[columnType] = [];
    }
    sortedColumnNames[columnType].push(label);
  });

  // Determine signatures to process
  // CGenericReflDataFile is a general container without specific column types
  // so it uses the generic signatures. Other types use their correctColumns qualifier.
  // Note: CMtzDataFile is handled earlier in selectMtzColumns() - it skips the dialog entirely.
  const signatures =
    item._class === "CGenericReflDataFile"
      ? [...GENERIC_SIGNATURES]
      : [...(item._qualifiers?.correctColumns || [])];

  const options: ColumnOptions = {};

  signatures.forEach((signature) => {
    const signatureOptions = createSignatureOptions(
      signature,
      sortedColumnNames
    );
    if (signatureOptions.length > 0) {
      options[signature] = signatureOptions;
    }
  });

  return options;
};

const getSignatureLabel = (signature: string): string => {
  return SIGNATURE_MAP[signature] || signature;
};

// Dialog Component (internal)
const MtzColumnDialogComponent: React.FC<MtzColumnDialogProps> = ({
  columnNames,
  item,
  onAccept,
  onCancel,
}) => {
  const [open, setOpen] = useState(true);
  const [selectedGroup, setSelectedGroup] = useState<string | null>(null);
  const [values, dispatch] = useReducer(valuesReducer, {});

  // Build column options on mount
  const columnOptions = useMemo(
    () => buildColumnOptions(columnNames, item),
    [columnNames, item]
  );

  // Initialize values and selected group
  useEffect(() => {
    if (Object.keys(columnOptions).length > 0) {
      const initialValues: ValuesState = {};
      Object.entries(columnOptions).forEach(([signature, signatureOptions]) => {
        initialValues[signature] = signatureOptions[0];
      });

      Object.entries(initialValues).forEach(([signature, value]) => {
        dispatch({ type: "SET_VALUE", signature, value });
      });

      setSelectedGroup(Object.keys(columnOptions)[0]);
    }
  }, [columnOptions]);

  const handleGroupChange = useCallback(
    (_event: SyntheticEvent<Element, Event>, newValue: string | null) => {
      setSelectedGroup(newValue);
    },
    []
  );

  const handleColumnChange = useCallback(
    (signature: string) =>
      (
        _event: SyntheticEvent<Element, Event>,
        value: string | null,
        _reason: AutocompleteChangeReason
      ) => {
        dispatch({
          type: "SET_VALUE",
          signature,
          value: value || "",
        });
      },
    []
  );

  const handleAccept = useCallback(() => {
    if (selectedGroup && values[selectedGroup]) {
      setOpen(false);
      onAccept(values[selectedGroup]);
    }
  }, [selectedGroup, values, onAccept]);

  const handleCancel = useCallback(() => {
    setOpen(false);
    onCancel();
  }, [onCancel]);

  const isAcceptDisabled = !selectedGroup || !values[selectedGroup];

  // If no options available, auto-cancel
  useEffect(() => {
    if (Object.keys(columnOptions).length === 0) {
      handleCancel();
    }
  }, [columnOptions, handleCancel]);

  return (
    <SimpleDialog
      open={open}
      onClose={handleCancel}
      slotProps={{
        paper: {
          sx: { width: "80%", maxWidth: "none" },
        },
      }}
    >
      <DialogTitle>{item._objectPath}</DialogTitle>

      <DialogContent>
        <RadioGroup value={selectedGroup} onChange={handleGroupChange}>
          {Object.entries(columnOptions).map(
            ([signature, options]) =>
              options.length > 0 && (
                <Stack key={signature} direction="row" spacing={2}>
                  <FormControlLabel
                    sx={{ minWidth: "20rem" }}
                    value={signature}
                    control={<Radio size="small" />}
                    label={getSignatureLabel(signature)}
                  />

                  <Autocomplete
                    options={options}
                    value={values[signature] || ""}
                    onChange={handleColumnChange(signature)}
                    renderInput={(params) => (
                      <TextField
                        sx={{ my: 2, minWidth: "20rem" }}
                        {...params}
                        label="Columns"
                      />
                    )}
                  />
                </Stack>
              )
          )}
        </RadioGroup>
      </DialogContent>

      <DialogActions>
        <Button disabled={isAcceptDisabled} onClick={handleAccept}>
          OK
        </Button>
        <Button onClick={handleCancel}>Cancel</Button>
      </DialogActions>
    </SimpleDialog>
  );
};

// Types for the public API (defined early for use in enhanced dialog)

/** Single column selection for one data type */
export interface ColumnSelection {
  /** The signature type (KMKM, GLGL, JQ, FQ, I) */
  signature: string;
  /** Column selector path (e.g., a path ending with [col1,col2,col3,col4]) */
  columnSelector: string;
  /** Content flag value (1-4 for obs data, 1 for FreeR) */
  contentFlag: number;
  /** File suffix for multi-file output (e.g., "_asIPAIR") */
  fileSuffix: string;
  /** Whether this is the primary selection (lowest contentFlag) */
  isPrimary: boolean;
}

/** Result from enhanced MTZ column dialog */
export interface MtzColumnResult {
  /** Primary column selector (for backward compatibility) */
  columnSelector: string;
  /** All reflection data selections (may include multiple representations) */
  reflectionSelections: ColumnSelection[];
  /** FreeR column selection if available and task needs it */
  freeRSelection?: ColumnSelection;
}

/** Legacy result type for backward compatibility */
export interface MtzColumnResultLegacy {
  columnSelector: string;
}

/** Sibling input info for context-aware dialog */
export interface SiblingInput {
  /** Object path of the sibling input */
  objectPath: string;
  /** Class name (e.g., "CFreeRDataFile", "CObsDataFile") */
  className: string;
}

// Enhanced Dialog Types
interface EnhancedMtzColumnDialogProps {
  columnNames: ColumnNames;
  item: MtzItem;
  /** Whether to show multi-selection mode (checkboxes) or single-selection (radios) */
  multiSelectMode: boolean;
  /** Whether to show FreeR column section */
  showFreeR: boolean;
  onAccept: (result: MtzColumnResult) => void;
  onCancel: () => void;
}

interface SelectionState {
  /** For single-select mode: which signature is selected */
  primarySignature: string | null;
  /** For multi-select mode: which signatures are checked */
  checkedSignatures: Set<string>;
  /** Column selector values for each signature */
  columnValues: ValuesState;
  /** FreeR column selector if applicable */
  freeRValue: string | null;
  /** Whether FreeR is checked (in multi-select mode) */
  freeRChecked: boolean;
}

type SelectionAction =
  | { type: "SET_PRIMARY"; signature: string }
  | { type: "TOGGLE_CHECKED"; signature: string }
  | { type: "SET_COLUMN_VALUE"; signature: string; value: string }
  | { type: "SET_FREER_VALUE"; value: string }
  | { type: "TOGGLE_FREER" }
  | { type: "INIT"; primarySignature: string; columnValues: ValuesState; freeRValue: string | null };

const selectionReducer = (state: SelectionState, action: SelectionAction): SelectionState => {
  switch (action.type) {
    case "SET_PRIMARY":
      return { ...state, primarySignature: action.signature };
    case "TOGGLE_CHECKED": {
      const newChecked = new Set(state.checkedSignatures);
      if (newChecked.has(action.signature)) {
        newChecked.delete(action.signature);
      } else {
        newChecked.add(action.signature);
      }
      return { ...state, checkedSignatures: newChecked };
    }
    case "SET_COLUMN_VALUE":
      return {
        ...state,
        columnValues: { ...state.columnValues, [action.signature]: action.value },
      };
    case "SET_FREER_VALUE":
      return { ...state, freeRValue: action.value };
    case "TOGGLE_FREER":
      return { ...state, freeRChecked: !state.freeRChecked };
    case "INIT":
      return {
        ...state,
        primarySignature: action.primarySignature,
        columnValues: action.columnValues,
        freeRValue: action.freeRValue,
        checkedSignatures: new Set([action.primarySignature]),
        freeRChecked: action.freeRValue !== null,
      };
    default:
      return state;
  }
};

/** Build FreeR column options from parsed column names */
const buildFreeROptions = (allColumnNames: ColumnNames): string[] => {
  const sortedColumnNames: Record<string, string[]> = {};

  Object.entries(allColumnNames).forEach(([label, columnType]) => {
    if (!sortedColumnNames[columnType]) {
      sortedColumnNames[columnType] = [];
    }
    sortedColumnNames[columnType].push(label);
  });

  // FreeR uses signature "I" (single integer column)
  return createSignatureOptions("I", sortedColumnNames);
};

/** Create ColumnSelection object from signature and column selector */
const createColumnSelection = (
  signature: string,
  columnSelector: string,
  isPrimary: boolean
): ColumnSelection => ({
  signature,
  columnSelector,
  contentFlag: SIGNATURE_TO_CONTENT_FLAG[signature] ?? 1,
  fileSuffix: SIGNATURE_TO_SUFFIX[signature] ?? "",
  isPrimary,
});

// Enhanced Dialog Component with multi-selection support
const EnhancedMtzColumnDialogComponent: React.FC<EnhancedMtzColumnDialogProps> = ({
  columnNames,
  item,
  multiSelectMode,
  showFreeR,
  onAccept,
  onCancel,
}) => {
  const [open, setOpen] = useState(true);

  const initialState: SelectionState = {
    primarySignature: null,
    checkedSignatures: new Set(),
    columnValues: {},
    freeRValue: null,
    freeRChecked: false,
  };
  const [state, dispatch] = useReducer(selectionReducer, initialState);

  // Build column options for reflections
  const columnOptions = useMemo(
    () => buildColumnOptions(columnNames, item),
    [columnNames, item]
  );

  // Build FreeR options if needed
  const freeROptions = useMemo(
    () => (showFreeR ? buildFreeROptions(columnNames) : []),
    [columnNames, showFreeR]
  );

  // Initialize state
  useEffect(() => {
    const signatures = Object.keys(columnOptions);
    if (signatures.length > 0) {
      const initialValues: ValuesState = {};
      signatures.forEach((sig) => {
        initialValues[sig] = columnOptions[sig][0];
      });

      const freeRInitial = freeROptions.length > 0 ? freeROptions[0] : null;

      dispatch({
        type: "INIT",
        primarySignature: signatures[0],
        columnValues: initialValues,
        freeRValue: freeRInitial,
      });
    }
  }, [columnOptions, freeROptions]);

  const handlePrimaryChange = useCallback(
    (_event: SyntheticEvent<Element, Event>, newValue: string | null) => {
      if (newValue) {
        dispatch({ type: "SET_PRIMARY", signature: newValue });
      }
    },
    []
  );

  const handleCheckboxChange = useCallback(
    (signature: string) => () => {
      dispatch({ type: "TOGGLE_CHECKED", signature });
    },
    []
  );

  const handleColumnChange = useCallback(
    (signature: string) =>
      (
        _event: SyntheticEvent<Element, Event>,
        value: string | null,
        _reason: AutocompleteChangeReason
      ) => {
        dispatch({ type: "SET_COLUMN_VALUE", signature, value: value || "" });
      },
    []
  );

  const handleFreeRChange = useCallback(
    (
      _event: SyntheticEvent<Element, Event>,
      value: string | null,
      _reason: AutocompleteChangeReason
    ) => {
      dispatch({ type: "SET_FREER_VALUE", value: value || "" });
    },
    []
  );

  const handleFreeRCheckboxChange = useCallback(() => {
    dispatch({ type: "TOGGLE_FREER" });
  }, []);

  const handleAccept = useCallback(() => {
    if (!state.primarySignature || !state.columnValues[state.primarySignature]) {
      return;
    }

    const reflectionSelections: ColumnSelection[] = [];

    if (multiSelectMode) {
      // Multi-select mode: collect all checked signatures
      // Sort by content flag to determine primary (lowest flag = primary)
      const checkedSigs = Array.from(state.checkedSignatures).sort(
        (a, b) => (SIGNATURE_TO_CONTENT_FLAG[a] ?? 99) - (SIGNATURE_TO_CONTENT_FLAG[b] ?? 99)
      );

      checkedSigs.forEach((sig, index) => {
        if (state.columnValues[sig]) {
          reflectionSelections.push(
            createColumnSelection(sig, state.columnValues[sig], index === 0)
          );
        }
      });
    } else {
      // Single-select mode: only primary signature
      reflectionSelections.push(
        createColumnSelection(
          state.primarySignature,
          state.columnValues[state.primarySignature],
          true
        )
      );
    }

    // Build result
    const result: MtzColumnResult = {
      columnSelector: state.columnValues[state.primarySignature],
      reflectionSelections,
    };

    // Add FreeR if checked and available
    if (showFreeR && state.freeRChecked && state.freeRValue) {
      result.freeRSelection = {
        signature: "I",
        columnSelector: state.freeRValue,
        contentFlag: 1, // FREER content flag
        fileSuffix: "_asFREER",
        isPrimary: false,
      };
    }

    setOpen(false);
    onAccept(result);
  }, [state, multiSelectMode, showFreeR, onAccept]);

  const handleCancel = useCallback(() => {
    setOpen(false);
    onCancel();
  }, [onCancel]);

  // Determine if accept is disabled
  const isAcceptDisabled = useMemo(() => {
    if (multiSelectMode) {
      // At least one signature must be checked with a valid value
      return !Array.from(state.checkedSignatures).some(
        (sig) => state.columnValues[sig]
      );
    } else {
      return !state.primarySignature || !state.columnValues[state.primarySignature];
    }
  }, [multiSelectMode, state]);

  // Auto-cancel if no options
  useEffect(() => {
    if (Object.keys(columnOptions).length === 0) {
      handleCancel();
    }
  }, [columnOptions, handleCancel]);

  return (
    <SimpleDialog
      open={open}
      onClose={handleCancel}
      slotProps={{
        paper: {
          sx: { width: "80%", maxWidth: "none" },
        },
      }}
    >
      <DialogTitle>{item._objectPath}</DialogTitle>

      <DialogContent>
        {/* Reflection columns section */}
        <Typography variant="subtitle2" sx={{ mb: 1, mt: 1 }}>
          Reflection Data
        </Typography>

        {multiSelectMode ? (
          // Multi-select mode with checkboxes
          <Stack spacing={1}>
            {Object.entries(columnOptions).map(
              ([signature, options]) =>
                options.length > 0 && (
                  <Stack key={signature} direction="row" spacing={2} alignItems="center">
                    <FormControlLabel
                      sx={{ minWidth: "20rem" }}
                      control={
                        <Checkbox
                          size="small"
                          checked={state.checkedSignatures.has(signature)}
                          onChange={handleCheckboxChange(signature)}
                        />
                      }
                      label={getSignatureLabel(signature)}
                    />
                    <Autocomplete
                      options={options}
                      value={state.columnValues[signature] || ""}
                      onChange={handleColumnChange(signature)}
                      disabled={!state.checkedSignatures.has(signature)}
                      sx={{ flexGrow: 1 }}
                      renderInput={(params) => (
                        <TextField
                          sx={{ my: 1, minWidth: "20rem" }}
                          {...params}
                          label="Columns"
                          size="small"
                        />
                      )}
                    />
                  </Stack>
                )
            )}
          </Stack>
        ) : (
          // Single-select mode with radios
          <RadioGroup value={state.primarySignature} onChange={handlePrimaryChange}>
            {Object.entries(columnOptions).map(
              ([signature, options]) =>
                options.length > 0 && (
                  <Stack key={signature} direction="row" spacing={2}>
                    <FormControlLabel
                      sx={{ minWidth: "20rem" }}
                      value={signature}
                      control={<Radio size="small" />}
                      label={getSignatureLabel(signature)}
                    />
                    <Autocomplete
                      options={options}
                      value={state.columnValues[signature] || ""}
                      onChange={handleColumnChange(signature)}
                      renderInput={(params) => (
                        <TextField
                          sx={{ my: 2, minWidth: "20rem" }}
                          {...params}
                          label="Columns"
                        />
                      )}
                    />
                  </Stack>
                )
            )}
          </RadioGroup>
        )}

        {/* FreeR section if applicable */}
        {showFreeR && freeROptions.length > 0 && (
          <>
            <Divider sx={{ my: 2 }} />
            <Typography variant="subtitle2" sx={{ mb: 1 }}>
              FreeR Flags (Optional)
            </Typography>
            <Stack direction="row" spacing={2} alignItems="center">
              <FormControlLabel
                sx={{ minWidth: "20rem" }}
                control={
                  <Checkbox
                    size="small"
                    checked={state.freeRChecked}
                    onChange={handleFreeRCheckboxChange}
                  />
                }
                label={getSignatureLabel("I")}
              />
              <Autocomplete
                options={freeROptions}
                value={state.freeRValue || ""}
                onChange={handleFreeRChange}
                disabled={!state.freeRChecked}
                sx={{ flexGrow: 1 }}
                renderInput={(params) => (
                  <TextField
                    sx={{ my: 1, minWidth: "20rem" }}
                    {...params}
                    label="FreeR Column"
                    size="small"
                  />
                )}
              />
            </Stack>
          </>
        )}
      </DialogContent>

      <DialogActions>
        <Button disabled={isAcceptDisabled} onClick={handleAccept}>
          OK
        </Button>
        <Button onClick={handleCancel}>Cancel</Button>
      </DialogActions>
    </SimpleDialog>
  );
};

export interface ParseMtzOptions {
  file: File;
  item: MtzItem;
  cootModule: any; // The coot WASM module (optional - no longer used for column parsing)
}

/**
 * Parse an MTZ file and get column names using the pure TypeScript parser.
 *
 * This is the primary parser - it's faster than WASM and works reliably
 * across all environments (including Azure where SharedArrayBuffer isn't
 * available for WASM threading).
 */
async function parseMtzColumnsNative(file: File): Promise<ColumnNames | null> {
  try {
    const { parseMtzFile, getColumnNamesByType } = await import("../../../lib/mtz-parser");
    const header = await parseMtzFile(file);
    const columns = getColumnNamesByType(header);

    if (Object.keys(columns).length === 0) {
      return null;
    }

    return columns;
  } catch (error) {
    console.error("Failed to parse MTZ file:", error);
    return null;
  }
}

/**
 * Parse an MTZ file and get column names.
 *
 * Uses a pure TypeScript implementation that parses the MTZ header directly.
 * This is faster and more reliable than the previous WASM-based approach,
 * and works in all environments without SharedArrayBuffer requirements.
 *
 * @param file - The MTZ file to parse
 * @param _cootModule - Deprecated, no longer used (kept for API compatibility)
 * @returns Column names map or null if parsing fails
 */
export async function parseMtzColumns(
  file: File,
  _cootModule?: any | null
): Promise<ColumnNames | null> {
  return parseMtzColumnsNative(file);
}

/**
 * Show the MTZ column selection dialog.
 * Returns a Promise that resolves with the selected column selector string,
 * or null if the user cancels.
 */
export function showMtzColumnDialog(
  columnNames: ColumnNames,
  item: MtzItem
): Promise<string | null> {
  return new Promise((resolve) => {
    // Create a container for the dialog
    const container = document.createElement("div");
    container.id = `mtz-dialog-${Date.now()}`;
    document.body.appendChild(container);

    const root = createRoot(container);

    const cleanup = () => {
      // Small delay to allow dialog close animation
      setTimeout(() => {
        root.unmount();
        container.remove();
      }, 300);
    };

    const handleAccept = (columnSelector: string) => {
      cleanup();
      resolve(columnSelector);
    };

    const handleCancel = () => {
      cleanup();
      resolve(null);
    };

    root.render(
      <MtzColumnDialogComponent
        columnNames={columnNames}
        item={item}
        onAccept={handleAccept}
        onCancel={handleCancel}
      />
    );
  });
}

/**
 * Complete flow: parse MTZ file and show column selection dialog if needed.
 * Returns the column selector string, or null if cancelled/failed.
 *
 * For non-MTZ files, returns a default column selector.
 * For CMtzDataFile (general container), skips column selection entirely.
 */
export async function selectMtzColumns(
  options: ParseMtzOptions
): Promise<string | null> {
  const { file, item, cootModule } = options;

  // Check if it's an MTZ file
  const isMtzFile = file.name.toLowerCase().endsWith(".mtz");

  if (!isMtzFile) {
    // For non-MTZ files, return a default
    return "/*/*/[FP,SIGFP]";
  }

  // CMtzDataFile is a general container that stores the intact reflection file
  // without separating out columns - no column selection needed
  if (item._class === "CMtzDataFile") {
    return "";  // Empty string signals "store whole file as-is"
  }

  // Parse the MTZ file
  const columnNames = await parseMtzColumns(file, cootModule);

  if (!columnNames) {
    // Parsing failed
    return null;
  }

  // Show the dialog and get user selection
  return showMtzColumnDialog(columnNames, item);
}

// Enhanced API with multi-selection support

/** Enhanced options for MTZ column selection with sibling awareness */
export interface ParseMtzOptionsEnhanced {
  file: File;
  item: MtzItem;
  cootModule: any;
  /** Sibling inputs from the same task (to detect if FreeR input exists) */
  siblingInputs?: SiblingInput[];
  /** Enable multi-selection mode for multiple reflection representations */
  multiSelectMode?: boolean;
}

/**
 * Show the enhanced MTZ column selection dialog with multi-selection support.
 * Returns a Promise that resolves with MtzColumnResult or null if cancelled.
 *
 * @param columnNames - Parsed column names from MTZ file
 * @param item - The MTZ item being configured
 * @param options - Dialog options including multi-select mode and FreeR display
 */
export function showMtzColumnDialogEnhanced(
  columnNames: ColumnNames,
  item: MtzItem,
  options: {
    multiSelectMode?: boolean;
    showFreeR?: boolean;
  } = {}
): Promise<MtzColumnResult | null> {
  const { multiSelectMode = false, showFreeR = false } = options;

  return new Promise((resolve) => {
    const container = document.createElement("div");
    container.id = `mtz-dialog-enhanced-${Date.now()}`;
    document.body.appendChild(container);

    const root = createRoot(container);

    const cleanup = () => {
      setTimeout(() => {
        root.unmount();
        container.remove();
      }, 300);
    };

    const handleAccept = (result: MtzColumnResult) => {
      cleanup();
      resolve(result);
    };

    const handleCancel = () => {
      cleanup();
      resolve(null);
    };

    root.render(
      <EnhancedMtzColumnDialogComponent
        columnNames={columnNames}
        item={item}
        multiSelectMode={multiSelectMode}
        showFreeR={showFreeR}
        onAccept={handleAccept}
        onCancel={handleCancel}
      />
    );
  });
}

/**
 * Enhanced MTZ column selection flow with multi-selection and FreeR support.
 *
 * Determines whether to show the enhanced dialog based on:
 * - multiSelectMode: If true, allows selecting multiple reflection representations
 * - siblingInputs: If a CFreeRDataFile sibling exists, shows FreeR column selection
 *
 * Returns MtzColumnResult with potentially multiple selections, or null if cancelled.
 */
export async function selectMtzColumnsEnhanced(
  options: ParseMtzOptionsEnhanced
): Promise<MtzColumnResult | null> {
  const {
    file,
    item,
    cootModule,
    siblingInputs = [],
    multiSelectMode = false,
  } = options;

  // Check if it's an MTZ file
  const isMtzFile = file.name.toLowerCase().endsWith(".mtz");

  if (!isMtzFile) {
    // For non-MTZ files, return a simple default result
    return {
      columnSelector: "/*/*/[FP,SIGFP]",
      reflectionSelections: [
        {
          signature: "FQ",
          columnSelector: "/*/*/[FP,SIGFP]",
          contentFlag: 4,
          fileSuffix: "",
          isPrimary: true,
        },
      ],
    };
  }

  // CMtzDataFile stores the intact file without column separation
  if (item._class === "CMtzDataFile") {
    return {
      columnSelector: "",
      reflectionSelections: [],
    };
  }

  // Parse the MTZ file
  const columnNames = await parseMtzColumns(file, cootModule);

  if (!columnNames) {
    return null;
  }

  // Determine if we should show FreeR section
  // Show FreeR if there's a sibling CFreeRDataFile input
  const hasFreeRSibling = siblingInputs.some(
    (sibling) => sibling.className === "CFreeRDataFile"
  );

  // Show the enhanced dialog
  return showMtzColumnDialogEnhanced(columnNames, item, {
    multiSelectMode,
    showFreeR: hasFreeRSibling,
  });
}
