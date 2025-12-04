import React, {
  ChangeEvent,
  KeyboardEvent,
  useCallback,
  useEffect,
  useMemo,
  useRef,
  useState,
} from "react";
import { TextField } from "@mui/material";

import { CCP4i2CSimpleElementProps } from "./csimple";
import { useJob, SetParameterResponse } from "../../../utils";
import { ErrorTrigger } from "./error-info";
import { useTaskInterface } from "../../../providers/task-provider";
import { usePopcorn } from "../../../providers/popcorn-provider";
import { useParameterChangeIntent } from "../../../providers/parameter-change-intent-provider";
import { inferFieldSize, getFieldSizeStyles } from "./field-sizes";
import { FieldWrapper } from "./field-wrapper";

// Types
type InputValue = string | number | boolean;
type InputType = "text" | "int" | "float" | "checkbox";

interface ProcessedItem {
  objectPath: string | null;
  value: InputValue;
  guiLabel: string;
  isMultiLine: boolean;
  tooltipText: string;
}

// Constants
const DEBOUNCE_DELAY = 1000;
const INPUT_TYPES = {
  TEXT: "text",
  INT: "int",
  FLOAT: "float",
  CHECKBOX: "checkbox",
} as const;

// Custom hooks
const useProcessedItem = (
  item: any,
  qualifiers: any,
  objectPath: string | null
): ProcessedItem => {
  return useMemo(() => {
    const guiLabel =
      qualifiers?.guiLabel || objectPath?.split(".").at(-1) || "";
    const isMultiLine = qualifiers?.guiMode === "multiLine";
    const tooltipText = qualifiers?.toolTip || objectPath || "";

    return {
      objectPath,
      value: item?._value ?? "",
      guiLabel,
      isMultiLine,
      tooltipText,
    };
  }, [item, qualifiers, objectPath]);
};

/**
 * Custom hook for form state with intent-aware syncing.
 *
 * Key improvement: Only syncs from server (initialValue) if the user
 * hasn't recently edited this field. This prevents the jarring experience
 * of having your edit overwritten by a container refetch.
 */
const useFormState = (
  initialValue: InputValue,
  type: InputType,
  objectPath: string | null,
  wasRecentlyChanged: (path: string, withinMs?: number) => boolean
) => {
  const [value, setValue] = useState<InputValue>(initialValue);
  const [isSubmitting, setIsSubmitting] = useState(false);
  const debounceRef = useRef<NodeJS.Timeout | null>(null);

  // Sync with prop changes - but skip if user recently edited this field
  useEffect(() => {
    // Don't overwrite local state if user just edited this field
    // Use a longer timeout (10 seconds) to account for API latency,
    // container mutation, and SWR revalidation
    if (objectPath && wasRecentlyChanged(objectPath, 10000)) {
      return;
    }
    setValue(initialValue);
  }, [initialValue, objectPath, wasRecentlyChanged]);

  // Cleanup debounce on unmount
  useEffect(() => {
    return () => {
      if (debounceRef.current) {
        clearTimeout(debounceRef.current);
      }
    };
  }, []);

  const setDebouncedValue = useCallback(
    (newValue: InputValue, callback: (value: InputValue) => void) => {
      if (debounceRef.current) {
        clearTimeout(debounceRef.current);
      }

      debounceRef.current = setTimeout(() => {
        callback(newValue);
        debounceRef.current = null;
      }, DEBOUNCE_DELAY);
    },
    []
  );

  return {
    value,
    setValue,
    isSubmitting,
    setIsSubmitting,
    setDebouncedValue,
  };
};

// Utility functions
const parseValueByType = (inputValue: string, type: InputType): InputValue => {
  switch (type) {
    case INPUT_TYPES.INT:
      // Only parse if it's a valid integer string
      return /^\d+$/.test(inputValue) ? parseInt(inputValue, 10) : inputValue;
    case INPUT_TYPES.FLOAT:
      // Only parse if it's a valid float string
      return /^-?\d*\.?\d*$/.test(inputValue) ? inputValue : inputValue;
    case INPUT_TYPES.TEXT:
    default:
      return inputValue;
  }
};

const isValueValid = (value: InputValue, type: InputType): boolean => {
  if (type === INPUT_TYPES.INT) {
    return typeof value === "string"
      ? /^\d+$/.test(value)
      : Number.isInteger(value);
  }
  if (type === INPUT_TYPES.FLOAT) {
    return typeof value === "string"
      ? /^-?\d*\.?\d*$/.test(value)
      : typeof value === "number" && !Number.isNaN(value);
  }
  return true;
};

// Main component
export const CSimpleTextFieldElement: React.FC<CCP4i2CSimpleElementProps> = ({
  itemName,
  job,
  type,
  sx,
  qualifiers,
  onChange,
  visibility,
  disabled: disabledProp,
  suppressMutations = false,
  size: sizeProp,
}) => {
  const inputRef = useRef<HTMLInputElement | null>(null);

  const {
    useTaskItem,
    getValidationColor,
    setParameter,
    setParameterNoMutate,
  } = useJob(job.id);
  const { item } = useTaskItem(itemName);
  const { inFlight, setInFlight } = useTaskInterface();
  const { setMessage } = usePopcorn();
  const { setIntentForPath, clearIntentForPath, wasRecentlyChanged } =
    useParameterChangeIntent();

  // Process item data
  const objectPath = useMemo(() => item?._objectPath || null, [item]);
  const processedItem = useProcessedItem(item, qualifiers, objectPath);

  // Use intent-aware form state
  const { value, setValue, isSubmitting, setIsSubmitting, setDebouncedValue } =
    useFormState(
      processedItem.value,
      type as InputType,
      objectPath,
      wasRecentlyChanged
    );

  // Computed properties
  const isVisible = useMemo(() => {
    if (typeof visibility === "function") return visibility();
    return visibility !== false;
  }, [visibility]);

  const isDisabled = useMemo(() => {
    if (typeof disabledProp === "function") {
      return disabledProp() || inFlight || isSubmitting || job.status !== 1;
    }
    return disabledProp || inFlight || isSubmitting || job.status !== 1;
  }, [disabledProp, inFlight, isSubmitting, job.status]);

  // Calculate field size - use explicit prop or infer from item/qualifiers
  const fieldSize = useMemo(
    () => sizeProp || inferFieldSize(item, qualifiers),
    [sizeProp, item, qualifiers]
  );

  const calculatedSx = useMemo(
    () => ({
      ...getFieldSizeStyles(fieldSize),
      ...sx,
    }),
    [fieldSize, sx]
  );

  const validationColor = useMemo(
    () => getValidationColor(item),
    [getValidationColor, item]
  );

  const hasError = useMemo(
    () =>
      validationColor === "error.light" ||
      !isValueValid(value, type as InputType),
    [validationColor, value, type]
  );

  const slotProps = useMemo(() => {
    const baseProps = {
      inputLabel: {
        shrink: true,
        disableAnimation: true,
      },
    };

    if (type === INPUT_TYPES.CHECKBOX) {
      return {
        ...baseProps,
        htmlInput: {
          checked: Boolean(value),
          sx: { my: 1 },
          "aria-label": processedItem.guiLabel,
        },
      };
    }

    return baseProps;
  }, [type, value, processedItem.guiLabel]);

  // Event handlers
  const handleParameterUpdate = useCallback(
    async (newValue: InputValue) => {
      if (!objectPath) {
        console.error("No object path available for parameter update");
        return;
      }

      let parsedValue = newValue;
      if (
        type === INPUT_TYPES.INT &&
        typeof newValue === "string" &&
        /^\d+$/.test(newValue)
      ) {
        parsedValue = parseInt(newValue, 10);
      }
      if (
        type === INPUT_TYPES.FLOAT &&
        typeof newValue === "string" &&
        /^-?\d*\.?\d+$/.test(newValue)
      ) {
        parsedValue = parseFloat(newValue);
      }

      const setParameterArg = {
        object_path: objectPath,
        value: parsedValue,
      };

      // Record intent BEFORE making the API call
      // This prevents the container refetch from overwriting our local state
      setIntentForPath({
        jobId: job.id,
        parameterPath: objectPath,
        reason: "UserEdit",
        previousValue: item?._value,
      });

      setInFlight(true);
      setIsSubmitting(true);

      try {
        const updateFn = suppressMutations
          ? setParameterNoMutate
          : setParameter;
        const result: SetParameterResponse | undefined =
          await updateFn(setParameterArg);

        if (result && !result.success) {
          setMessage(`Unacceptable value provided: "${newValue}"`);
          setValue(item?._value ?? ""); // Revert to original value
          // Clear intent on failure so next sync can happen
          clearIntentForPath(objectPath);
        } else if (result?.success && result.data?.updated_item && onChange) {
          // Clear intent after successful update and before calling onChange
          // This allows derived updates to trigger proper re-renders
          clearIntentForPath(objectPath);
          await onChange(result.data.updated_item);
        } else if (result?.success) {
          // Successful but no onChange handler - still clear intent
          clearIntentForPath(objectPath);
        }
      } catch (error) {
        const errorMessage =
          error instanceof Error ? error.message : String(error);
        setMessage(`Error updating parameter: ${errorMessage}`);
        console.error("Parameter update failed:", error);
        setValue(item?._value ?? ""); // Revert to original value
        // Clear intent on error so next sync can happen
        clearIntentForPath(objectPath);
      } finally {
        setInFlight(false);
        setIsSubmitting(false);
      }
    },
    [
      objectPath,
      suppressMutations,
      setParameterNoMutate,
      setParameter,
      setInFlight,
      setIsSubmitting,
      setMessage,
      onChange,
      item,
      type,
      job.id,
      setIntentForPath,
      clearIntentForPath,
    ]
  );

  const handleChange = useCallback(
    (event: ChangeEvent<HTMLInputElement | HTMLTextAreaElement>) => {
      if (type === INPUT_TYPES.CHECKBOX) {
        const newValue = (event.target as HTMLInputElement).checked;
        setValue(newValue);
        // Immediately update checkbox values
        handleParameterUpdate(newValue);
      } else {
        const inputValue = event.target.value;
        setValue(inputValue); // Always store as string while editing

        // Record intent immediately when user starts typing
        // This prevents any in-flight container refetch from overwriting
        if (objectPath) {
          setIntentForPath({
            jobId: job.id,
            parameterPath: objectPath,
            reason: "UserEdit",
            previousValue: item?._value,
          });
        }

        // Debounce updates for text inputs
        setDebouncedValue(inputValue, handleParameterUpdate);
      }
    },
    [
      type,
      handleParameterUpdate,
      setDebouncedValue,
      objectPath,
      job.id,
      item?._value,
      setIntentForPath,
    ]
  );

  const handleKeyDown = useCallback(
    (event: KeyboardEvent<HTMLDivElement>) => {
      if (event.key === "Enter" && value !== item?._value) {
        handleParameterUpdate(value);
      }
    },
    [value, item, handleParameterUpdate]
  );

  const handleBlur = useCallback(() => {
    if (value !== item?._value) {
      handleParameterUpdate(value);
    }
  }, [value, item, handleParameterUpdate]);

  // Sync checkbox state with ref
  useEffect(() => {
    if (type === INPUT_TYPES.CHECKBOX && inputRef.current) {
      inputRef.current.checked = Boolean(item?._value);
    }
  }, [type, item]);

  // Early return for invisible components
  if (!isVisible) {
    return null;
  }

  return (
    <FieldWrapper ariaLabel={`${processedItem.guiLabel} input`}>
      <TextField
        inputRef={inputRef}
        disabled={isDisabled}
        multiline={processedItem.isMultiLine}
        size="small"
        sx={calculatedSx}
        slotProps={slotProps}
        type={type === INPUT_TYPES.CHECKBOX ? "checkbox" : "text"}
        value={type === INPUT_TYPES.CHECKBOX ? undefined : value}
        label={processedItem.guiLabel}
        title={processedItem.tooltipText}
        onChange={handleChange}
        onKeyDown={handleKeyDown}
        onBlur={handleBlur}
        error={hasError}
        inputProps={{
          "aria-describedby": hasError ? `${itemName}-error` : undefined,
          "aria-invalid": hasError,
        }}
      />
      <ErrorTrigger item={item} job={job} />
    </FieldWrapper>
  );
};
