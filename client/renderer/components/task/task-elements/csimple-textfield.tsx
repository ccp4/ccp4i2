import React, {
  ChangeEvent,
  KeyboardEvent,
  useCallback,
  useEffect,
  useMemo,
  useRef,
} from "react";
import { Checkbox, FormControlLabel, TextField } from "@mui/material";

import { CCP4i2CSimpleElementProps } from "./csimple";
import { ErrorTrigger } from "./error-info";
import {
  FULL_WIDTH_FIELD_STYLES,
  getContainerSizeStyles,
} from "./field-sizes";
import { FieldWrapper } from "./field-wrapper";
import {
  useContainerField,
  useSyncedLocalValue,
} from "./hooks/useContainerField";

type InputValue = string | number | boolean;
type InputType = "text" | "int" | "float" | "checkbox";

const DEBOUNCE_DELAY = 1000;

const isValueValid = (value: InputValue, type: InputType): boolean => {
  if (type === "int") {
    return typeof value === "string"
      ? /^\d+$/.test(value)
      : Number.isInteger(value);
  }
  if (type === "float") {
    return typeof value === "string"
      ? /^-?\d*\.?\d*$/.test(value)
      : typeof value === "number" && !Number.isNaN(value);
  }
  return true;
};

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
  const debounceRef = useRef<NodeJS.Timeout | null>(null);

  const {
    item,
    objectPath,
    serverValue,
    isVisible,
    isDisabled,
    validationColor,
    commit,
  } = useContainerField<InputValue>({
    job,
    itemName,
    visibility,
    disabled: disabledProp,
    suppressMutations,
    onChange,
  });

  const guiLabel =
    qualifiers?.guiLabel || objectPath?.split(".").at(-1) || "";
  const isMultiLine = qualifiers?.guiMode === "multiLine";
  const tooltipText = qualifiers?.toolTip || objectPath || "";

  // An unset primitive arrives with a sentinel _value (0 for numbers), so it
  // renders blank rather than as a typed 0; the qualifier default, when there
  // is one, is shown as placeholder text so the user can see what blank means.
  const isUnset = item?._valueState === "NOT_SET";
  const displayedServerValue: InputValue = isUnset ? "" : serverValue ?? "";
  const defaultQualifier = qualifiers?.default;
  const placeholder =
    defaultQualifier !== undefined && defaultQualifier !== null
      ? String(defaultQualifier)
      : undefined;

  const [value, setValue] = useSyncedLocalValue<InputValue>(
    displayedServerValue
  );

  useEffect(() => {
    return () => {
      if (debounceRef.current) clearTimeout(debounceRef.current);
    };
  }, []);

  const calculatedSx = useMemo(
    () => ({
      ...(sizeProp ? getContainerSizeStyles(sizeProp) : FULL_WIDTH_FIELD_STYLES),
      ...sx,
    }),
    [sizeProp, sx]
  );

  const hasError =
    validationColor === "error.light" ||
    !isValueValid(value, type as InputType);

  const slotProps = useMemo(
    () => ({ inputLabel: { shrink: true, disableAnimation: true } }),
    []
  );

  const handleCommit = useCallback(
    async (newValue: InputValue) => {
      let parsed: InputValue | null = newValue;
      if ((type === "int" || type === "float") && typeof newValue === "string") {
        const trimmed = newValue.trim();
        if (trimmed === "") {
          // Clearing a numeric field returns the parameter to unset (the
          // server maps null to unSet()), so the wrapper's default applies.
          if (isUnset) return;
          parsed = null;
        } else if (type === "int" && /^\d+$/.test(trimmed)) {
          parsed = parseInt(trimmed, 10);
        } else if (type === "float" && /^-?\d*\.?\d+$/.test(trimmed)) {
          parsed = parseFloat(trimmed);
        }
      }

      const result = await commit(parsed);
      if (result && !result.success) {
        setValue(displayedServerValue);
      }
    },
    [commit, type, isUnset, displayedServerValue, setValue]
  );

  const handleChange = useCallback(
    (event: ChangeEvent<HTMLInputElement | HTMLTextAreaElement>) => {
      if (type === "checkbox") {
        const newValue = (event.target as HTMLInputElement).checked;
        setValue(newValue);
        handleCommit(newValue);
        return;
      }

      const inputValue = event.target.value;
      setValue(inputValue);

      if (debounceRef.current) clearTimeout(debounceRef.current);
      debounceRef.current = setTimeout(() => {
        handleCommit(inputValue);
        debounceRef.current = null;
      }, DEBOUNCE_DELAY);
    },
    [type, handleCommit, setValue]
  );

  const handleKeyDown = useCallback(
    (event: KeyboardEvent<HTMLDivElement>) => {
      if (event.key === "Enter" && value !== displayedServerValue) {
        handleCommit(value);
      }
    },
    [value, displayedServerValue, handleCommit]
  );

  const handleBlur = useCallback(() => {
    if (value !== displayedServerValue) handleCommit(value);
  }, [value, displayedServerValue, handleCommit]);

  if (!isVisible) return null;

  if (type === "checkbox") {
    return (
      <FieldWrapper ariaLabel={`${guiLabel} input`}>
        <FormControlLabel
          disabled={isDisabled}
          title={tooltipText}
          control={
            <Checkbox
              inputRef={inputRef}
              checked={Boolean(value)}
              onChange={handleChange}
              size="small"
            />
          }
          label={guiLabel}
        />
        <ErrorTrigger item={item} job={job} />
      </FieldWrapper>
    );
  }

  return (
    <FieldWrapper ariaLabel={`${guiLabel} input`}>
      <TextField
        inputRef={inputRef}
        disabled={isDisabled}
        multiline={isMultiLine}
        size="small"
        sx={calculatedSx}
        slotProps={slotProps}
        type="text"
        value={value}
        placeholder={placeholder}
        label={guiLabel}
        title={tooltipText}
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
