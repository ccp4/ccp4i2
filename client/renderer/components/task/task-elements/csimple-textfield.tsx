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
  getFieldSizeStyles,
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

  const [value, setValue] = useSyncedLocalValue<InputValue>(serverValue ?? "");

  useEffect(() => {
    return () => {
      if (debounceRef.current) clearTimeout(debounceRef.current);
    };
  }, []);

  const calculatedSx = useMemo(
    () => ({
      ...(sizeProp ? getFieldSizeStyles(sizeProp) : FULL_WIDTH_FIELD_STYLES),
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
      let parsed: InputValue = newValue;
      if (type === "int" && typeof newValue === "string") {
        if (newValue.trim() === "") return;
        if (/^\d+$/.test(newValue)) parsed = parseInt(newValue, 10);
      }
      if (type === "float" && typeof newValue === "string") {
        if (newValue.trim() === "") return;
        if (/^-?\d*\.?\d+$/.test(newValue)) parsed = parseFloat(newValue);
      }

      const result = await commit(parsed);
      if (result && !result.success) {
        setValue(serverValue ?? "");
      }
    },
    [commit, type, serverValue, setValue]
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
      if (event.key === "Enter" && value !== serverValue) {
        handleCommit(value);
      }
    },
    [value, serverValue, handleCommit]
  );

  const handleBlur = useCallback(() => {
    if (value !== serverValue) handleCommit(value);
  }, [value, serverValue, handleCommit]);

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
