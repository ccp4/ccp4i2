import React, {
  SyntheticEvent,
  useCallback,
  useMemo,
} from "react";
import {
  Autocomplete,
  AutocompleteChangeReason,
  FormControlLabel,
  Radio,
  RadioGroup,
  TextField,
} from "@mui/material";

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

type OptionValue = string | number;

export const CSimpleAutocompleteElement: React.FC<
  CCP4i2CSimpleElementProps
> = ({
  itemName,
  job,
  sx,
  qualifiers,
  onChange,
  visibility,
  disabled: disabledProp,
  suppressMutations = false,
  size: sizeProp,
}) => {
  const {
    item,
    objectPath,
    serverValue,
    isVisible,
    isDisabled,
    validationColor,
    commit,
  } = useContainerField<OptionValue>({
    job,
    itemName,
    visibility,
    disabled: disabledProp,
    suppressMutations,
    onChange,
  });

  const [localValue, setLocalValue] = useSyncedLocalValue<OptionValue>(
    serverValue ?? ""
  );

  const { enumerators, labels, guiLabel, guiMode, onlyEnumerators } = useMemo(() => {
    const rawEnumerators =
      qualifiers?.enumerators?.map((element: any) => {
        if (typeof element === "string" || element instanceof String) {
          return element.trim();
        }
        return element;
      }) || [];

    const processedEnumerators = [...rawEnumerators];
    if (serverValue && !processedEnumerators.includes(serverValue)) {
      processedEnumerators.push(serverValue);
    }

    let processedLabels = processedEnumerators.map((v) => `${v}`);
    if (
      qualifiers?.menuText &&
      qualifiers.menuText.length === processedEnumerators.length
    ) {
      processedLabels = qualifiers.menuText.map((text: string) => text.trim());
    }

    const processedGuiLabel =
      qualifiers?.guiLabel || objectPath?.split(".").at(-1) || "Select option";
    const processedGuiMode = qualifiers?.guiMode || "autocomplete";
    const processedOnlyEnumerators = qualifiers?.onlyEnumerators === true;

    return {
      enumerators: processedEnumerators,
      labels: processedLabels,
      guiLabel: processedGuiLabel,
      guiMode: processedGuiMode,
      onlyEnumerators: processedOnlyEnumerators,
    };
  }, [qualifiers, serverValue, objectPath]);

  const isRadioMode = guiMode === "multiLineRadio" || guiMode === "radio";

  const calculatedSx = useMemo(
    () => ({
      ...(sizeProp ? getFieldSizeStyles(sizeProp) : FULL_WIDTH_FIELD_STYLES),
      marginTop: 0,
      verticalAlign: "top",
      ...sx,
    }),
    [sizeProp, sx]
  );

  const handleCommit = useCallback(
    async (newValue: OptionValue) => {
      const result = await commit(newValue);
      if (result && !result.success) {
        setLocalValue(serverValue ?? "");
      }
    },
    [commit, serverValue, setLocalValue]
  );

  const handleAutocompleteChange = useCallback(
    async (
      _event: SyntheticEvent<Element, Event>,
      newValue: OptionValue | null,
      _reason: AutocompleteChangeReason
    ) => {
      if (newValue !== null && newValue !== localValue) {
        setLocalValue(newValue);
        await handleCommit(newValue);
      }
    },
    [localValue, handleCommit, setLocalValue]
  );

  const handleRadioChange = useCallback(
    async (event: React.ChangeEvent<HTMLInputElement>) => {
      const newValue = event.currentTarget.value;
      if (newValue && newValue !== localValue) {
        setLocalValue(newValue);
        await handleCommit(newValue);
      }
    },
    [localValue, handleCommit, setLocalValue]
  );

  const getOptionLabel = useCallback(
    (option: OptionValue): string => {
      if (option === "" || option === null || option === undefined) return "";
      const index = enumerators.indexOf(option);
      const label = labels[index];
      if (!label) {
        console.warn(`No label found for option: ${option}`);
        return String(option);
      }
      return label;
    },
    [enumerators, labels]
  );

  if (!isVisible || !enumerators.length || !labels.length) return null;

  if (isRadioMode) {
    return (
      <FieldWrapper ariaLabel={`${guiLabel} selection`}>
        <RadioGroup
          row={guiMode === "radio"}
          value={localValue}
          onChange={handleRadioChange}
          sx={calculatedSx}
          name={`radio-group-${itemName}`}
        >
          <FormControlLabel control={<></>} label={guiLabel} sx={{ mr: 2 }} />
          {enumerators.map((enumerator, index) => (
            <FormControlLabel
              key={`${enumerator}-${index}`}
              value={enumerator}
              control={
                <Radio
                  size="small"
                  disabled={isDisabled}
                  inputProps={{ "aria-label": getOptionLabel(enumerator) }}
                />
              }
              label={getOptionLabel(enumerator)}
            />
          ))}
        </RadioGroup>
        <ErrorTrigger item={item} job={job} />
      </FieldWrapper>
    );
  }

  return (
    <FieldWrapper ariaLabel={`${guiLabel} selection`}>
      <Autocomplete
        disabled={isDisabled}
        sx={calculatedSx}
        value={localValue}
        onChange={handleAutocompleteChange}
        onInputChange={(_event, newInputValue, reason) => {
          if (
            reason === "input" &&
            newInputValue !== localValue &&
            !onlyEnumerators
          ) {
            setLocalValue(newInputValue);
          }
        }}
        onBlur={async () => {
          if (!onlyEnumerators && localValue !== serverValue) {
            await handleCommit(localValue);
          }
        }}
        getOptionLabel={getOptionLabel}
        options={enumerators}
        size="small"
        freeSolo={!onlyEnumerators}
        disableClearable={onlyEnumerators}
        isOptionEqualToValue={(option, value) => option === value}
        renderInput={(params) => (
          <TextField
            {...params}
            error={validationColor === "error.light"}
            label={guiLabel}
            size="small"
            slotProps={{
              inputLabel: { shrink: true, disableAnimation: true },
            }}
            inputProps={{
              ...params.inputProps,
              "aria-label": guiLabel,
            }}
          />
        )}
      />
      <ErrorTrigger item={item} job={job} />
    </FieldWrapper>
  );
};
