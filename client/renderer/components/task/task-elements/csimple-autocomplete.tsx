import React, {
  SyntheticEvent,
  useCallback,
  useMemo,
  useState,
  useEffect,
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
import { useJob, SetParameterResponse } from "../../../utils";
import { ErrorTrigger } from "./error-info";
import { useTaskInterface } from "../../../providers/task-provider";
import { usePopcorn } from "../../../providers/popcorn-provider";
import {
  FULL_WIDTH_FIELD_STYLES,
  getFieldSizeStyles,
  FieldSize,
} from "./field-sizes";
import { FieldWrapper } from "./field-wrapper";

type OptionValue = string | number;
type GuiMode = "multiLineRadio" | "radio" | "autocomplete";

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
    useTaskItem,
    getValidationColor,
    setParameter,
    setParameterNoMutate,
  } = useJob(job.id);
  const { item } = useTaskItem(itemName);
  const { setMessage } = usePopcorn();
  const { inFlight, setInFlight } = useTaskInterface();

  // Local state
  const [localValue, setLocalValue] = useState<OptionValue>(item?._value || "");
  const [isSubmitting, setIsSubmitting] = useState(false);

  // Get object path
  const objectPath = useMemo(() => item?._objectPath || null, [item]);

  // Sync local state with item value changes directly (no intent guards needed with local patching)
  useEffect(() => {
    if (item?._value !== undefined) {
      setLocalValue(item._value);
    }
  }, [item?._value]);

  // Process enumerators and labels
  const { enumerators, labels, guiLabel, guiMode, onlyEnumerators } =
    useMemo(() => {
      const rawEnumerators =
        qualifiers?.enumerators?.map((element: any) => {
          if (typeof element === "string" || element instanceof String) {
            return element.trim();
          }
          return element;
        }) || [];

      // Include current value if not in enumerators
      const processedEnumerators = [...rawEnumerators];
      if (item?._value && !processedEnumerators.includes(item._value)) {
        processedEnumerators.push(item._value);
      }

      // Process labels
      let processedLabels = processedEnumerators.map((item) => `${item}`);
      if (
        qualifiers?.menuText &&
        qualifiers.menuText.length === processedEnumerators.length
      ) {
        processedLabels = qualifiers.menuText.map((text: string) =>
          text.trim()
        );
      }

      // Process GUI label
      const processedGuiLabel =
        qualifiers?.guiLabel ||
        item?._objectPath?.split(".").at(-1) ||
        "Select option";

      // Determine GUI mode
      const processedGuiMode = qualifiers?.guiMode || "autocomplete";

      // Check onlyEnumerators setting
      const processedOnlyEnumerators = qualifiers?.onlyEnumerators === true;

      return {
        enumerators: processedEnumerators,
        labels: processedLabels,
        guiLabel: processedGuiLabel,
        guiMode: processedGuiMode,
        onlyEnumerators: processedOnlyEnumerators,
      };
    }, [item, qualifiers]);

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

  const validationColor = useMemo(
    () => getValidationColor(item),
    [getValidationColor, item]
  );

  const isRadioMode = useMemo(
    () => guiMode === "multiLineRadio" || guiMode === "radio",
    [guiMode]
  );

  // Full-width by default - containers control the width
  // Only use explicit sizeProp for legacy compatibility
  const calculatedSx = useMemo(
    () => ({
      ...(sizeProp ? getFieldSizeStyles(sizeProp) : FULL_WIDTH_FIELD_STYLES),
      // Ensure Autocomplete aligns with TextField in grid/flex layouts
      // MUI Autocomplete adds extra wrapper that can cause vertical offset
      marginTop: 0,
      verticalAlign: "top",
      ...sx,
    }),
    [sizeProp, sx]
  );

  // Parameter update handler
  const handleParameterUpdate = useCallback(
    async (newValue: OptionValue) => {
      if (!objectPath) {
        console.error("No object path available for parameter update");
        return;
      }

      const setParameterArg = {
        object_path: objectPath,
        value: newValue,
      };

      setInFlight(true);
      setIsSubmitting(true);

      try {
        const updateFn = suppressMutations
          ? setParameterNoMutate
          : setParameter;
        const result: SetParameterResponse | undefined = await updateFn(
          setParameterArg
        );

        if (result && !result.success) {
          setMessage(`Unacceptable value provided: "${newValue}"`);
          setLocalValue(item._value); // Revert to original value
        } else if (result?.success && result.data?.updated_item && onChange) {
          await onChange(result.data.updated_item);
        }
      } catch (error) {
        const errorMessage =
          error instanceof Error ? error.message : String(error);
        setMessage(`Error updating parameter: ${errorMessage}`);
        console.error("Parameter update failed:", error);
        setLocalValue(item._value); // Revert to original value
      } finally {
        setInFlight(false);
        setIsSubmitting(false);
      }
    },
    [
      objectPath,
      item,
      suppressMutations,
      setParameterNoMutate,
      setParameter,
      setInFlight,
      setIsSubmitting,
      setMessage,
      onChange,
    ]
  );

  // Event handlers
  const handleAutocompleteChange = useCallback(
    async (
      event: SyntheticEvent<Element, Event>,
      newValue: OptionValue | null,
      reason: AutocompleteChangeReason
    ) => {
      if (newValue !== null && newValue !== localValue) {
        setLocalValue(newValue);
        await handleParameterUpdate(newValue);
      }
    },
    [localValue, handleParameterUpdate]
  );

  const handleRadioChange = useCallback(
    async (event: React.ChangeEvent<HTMLInputElement>) => {
      const newValue = event.currentTarget.value;
      if (newValue && newValue !== localValue) {
        setLocalValue(newValue);
        await handleParameterUpdate(newValue);
      }
    },
    [localValue, handleParameterUpdate]
  );

  const getOptionLabel = useCallback(
    (option: OptionValue): string => {
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

  // Early return if not visible or no options
  if (!isVisible || !enumerators.length || !labels.length) {
    return null;
  }

  // Render radio group
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

  // Render autocomplete
  return (
    <FieldWrapper ariaLabel={`${guiLabel} selection`}>
      <Autocomplete
        disabled={isDisabled}
        sx={calculatedSx}
        value={localValue}
        onChange={handleAutocompleteChange}
        onInputChange={(_event, newInputValue, reason) => {
          // Handle free text input only if not restricted to enumerators
          if (
            reason === "input" &&
            newInputValue !== localValue &&
            !onlyEnumerators
          ) {
            setLocalValue(newInputValue);
          }
        }}
        onBlur={async () => {
          // Update parameter when user finishes typing (only for freeSolo mode)
          if (!onlyEnumerators && localValue !== item?._value) {
            await handleParameterUpdate(localValue);
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
              inputLabel: {
                shrink: true,
                disableAnimation: true,
              },
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
