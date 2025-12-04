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
import { useParameterChangeIntent } from "../../../providers/parameter-change-intent-provider";
import { inferFieldSize, getFieldSizeStyles } from "./field-sizes";
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
  const { setIntentForPath, clearIntentForPath, wasRecentlyChanged } =
    useParameterChangeIntent();

  // Local state
  const [localValue, setLocalValue] = useState<OptionValue>(item?._value || "");
  const [isSubmitting, setIsSubmitting] = useState(false);

  // Get object path for intent tracking
  const objectPath = useMemo(() => item?._objectPath || null, [item]);

  // Sync local state with item value changes - but skip if recently edited
  useEffect(() => {
    if (item?._value !== undefined) {
      // Don't overwrite local state if user just edited this field
      // Use a longer timeout (10 seconds) to account for:
      // - API call latency
      // - Container mutation and SWR revalidation
      // - Multiple cascading updates from parent components
      if (objectPath && wasRecentlyChanged(objectPath, 10000)) {
        return;
      }
      setLocalValue(item._value);
    }
  }, [item?._value, objectPath, wasRecentlyChanged]);

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

  // Calculate field size - use explicit prop or infer from item/qualifiers
  const fieldSize = useMemo(
    () => sizeProp || inferFieldSize(item, qualifiers),
    [sizeProp, item, qualifiers]
  );

  const calculatedSx = useMemo(
    () => ({
      ...getFieldSizeStyles(fieldSize),
      // Ensure Autocomplete aligns with TextField in grid/flex layouts
      // MUI Autocomplete adds extra wrapper that can cause vertical offset
      marginTop: 0,
      verticalAlign: "top",
      ...sx,
    }),
    [fieldSize, sx]
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

      // Record intent BEFORE making the API call
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
        const result: SetParameterResponse | undefined = await updateFn(
          setParameterArg
        );

        if (result && !result.success) {
          setMessage(`Unacceptable value provided: "${newValue}"`);
          setLocalValue(item._value); // Revert to original value
          clearIntentForPath(objectPath);
        } else if (result?.success && result.data?.updated_item && onChange) {
          // Don't clear intent here - let the timeout handle cleanup
          // This prevents bounce-back when container mutation triggers re-render
          await onChange(result.data.updated_item);
        }
        // On success without onChange, don't clear intent - let timeout handle it
      } catch (error) {
        const errorMessage =
          error instanceof Error ? error.message : String(error);
        setMessage(`Error updating parameter: ${errorMessage}`);
        console.error("Parameter update failed:", error);
        setLocalValue(item._value); // Revert to original value
        clearIntentForPath(objectPath);
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
      job.id,
      setIntentForPath,
      clearIntentForPath,
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
            // Record intent when typing to prevent container refetch overwrites
            if (objectPath) {
              setIntentForPath({
                jobId: job.id,
                parameterPath: objectPath,
                reason: "UserEdit",
                previousValue: item?._value,
              });
            }
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
