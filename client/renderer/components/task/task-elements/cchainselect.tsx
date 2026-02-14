import React, { useCallback, useEffect, useMemo, useState } from "react";
import { Autocomplete, Chip, TextField, SxProps, Theme } from "@mui/material";
import { Job } from "../../../types/models";
import { useJob, SetParameterResponse } from "../../../utils";
import { useTaskInterface } from "../../../providers/task-provider";
import { usePopcorn } from "../../../providers/popcorn-provider";
import { FULL_WIDTH_FIELD_STYLES } from "./field-sizes";
import { FieldWrapper } from "./field-wrapper";
import { ErrorTrigger } from "./error-info";

interface CChainSelectProps {
  job: Job;
  /** Dotted path to the CString parameter, e.g. "prosmartProtein.CHAINLIST_1" */
  itemName: string;
  /** Available chain IDs to offer as options (from digest composition) */
  options: string[];
  /** Label for the field */
  label?: string;
  sx?: SxProps<Theme>;
  visibility?: boolean | (() => boolean);
  disabled?: boolean | (() => boolean);
}

/** Parse a comma-separated chain string into an array of trimmed chain IDs */
const parseChains = (value: any): string[] => {
  if (!value || typeof value !== "string") return [];
  return value
    .split(",")
    .map((s) => s.trim())
    .filter(Boolean);
};

/**
 * Multi-select chain selector that stores its value as a comma-separated CString.
 *
 * Renders an MUI Autocomplete with `multiple` and chip tags. Options are
 * populated from the parent (typically from the XYZIN digest composition).
 */
export const CChainSelectElement: React.FC<CChainSelectProps> = ({
  job,
  itemName,
  options,
  label = "Chains",
  sx,
  visibility,
  disabled: disabledProp,
}) => {
  const { useTaskItem, setParameter, getValidationColor } = useJob(job.id);
  const { item } = useTaskItem(itemName);
  const { setMessage } = usePopcorn();
  const { inFlight, setInFlight } = useTaskInterface();

  const [localValue, setLocalValue] = useState<string[]>(() =>
    parseChains(item?._value)
  );
  const [isSubmitting, setIsSubmitting] = useState(false);

  const objectPath = useMemo(() => item?._objectPath || null, [item]);

  // Sync local state when the server value changes
  useEffect(() => {
    if (item?._value !== undefined) {
      setLocalValue(parseChains(item._value));
    }
  }, [item?._value]);

  const isVisible = useMemo(() => {
    if (typeof visibility === "function") return visibility();
    return visibility !== false;
  }, [visibility]);

  const isDisabled = useMemo(() => {
    const base =
      typeof disabledProp === "function" ? disabledProp() : disabledProp;
    return base || inFlight || isSubmitting || job.status !== 1;
  }, [disabledProp, inFlight, isSubmitting, job.status]);

  const calculatedSx = useMemo(
    () => ({ ...FULL_WIDTH_FIELD_STYLES, marginTop: 0, ...sx }),
    [sx]
  );

  const handleChange = useCallback(
    async (_event: React.SyntheticEvent, newValue: string[]) => {
      setLocalValue(newValue);

      if (!objectPath) return;

      const csvValue = newValue.join(",");

      setInFlight(true);
      setIsSubmitting(true);
      try {
        const result: SetParameterResponse | undefined = await setParameter({
          object_path: objectPath,
          value: csvValue,
        });
        if (result && !result.success) {
          setMessage(`Unacceptable value: "${csvValue}"`);
          setLocalValue(parseChains(item._value));
        }
      } catch (error) {
        const msg = error instanceof Error ? error.message : String(error);
        setMessage(`Error updating chains: ${msg}`);
        setLocalValue(parseChains(item._value));
      } finally {
        setInFlight(false);
        setIsSubmitting(false);
      }
    },
    [objectPath, item, setParameter, setInFlight, setMessage]
  );

  if (!isVisible) return null;

  return (
    <FieldWrapper ariaLabel={`${label} chain selection`}>
      <Autocomplete
        multiple
        disabled={isDisabled}
        sx={calculatedSx}
        value={localValue}
        onChange={handleChange}
        options={options}
        size="small"
        disableCloseOnSelect
        renderTags={(value, getTagProps) =>
          value.map((chain, index) => {
            const { key, ...tagProps } = getTagProps({ index });
            return (
              <Chip key={key} label={chain} size="small" {...tagProps} />
            );
          })
        }
        renderInput={(params) => (
          <TextField
            {...params}
            label={label}
            size="small"
            slotProps={{
              inputLabel: { shrink: true, disableAnimation: true },
            }}
          />
        )}
      />
      <ErrorTrigger item={item} job={job} />
    </FieldWrapper>
  );
};

export default CChainSelectElement;
