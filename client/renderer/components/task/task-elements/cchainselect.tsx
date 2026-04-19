import React, { useCallback, useMemo } from "react";
import { Autocomplete, Chip, TextField, SxProps, Theme } from "@mui/material";

import { Job } from "../../../types/models";
import { FULL_WIDTH_FIELD_STYLES } from "./field-sizes";
import { FieldWrapper } from "./field-wrapper";
import { ErrorTrigger } from "./error-info";
import {
  useContainerField,
  useSyncedLocalValue,
} from "./hooks/useContainerField";

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

const parseChains = (value: any): string[] => {
  if (!value || typeof value !== "string") return [];
  return value.split(",").map((s) => s.trim()).filter(Boolean);
};

/**
 * Multi-select chain selector that stores its value as a comma-separated CString.
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
  const { item, serverValue, isVisible, isDisabled, commit } =
    useContainerField<string>({
      job,
      itemName,
      visibility,
      disabled: disabledProp,
    });

  const serverChains = useMemo(() => parseChains(serverValue), [serverValue]);
  const [localValue, setLocalValue] =
    useSyncedLocalValue<string[]>(serverChains);

  const calculatedSx = useMemo(
    () => ({ ...FULL_WIDTH_FIELD_STYLES, marginTop: 0, ...sx }),
    [sx]
  );

  const handleChange = useCallback(
    async (_event: React.SyntheticEvent, newValue: string[]) => {
      setLocalValue(newValue);
      const csvValue = newValue.join(",");
      const result = await commit(csvValue);
      if (result && !result.success) {
        setLocalValue(parseChains(serverValue));
      }
    },
    [commit, serverValue, setLocalValue]
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
            return <Chip key={key} label={chain} size="small" {...tagProps} />;
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
