import { useEffect, useState } from "react";
import { Autocomplete, TextField } from "@mui/material";

import { CCP4i2TaskElementProps } from "./task-element";
import { SpaceGroup, spaceGroups } from "../../../spacegroups";
import { useTaskInterface } from "../../../providers/task-provider";
import { useContainerField } from "./hooks/useContainerField";
import { FieldShell } from "./field-shell";

export const CAltSpaceGroupElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => {
  const { job, itemName, qualifiers, visibility, disabled, onChange } = props;

  const {
    item,
    serverValue,
    isVisible,
    isDisabled,
    validationColor,
    commit,
  } = useContainerField<string>({
    job,
    itemName,
    visibility,
    disabled,
    onChange,
  });
  const { inFlight } = useTaskInterface();

  // Empty until the job has a space group. (It started on the list's first
  // entry, so an unset space group read "P 1", as if P 1 had been chosen.)
  const [value, setValue] = useState<SpaceGroup | null>(null);

  useEffect(() => {
    if (typeof serverValue === "string" && serverValue) {
      const normalized = serverValue.replace(/\s+/g, "");
      setValue(
        spaceGroups.find(
          (sg: SpaceGroup) =>
            sg.name === serverValue ||
            sg.name.replace(/\s+/g, "") === normalized
        ) ?? null
      );
    } else {
      setValue(null);
    }
  }, [serverValue]);

  const handleInputChanged = async (arg: SpaceGroup) => {
    setValue(arg);
    await commit(arg.name);
  };

  if (!isVisible) return null;

  return (
    <FieldShell title={qualifiers?.guiLabel} borderColor={validationColor}>
      {item && (
        <Autocomplete
          sx={{
            backgroundColor: inFlight ? "warning.light" : "background.paper",
          }}
          id="autocomplete-spacegroup"
          size="small"
          disabled={isDisabled}
          multiple={false}
          options={spaceGroups}
          getOptionLabel={(option: SpaceGroup) => option.name}
          getOptionKey={(option: SpaceGroup) => option.name}
          value={value}
          style={{ minWidth: "15rem" }}
          onChange={(
            _event: React.SyntheticEvent<Element, Event>,
            newValue: SpaceGroup | null
          ) => {
            if (newValue) handleInputChanged(newValue);
          }}
          renderInput={(params: any) => (
            <TextField {...params} label="Space groups" size="small" />
          )}
        />
      )}
    </FieldShell>
  );
};
