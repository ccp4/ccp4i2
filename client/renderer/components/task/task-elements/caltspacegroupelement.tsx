import { useEffect, useState } from "react";
import {
  Autocomplete,
  Card,
  CardContent,
  CardHeader,
  TextField,
} from "@mui/material";

import { CCP4i2TaskElementProps } from "./task-element";
import { SpaceGroup, spaceGroups } from "../../../spacegroups";
import { useTaskInterface } from "../../../providers/task-provider";
import { useContainerField } from "./hooks/useContainerField";

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

  const [value, setValue] = useState<SpaceGroup | undefined>(spaceGroups[0]);

  useEffect(() => {
    if (typeof serverValue === "string" && serverValue) {
      const normalized = serverValue.replace(/\s+/g, "");
      setValue(
        spaceGroups.find(
          (sg: SpaceGroup) =>
            sg.name === serverValue ||
            sg.name.replace(/\s+/g, "") === normalized
        )
      );
    }
  }, [serverValue]);

  const handleInputChanged = async (arg: SpaceGroup) => {
    setValue(arg);
    await commit(arg.name);
  };

  if (!isVisible) return null;

  return (
    <Card sx={{ border: "3px solid", borderColor: validationColor }}>
      <CardHeader title={qualifiers?.guiLabel} sx={{ borderColor: validationColor }} />
      <CardContent sx={{ my: 0, py: 0 }}>
        {item && (
          <Autocomplete
            sx={{
              mt: 1,
              backgroundColor: inFlight ? "#ffeebe" : "palette.common.white",
            }}
            id="autocomplete-spacegroup"
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
              <TextField {...params} label="Space groups" />
            )}
          />
        )}
      </CardContent>
    </Card>
  );
};
