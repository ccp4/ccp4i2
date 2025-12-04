import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { useJob } from "../../../utils";
import {
  Card,
  CardContent,
  CardHeader,
  Grid2,
  Paper,
  Stack,
} from "@mui/material";
import { ErrorInfo } from "./error-info";
import { useMemo } from "react";
import { CSimpleDataFileElement } from "./csimpledatafile";

export const CPdbEnsembleItemElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => {
  const { job, itemName } = props;
  const { useTaskItem, getValidationColor } = useJob(job.id);
  const { item } = useTaskItem(itemName);

  const inferredVisibility = useMemo(() => {
    if (!props.visibility) return true;
    if (typeof props.visibility === "function") {
      return props.visibility();
    }
    return props.visibility;
  }, [props.visibility]);

  return inferredVisibility && item ? (
    <CSimpleDataFileElement
      {...props}
      sx={{
        border: "3px solid",
        borderColor: getValidationColor(item),
      }}
      itemName={`${item._objectPath}.structure`}
      qualifiers={{
        guiLabel: "Coordinates",
      }}
    >
      <Stack direction="row" spacing={2}>
        <CCP4i2TaskElement
          {...props}
          sx={{ my: 0, py: 0, minWidth: "10rem" }}
          itemName={`${item._objectPath}.identity_to_target`}
          qualifiers={{ guiLabel: "Identity" }}
        />
        <CCP4i2TaskElement
          {...props}
          sx={{ my: 0, py: 0, minWidth: "10rem" }}
          itemName={`${item._objectPath}.rms_to_target`}
          qualifiers={{
            guiLabel: "Rms",
          }}
        />
      </Stack>
    </CSimpleDataFileElement>
  ) : null;
};
