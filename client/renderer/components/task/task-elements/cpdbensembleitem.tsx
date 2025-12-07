import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { useJob } from "../../../utils";
import { Box, Stack, Typography } from "@mui/material";
import { CSimpleDataFileElement } from "./csimpledatafile";
import { useInferredVisibility } from "./hooks/useInferredVisibility";

export const CPdbEnsembleItemElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => {
  const { job, itemName } = props;
  const { useTaskItem, getValidationColor } = useJob(job.id);
  const { item } = useTaskItem(itemName);

  const isVisible = useInferredVisibility(props.visibility);

  // Get the validation color for the entire CPdbEnsembleItem (not just the structure child)
  const itemValidationColor = getValidationColor(item);
  const hasValidationError = itemValidationColor === "error.light";

  return isVisible && item ? (
    <Box
      sx={{
        border: "3px solid",
        borderColor: itemValidationColor,
        borderRadius: "0.5rem",
        p: 0.5,
      }}
    >
      <CSimpleDataFileElement
        {...props}
        itemName={`${item._objectPath}.structure`}
        qualifiers={{
          guiLabel: "Coordinates",
        }}
        forceExpanded={hasValidationError}
      >
        <Stack spacing={1}>
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
          {hasValidationError && (
            <Typography variant="caption" color="error">
              At least one of Identity or RMS must be set to a non-zero value
            </Typography>
          )}
        </Stack>
      </CSimpleDataFileElement>
    </Box>
  ) : null;
};
