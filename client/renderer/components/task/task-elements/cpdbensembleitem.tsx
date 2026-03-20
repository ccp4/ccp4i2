import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { useJob } from "../../../utils";
import { Box, Typography } from "@mui/material";
import { CSimpleDataFileElement } from "./csimpledatafile";
import { useInferredVisibility } from "./hooks/useInferredVisibility";
import { FieldRow } from "./field-row";

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
        borderLeft: hasValidationError ? "3px solid" : "2px solid",
        borderColor: itemValidationColor,
        borderRadius: 0.5,
        pl: 1,
        py: 0.5,
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
        <Box sx={{ display: "flex", flexDirection: "column", gap: 0.5 }}>
          <FieldRow equalWidth={false} size="sm">
            <CCP4i2TaskElement
              {...props}
              sx={{ my: 0, py: 0 }}
              itemName={`${item._objectPath}.identity_to_target`}
              qualifiers={{ guiLabel: "Identity" }}
            />
            <CCP4i2TaskElement
              {...props}
              sx={{ my: 0, py: 0 }}
              itemName={`${item._objectPath}.rms_to_target`}
              qualifiers={{
                guiLabel: "Rms",
              }}
            />
          </FieldRow>
          {hasValidationError && (
            <Typography variant="caption" color="error">
              At least one of Identity or RMS must be set to a non-zero value
            </Typography>
          )}
        </Box>
      </CSimpleDataFileElement>
    </Box>
  ) : null;
};
