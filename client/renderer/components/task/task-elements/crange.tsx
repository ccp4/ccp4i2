import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { Box, Typography } from "@mui/material";
import { FIELD_SIZES } from "./field-sizes";
import { useInferredVisibility } from "./hooks/useInferredVisibility";

/**
 * CRangeElement — a compact "[label] [start] to [end]" inline row (e.g. a
 * resolution range).
 *
 * Renders as a plain inline row (no bordered container). The start/end CFloat
 * fields are full-width by default, so they MUST be wrapped in fixed-width
 * boxes here — otherwise two width:100% fields can't share a flex row and the
 * "end" field wraps onto its own line inconsistently depending on panel width.
 * Each unit is flexShrink:0 so a narrow panel wraps between whole units rather
 * than splitting a field from its label.
 */
export const CRangeElement: React.FC<CCP4i2TaskElementProps> = (props) => {
  const isVisible = useInferredVisibility(props.visibility);
  if (!isVisible) return null;

  return (
    <Box
      sx={{
        display: "flex",
        alignItems: "center",
        gap: 1,
        flexWrap: "wrap",
      }}
    >
      <Typography variant="body2" sx={{ flexShrink: 0 }}>
        {props.qualifiers?.guiLabel || "Resolution range"}
      </Typography>
      <Box sx={{ width: FIELD_SIZES.xs, flexShrink: 0 }}>
        <CCP4i2TaskElement
          job={props.job}
          itemName={`${props.itemName}.start`}
          qualifiers={{ guiLabel: " " }}
        />
      </Box>
      <Typography variant="body2" sx={{ flexShrink: 0 }}>
        to
      </Typography>
      <Box sx={{ width: FIELD_SIZES.xs, flexShrink: 0 }}>
        <CCP4i2TaskElement
          job={props.job}
          itemName={`${props.itemName}.end`}
          qualifiers={{ guiLabel: " " }}
        />
      </Box>
    </Box>
  );
};
