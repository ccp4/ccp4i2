import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { Box, Typography } from "@mui/material";
import { FIELD_SIZES } from "./field-sizes";
import { useInferredVisibility } from "./hooks/useInferredVisibility";

/**
 * CDmDomainElement — one rigid body (dm "domain") as a compact inline row:
 *
 *   segments [ cyclin:10-95,CDK:45-60 ]   mode [average ▾]
 *
 * A rigid body is a set of residue-range segments that move together, plus an
 * averaging mode. Segments are "role:first-last" ranges; a bare range uses the
 * single implicit role (homomer), and mixing roles expresses a cross-chain body
 * (the CDK C-helix travelling with the cyclin N-lobe). The role→chain mapping
 * for each NCS copy lives in the task's ASSEMBLY list, not here. Each sub-field
 * is wrapped in a fixed-width Box so a narrow panel wraps between whole units
 * rather than splitting a field from its label.
 */
export const CDmDomainElement: React.FC<CCP4i2TaskElementProps> = (props) => {
  const isVisible = useInferredVisibility(props.visibility);
  if (!isVisible) return null;

  const { job, itemName } = props;
  const sub = (name: string, label: string, width: string) => (
    <>
      <Typography variant="body2" sx={{ flexShrink: 0 }}>
        {label}
      </Typography>
      <Box sx={{ width, flexShrink: 0 }}>
        <CCP4i2TaskElement
          job={job}
          itemName={`${itemName}.${name}`}
          qualifiers={{ guiLabel: " " }}
        />
      </Box>
    </>
  );

  return (
    <Box sx={{ display: "flex", alignItems: "center", gap: 1, flexWrap: "wrap" }}>
      {sub("segments", "segments", FIELD_SIZES.lg)}
      {sub("mode", "mode", FIELD_SIZES.sm)}
    </Box>
  );
};
