import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { Box, Typography } from "@mui/material";
import { FIELD_SIZES } from "./field-sizes";
import { useInferredVisibility } from "./hooks/useInferredVisibility";

/**
 * CDmDomainElement — one multi-domain-NCS domain as a compact inline row:
 *
 *   Chain [A]  residues [340] to [485]  mode [average ▾]
 *
 * Composed from the inherited CResidueRange fields (chainId / firstRes /
 * lastRes) plus the averaging `mode` enumerator. Each sub-field is wrapped in a
 * fixed-width Box (its CString/COneWord renderer is full-width by default), so a
 * narrow panel wraps between whole units rather than splitting a field from its
 * label -- same idiom as CRangeElement.
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
      {sub("chainId", "Chain", FIELD_SIZES.xs)}
      {sub("firstRes", "residues", FIELD_SIZES.xs)}
      <Typography variant="body2" sx={{ flexShrink: 0 }}>
        to
      </Typography>
      {sub("lastRes", "", FIELD_SIZES.xs)}
      {sub("mode", "mode", FIELD_SIZES.sm)}
    </Box>
  );
};
