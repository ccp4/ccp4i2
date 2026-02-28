import React, { PropsWithChildren } from "react";
import { Box, Typography, SxProps, Theme } from "@mui/material";

interface InlineFieldProps {
  /** Label text shown before the widget */
  label?: string;
  /** Italic hint shown after the widget — string or ReactNode */
  hint?: React.ReactNode;
  /** Width of the widget container (default "8rem") */
  width?: string;
  /** Additional content rendered after the hint (e.g. a second widget) */
  after?: React.ReactNode;
  /** Override sx on the outer flex container */
  sx?: SxProps<Theme>;
}

/**
 * A horizontal row with: [label] [fixed-width child] [hint] [after]
 *
 * Replaces the common pattern of:
 *   <Box sx={{ display: "flex", alignItems: "center", gap: 1, flexWrap: "wrap" }}>
 *     <Typography variant="body1">Label</Typography>
 *     <Box sx={{ width: "8rem" }}>
 *       <CCP4i2TaskElement ... />
 *     </Box>
 *     <Typography variant="body2" sx={{ fontStyle: "italic" }}>[hint]</Typography>
 *   </Box>
 */
export const InlineField: React.FC<PropsWithChildren<InlineFieldProps>> = ({
  label,
  hint,
  width = "8rem",
  after,
  sx,
  children,
}) => (
  <Box
    sx={{
      display: "flex",
      alignItems: "center",
      gap: 1,
      flexWrap: "wrap",
      ...sx,
    }}
  >
    {label && <Typography variant="body1">{label}</Typography>}
    <Box sx={{ width }}>{children}</Box>
    {hint &&
      (typeof hint === "string" ? (
        <Typography variant="body2" sx={{ fontStyle: "italic" }}>
          {hint}
        </Typography>
      ) : (
        hint
      ))}
    {after}
  </Box>
);

export default InlineField;
