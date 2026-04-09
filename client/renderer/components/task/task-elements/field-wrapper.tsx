import React, { PropsWithChildren } from "react";
import { Box, SxProps, Theme } from "@mui/material";

interface FieldWrapperProps {
  /** Additional sx props */
  sx?: SxProps<Theme>;
  /** ARIA label for accessibility */
  ariaLabel?: string;
}

/**
 * Minimal wrapper for form field components.
 *
 * Spacing is handled by the parent container's gap properties,
 * so this wrapper just provides alignment and accessibility.
 */
export const FieldWrapper: React.FC<PropsWithChildren<FieldWrapperProps>> = ({
  children,
  sx,
  ariaLabel,
}) => (
  <Box
    sx={{
      display: "flex",
      alignItems: "flex-start",
      ...sx,
    }}
    role="group"
    aria-label={ariaLabel}
  >
    {children}
  </Box>
);

export default FieldWrapper;
