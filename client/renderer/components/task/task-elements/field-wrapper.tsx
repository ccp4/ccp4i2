/*
 * Copyright (C) 2025 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
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
