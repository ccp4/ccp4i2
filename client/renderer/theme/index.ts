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
"use client";
import { createTheme, ThemeOptions } from "@mui/material/styles";
import { paletteOptions, customColors } from "./palette";
import { typographyOptions } from "./typography";
import { createComponentVariants } from "./variants";

// Base theme configuration (legacy - kept for backward compatibility)
const baseThemeOptions: ThemeOptions = {
  palette: paletteOptions,
  typography: typographyOptions,
};

// Create the base theme first
const baseTheme = createTheme(baseThemeOptions);

// Create the final theme with component variants
export const theme = createTheme({
  ...baseTheme,
  components: createComponentVariants(baseTheme),
});

// Export custom colors for use in components
export { customColors };

// Export theme provider and hooks
export { CCP4i2ThemeProvider, useTheme } from "./theme-provider";
export type { ThemeMode } from "./theme-provider";

export default theme;
