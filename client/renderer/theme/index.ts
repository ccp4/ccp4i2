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
