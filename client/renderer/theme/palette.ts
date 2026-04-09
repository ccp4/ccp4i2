import { PaletteOptions } from "@mui/material/styles";

// Light mode palette
export const lightPaletteOptions: PaletteOptions = {
  mode: "light",
  primary: {
    main: "rgb(92, 149, 165)",
    dark: "rgb(39, 62, 69)",
    contrastText: "rgb(0, 0, 0)",
  },
  secondary: {
    main: "rgb(156, 39, 176)",
    light: "rgb(225, 190, 231)",
    dark: "rgb(123, 31, 162)",
    contrastText: "rgb(255, 255, 255)",
  },
  success: {
    main: "rgb(70, 164, 75)",
    light: "rgb(193, 226, 214)",
    dark: "rgb(27, 94, 32)",
    contrastText: "rgb(255, 255, 255)",
  },
  background: {
    default: "#ffffff",
    paper: "#ffffff",
  },
  text: {
    primary: "rgba(0, 0, 0, 0.87)",
    secondary: "rgba(0, 0, 0, 0.6)",
  },
};

// Dark mode palette
export const darkPaletteOptions: PaletteOptions = {
  mode: "dark",
  primary: {
    main: "rgb(129, 199, 212)",
    dark: "rgb(79, 144, 154)",
    light: "rgb(179, 229, 252)",
    contrastText: "rgba(0, 0, 0, 0.87)",
  },
  secondary: {
    main: "rgb(206, 147, 216)",
    light: "rgb(245, 208, 250)",
    dark: "rgb(171, 71, 188)",
    contrastText: "rgba(0, 0, 0, 0.87)",
  },
  success: {
    main: "rgb(129, 199, 132)",
    light: "rgb(200, 230, 201)",
    dark: "rgb(51, 105, 30)",
    contrastText: "rgba(0, 0, 0, 0.87)",
  },
  background: {
    default: "#121212",
    paper: "#1e1e1e",
  },
  text: {
    primary: "#ffffff",
    secondary: "rgba(255, 255, 255, 0.7)",
  },
};

// Light mode custom colors
export const lightCustomColors = {
  // Amino acid colors for sequence alignment viewer (Clustal color scheme)
  aminoAcids: {
    A: "#80a0f0", // Alanine
    R: "#f01505", // Arginine
    N: "#00ff00", // Asparagine
    D: "#c048c0", // Aspartic acid
    C: "#f08080", // Cysteine
    Q: "#00ff00", // Glutamine
    E: "#c048c0", // Glutamic acid
    G: "#f09048", // Glycine
    H: "#15a4a4", // Histidine
    I: "#80a0f0", // Isoleucine
    L: "#80a0f0", // Leucine
    K: "#f01505", // Lysine
    M: "#80a0f0", // Methionine
    F: "#80a0f0", // Phenylalanine
    P: "#ffff00", // Proline
    S: "#00ff00", // Serine
    T: "#00ff00", // Threonine
    W: "#80a0f0", // Tryptophan
    Y: "#15a4a4", // Tyrosine
    V: "#80a0f0", // Valine
    "-": "#ffffff", // Gap
    ".": "#ffffff", // Gap
    " ": "#ffffff", // Space
  },
  // UI background colors
  ui: {
    lightGray: "#f9f9f9",
    veryLightGray: "#f5f5f5",
    mediumGray: "#e0e0e0",
    lightBlue: "#1976d2",
    white: "#ffffff",
  },
};

// Dark mode custom colors
export const darkCustomColors = {
  // Amino acid colors for sequence alignment viewer (adjusted for dark backgrounds)
  aminoAcids: {
    A: "#6080d0", // Alanine - slightly darker
    R: "#d01505", // Arginine - slightly darker
    N: "#00e000", // Asparagine - slightly darker
    D: "#a040a0", // Aspartic acid - slightly darker
    C: "#d08080", // Cysteine - slightly darker
    Q: "#00e000", // Glutamine - slightly darker
    E: "#a040a0", // Glutamic acid - slightly darker
    G: "#d09048", // Glycine - slightly darker
    H: "#158484", // Histidine - slightly darker
    I: "#6080d0", // Isoleucine - slightly darker
    L: "#6080d0", // Leucine - slightly darker
    K: "#d01505", // Lysine - slightly darker
    M: "#6080d0", // Methionine - slightly darker
    F: "#6080d0", // Phenylalanine - slightly darker
    P: "#e0e000", // Proline - slightly darker
    S: "#00e000", // Serine - slightly darker
    T: "#00e000", // Threonine - slightly darker
    W: "#6080d0", // Tryptophan - slightly darker
    Y: "#158484", // Tyrosine - slightly darker
    V: "#6080d0", // Valine - slightly darker
    "-": "#404040", // Gap - dark gray for dark mode
    ".": "#404040", // Gap - dark gray for dark mode
    " ": "#404040", // Space - dark gray for dark mode
  },
  // UI background colors for dark mode
  ui: {
    lightGray: "#2a2a2a",
    veryLightGray: "#333333",
    mediumGray: "#555555",
    lightBlue: "#42a5f5",
    white: "#ffffff",
  },
};

// Legacy exports for backward compatibility
export const paletteOptions = lightPaletteOptions;
export const customColors = lightCustomColors;
