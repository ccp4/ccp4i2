/**
 * #613, the task Validation tab:
 * - its summary bar was a fixed light grey (#fafafa) in dark mode too;
 * - warnings (severity 1 from the validation parser) were filed under a
 *   collapsed "Info" group, and the "Warnings" group keyed on a severity 3
 *   that nothing produces.
 */
import React from "react";
import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import { ThemeProvider, createTheme } from "@mui/material/styles";
import {
  darkCustomColors,
  darkPaletteOptions,
  lightPaletteOptions,
} from "../theme/palette";

let themeMode: "light" | "dark" = "dark";

vi.mock("../utils", () => ({
  useJob: () => ({
    validation: {
      "refmac.container.inputData.XYZIN": {
        messages: ["XYZIN: Data has undefined value"],
        maxSeverity: 2,
      },
      "refmac.container.inputData.FREERFLAG": {
        messages: ["FREERFLAG: Free R flag is strongly recommended"],
        maxSeverity: 1,
      },
    },
  }),
}));
vi.mock("../api", () => ({
  useApi: () => ({
    get_pretty_endpoint_xml: () => ({
      data: "<xml/>",
      error: undefined,
      mutate: vi.fn(),
    }),
  }),
}));
vi.mock("../app-context", () => ({
  useCCP4i2Window: () => ({ devMode: false }),
}));
vi.mock("../theme/theme-provider", async () => {
  const palette = await import("../theme/palette");
  return {
    useTheme: () => ({
      mode: themeMode,
      customColors:
        themeMode === "dark" ? palette.darkCustomColors : palette.lightCustomColors,
    }),
  };
});
vi.mock("@monaco-editor/react", () => ({
  Editor: ({ value }: { value: string }) => <pre>{value}</pre>,
}));

import { ValidationViewer } from "../components/validation-viewer";

const renderIn = (mode: "light" | "dark") => {
  themeMode = mode;
  const theme = createTheme({
    palette: mode === "dark" ? darkPaletteOptions : lightPaletteOptions,
  });
  return render(
    <ThemeProvider theme={theme}>
      <ValidationViewer job={{ id: 1 } as any} />
    </ThemeProvider>
  );
};

describe("ValidationViewer", () => {
  it("dark mode: the issues bar takes the dark surface, not light grey", () => {
    renderIn("dark");
    const heading = screen.getByText(/2 validation issues/);
    const bar = heading.closest(".MuiStack-root")!.parentElement!;
    const bg = getComputedStyle(bar).backgroundColor;
    expect(bg).not.toBe("rgb(250, 250, 250)");
    expect(darkCustomColors.ui.veryLightGray.toLowerCase()).toBe("#333333");
    expect(bg).toBe("rgb(51, 51, 51)");
  });

  it("files a warning under Warnings, expanded, beside Errors", () => {
    renderIn("light");
    expect(screen.getByText("Errors")).toBeTruthy();
    expect(screen.getByText("Warnings")).toBeTruthy();
    expect(screen.queryByText("Info")).toBeNull();
    // Both groups open by default: their messages are on the page
    expect(screen.getByText(/Free R flag is strongly recommended/)).toBeTruthy();
    expect(screen.getByText(/Not defined/)).toBeTruthy();
  });
});
