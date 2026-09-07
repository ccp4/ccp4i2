/**
 * Alpha tracker: "Error styling in dark mode: the error details shows white
 * text on a light background." The Details and stack-trace <pre> blocks were
 * painted `grey.100` — a fixed light grey — under `text.primary`, which is
 * white in the dark palette. Every surface in the panel now resolves against
 * the active theme; these tests pin that down in both modes.
 *
 * jsdom resolves emotion's injected <style> rules through getComputedStyle,
 * so the assertions are on the colours actually applied to the elements, not
 * on the sx objects.
 */
import React from "react";
import { describe, it, expect } from "vitest";
import { render } from "@testing-library/react";
import { ThemeProvider, createTheme } from "@mui/material/styles";
import Diagnostic, { diagnosticSurfaceColors } from "../components/diagnostic";
import { darkPaletteOptions, lightPaletteOptions } from "../theme/palette";

const LIGHT_GREY_100 = "rgb(245, 245, 245)";

const XML = `<?xml version="1.0"?>
<ccp4i2>
  <ccp4i2_header><pluginName>cmapcoeff</pluginName></ccp4i2_header>
  <ccp4i2_body>
    <errorReportList>
      <errorReport>
        <className>cmapcoeff</className>
        <code>45</code>
        <details>Error in processing output files</details>
        <stack>Traceback (most recent call last):
  File "x.py", line 1
RuntimeError: boom</stack>
        <severity>ERROR</severity>
      </errorReport>
    </errorReportList>
  </ccp4i2_body>
</ccp4i2>`;

const lightTheme = createTheme({ palette: lightPaletteOptions });
const darkTheme = createTheme({ palette: darkPaletteOptions });

const renderPres = (theme: typeof lightTheme) => {
  const { container } = render(
    <ThemeProvider theme={theme}>
      <Diagnostic xmlDocument={XML} />
    </ThemeProvider>
  );
  const pres = Array.from(container.querySelectorAll("pre"));
  // Details + stack trace: both must be present or the assertions are vacuous.
  expect(pres).toHaveLength(2);
  return pres.map((pre) => getComputedStyle(pre));
};

describe("Diagnostic <pre> blocks", () => {
  it("dark mode: no light grey background, and text contrasts with it", () => {
    for (const style of renderPres(darkTheme)) {
      expect(style.backgroundColor).not.toBe(LIGHT_GREY_100);
      expect(style.backgroundColor).toBe("rgb(51, 51, 51)");
      expect(style.color).toBe("rgb(255, 255, 255)");
      expect(style.color).not.toBe(style.backgroundColor);
    }
  });

  it("light mode: light grey background under dark text", () => {
    for (const style of renderPres(lightTheme)) {
      expect(style.backgroundColor).toBe(LIGHT_GREY_100);
      expect(style.color).toBe("rgba(0, 0, 0, 0.87)");
      expect(style.color).not.toBe(style.backgroundColor);
    }
  });
});

describe("diagnosticSurfaceColors", () => {
  it.each([
    ["light", lightTheme],
    ["dark", darkTheme],
  ])("%s: error rows are a translucent tint, never a solid red", (_, theme) => {
    const { errorSummary } = diagnosticSurfaceColors(theme);
    for (const colour of Object.values(errorSummary)) {
      // alpha() yields rgba(...) with an alpha strictly below 1, so the row
      // keeps the paper colour underneath and text.primary stays legible.
      const match = colour.match(/^rgba\(\d+, \d+, \d+, (0\.\d+)\)$/);
      expect(match, colour).not.toBeNull();
      expect(Number(match![1])).toBeLessThan(0.5);
    }
    // The expanded row must read as distinct from the collapsed one.
    expect(errorSummary.expanded).not.toBe(errorSummary.idle);
  });

  it("swaps the pre surface between modes rather than fixing it light", () => {
    const light = diagnosticSurfaceColors(lightTheme).pre;
    const dark = diagnosticSurfaceColors(darkTheme).pre;
    expect(light.backgroundColor).not.toBe(dark.backgroundColor);
    expect(light.color).not.toBe(dark.color);
  });
});
