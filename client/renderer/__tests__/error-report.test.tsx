/**
 * "Copy error report": a job's errors as text, with the versions and platform
 * that make them actionable. A freerflag UnicodeEncodeError (2026-10-22) could
 * be read in the Diagnostics panel but not easily copied out of it.
 */
import React from "react";
import { describe, it, expect, vi, afterEach } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import Diagnostic from "../components/diagnostic";
import { formatErrorReport } from "../lib/error-report";

const XML = `<ccp4i2><ccp4i2_header><pluginName>freerflag</pluginName>
<projectName>toxd</projectName><jobNumber>3</jobNumber></ccp4i2_header>
<errorReportList><errorReport><class>freerflag</class><code>39</code>
<severity>ERROR</severity><details>UnicodeEncodeError: 'charmap' codec</details>
<stack>Traceback (most recent call last):\n  File "freerflag.py", line 354</stack>
</errorReport></errorReportList></ccp4i2>`;

const SERVER = {
  version: "3.1.0a104",
  ccp4_version: "10.0.000",
  python_version: "3.11.9",
  platform: "Windows-10-10.0.26200-SP0",
};

describe("formatErrorReport", () => {
  it("puts the environment first, then every error with its traceback", () => {
    const text = formatErrorReport(
      { pluginName: "freerflag", projectName: "toxd", jobNumber: "3" },
      [{ className: "freerflag", code: "39", name: "", description: "",
         severity: "ERROR", details: "UnicodeEncodeError", stack: "Traceback ..." }],
      SERVER,
      "Mozilla/5.0 Electron/37"
    );
    expect(text).toContain("CCP4i2: 3.1.0a104");
    expect(text).toContain("CCP4: 10.0.000");
    expect(text).toContain("Server OS: Windows-10-10.0.26200-SP0");
    expect(text).toContain("Client: Mozilla/5.0 Electron/37");
    expect(text).toContain("Task: freerflag");
    expect(text).toContain("--- 1. ERROR: freerflag — 39");
    expect(text.indexOf("CCP4i2:")).toBeLessThan(text.indexOf("UnicodeEncodeError"));
    expect(text).toContain("Traceback ...");
  });

  it("still reports the errors when the server cannot be asked", () => {
    const text = formatErrorReport({}, [], null, "ua");
    expect(text).toContain("version information unavailable");
  });
});

describe("Diagnostic's Copy error report button", () => {
  afterEach(() => vi.restoreAllMocks());

  it("copies the report, traceback included", async () => {
    vi.spyOn(globalThis, "fetch").mockResolvedValue(
      new Response(JSON.stringify(SERVER), { status: 200 })
    );
    const writeText = vi.fn().mockResolvedValue(undefined);
    Object.defineProperty(navigator, "clipboard", { value: { writeText }, configurable: true });

    render(<Diagnostic xmlDocument={XML} embedded />);
    fireEvent.click(screen.getByRole("button", { name: /copy error report/i }));

    await waitFor(() => expect(writeText).toHaveBeenCalledOnce());
    const text = writeText.mock.calls[0][0] as string;
    expect(text).toContain("line 354");
    expect(text).toContain("Project: toxd");
    await screen.findByText("Copied");
  });
});
