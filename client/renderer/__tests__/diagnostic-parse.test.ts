// @vitest-environment jsdom
import { describe, it, expect } from "vitest";
import {
  describeReport,
  parseDiagnosticXml,
} from "../lib/diagnostic-parse";

/**
 * The panel's half of the diagnostic.xml contract. The backend's half is
 * server/ccp4i2/tests/unit/validation/test_diagnostic_xml_contract.py, and the
 * two name the same tags on purpose.
 */

// Exactly what the backend writes -- taken from a real failing job.
const REAL_DIAGNOSTIC = `<errorReportList>
  <errorReport>
    <class>refmac</class>
    <code>201</code>
    <details>Exit code: 5</details>
    <name>XYZIN</name>
    <description>Ligand has no geometry dictionary</description>
    <stack>Traceback (most recent call last):
  File "refmac.py", line 189</stack>
    <severity>4</severity>
    <severityName>ERROR</severityName>
  </errorReport>
</errorReportList>`;

// A diagnostic.xml written before description and stack existed. Every
// finished job on disk looks like this, and must still render.
const HISTORICAL_DIAGNOSTIC = `<errorReportList>
  <errorReport>
    <class>chainsaw</class>
    <code>993</code>
    <details>Python exception in processOutputFiles()</details>
    <name>processOutputFiles</name>
    <severity>4</severity>
    <severityName>ERROR</severityName>
  </errorReport>
</errorReportList>`;

describe("parseDiagnosticXml", () => {
  it("reads every field the backend writes", () => {
    const { errorReports } = parseDiagnosticXml(REAL_DIAGNOSTIC);
    expect(errorReports).toHaveLength(1);
    const report = errorReports[0];
    expect(report.className).toBe("refmac");
    expect(report.code).toBe("201");
    expect(report.name).toBe("XYZIN");
    expect(report.description).toBe("Ligand has no geometry dictionary");
    expect(report.details).toBe("Exit code: 5");
    expect(report.stack).toContain("Traceback");
    expect(report.severity).toBe("ERROR");
  });

  it("names the task from <class>, which is what is actually written", () => {
    // The bug: this parser read <className>, so the heading was always empty.
    const { errorReports } = parseDiagnosticXml(HISTORICAL_DIAGNOSTIC);
    expect(errorReports[0].className).toBe("chainsaw");
  });

  it("renders a job written before description and stack existed", () => {
    const { errorReports } = parseDiagnosticXml(HISTORICAL_DIAGNOSTIC);
    const report = errorReports[0];
    expect(report.description).toBe("");
    expect(report.stack).toBe("");
    expect(report.details).toContain("processOutputFiles");
  });

  it("accepts <className> too, so neither spelling is silently dropped", () => {
    const { errorReports } = parseDiagnosticXml(
      `<errorReportList><errorReport><className>servalcat</className></errorReport></errorReportList>`
    );
    expect(errorReports[0].className).toBe("servalcat");
  });

  it("reads the header block", () => {
    const { header } = parseDiagnosticXml(
      `<errorReportList><ccp4i2_header><pluginName>refmac</pluginName><jobNumber>3</jobNumber></ccp4i2_header></errorReportList>`
    );
    expect(header.pluginName).toBe("refmac");
    expect(header.jobNumber).toBe("3");
  });
});

describe("describeReport", () => {
  it("names the task and the parameter", () => {
    const { errorReports } = parseDiagnosticXml(REAL_DIAGNOSTIC);
    expect(describeReport(errorReports[0])).toBe("refmac · XYZIN — 201");
  });

  it("never renders a heading of punctuation alone", () => {
    // What the panel showed for years: " - 993".
    expect(
      describeReport({
        className: "",
        code: "993",
        name: "",
        description: "",
        severity: "ERROR",
        details: "",
        stack: "",
      })
    ).toBe("Error — 993");
  });
});
