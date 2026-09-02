// @vitest-environment jsdom
import { describe, it, expect } from "vitest";
import { errorReportListXml } from "../components/report/CCP4i2ReportErrorReports";
import { parseDiagnosticXml } from "../lib/diagnostic-parse";

/**
 * A report that could not be rendered arrives as an `errorReportList` — the
 * same payload a job's diagnostic.xml carries — embedded in the report tree.
 * This is the seam: it has to come back out of that tree and through the
 * Diagnostics parser with every field intact, or the panel draws empty
 * headings, which is exactly the bug the contract test exists to prevent.
 */
const REPORT_ELEMENT = `<CCP4i2ReportErrorReports key="Failure_1" class="" style="">
  <errorReportList>
    <errorReport>
      <class>cmapcoeff report</class>
      <code>PROGRAM_XML_PARSE_ERROR</code>
      <details>File: /projects/TcECH/CCP4_JOBS/job_25/program.xml</details>
      <name>ccp4i2/lib/utils/reporting/i2_report.py:612 in generate_job_report()</name>
      <description>The task's program.xml is not well-formed XML.</description>
      <stack>Traceback (most recent call last):
  File "i2_report.py", line 612
ParseError: no element found</stack>
      <severity>4</severity>
      <severityName>ERROR</severityName>
    </errorReport>
    <errorReport>
      <class>cmapcoeff</class>
      <code>45</code>
      <details>Error in processing output files</details>
      <severityName>ERROR</severityName>
    </errorReport>
  </errorReportList>
</CCP4i2ReportErrorReports>`;

const element = (xml: string): Element =>
  new DOMParser().parseFromString(xml, "text/xml").documentElement;

describe("errorReportListXml", () => {
  it("brings every field back out through the Diagnostics parser", () => {
    const { errorReports } = parseDiagnosticXml(
      errorReportListXml(element(REPORT_ELEMENT))
    );

    expect(errorReports).toHaveLength(2);
    const [failure] = errorReports;
    expect(failure.className).toBe("cmapcoeff report");
    expect(failure.code).toBe("PROGRAM_XML_PARSE_ERROR");
    expect(failure.name).toContain("i2_report.py:612");
    expect(failure.description).toContain("not well-formed");
    expect(failure.details).toContain("program.xml");
    expect(failure.stack).toContain("ParseError");
    expect(failure.severity).toBe("ERROR");
  });

  it("keeps the rendering failure ahead of the job's own report", () => {
    const { errorReports } = parseDiagnosticXml(
      errorReportListXml(element(REPORT_ELEMENT))
    );
    expect(errorReports[0].code).toBe("PROGRAM_XML_PARSE_ERROR");
    expect(errorReports[1].code).toBe("45");
  });

  it("returns nothing when the element carries no list", () => {
    expect(errorReportListXml(element("<CCP4i2ReportErrorReports />"))).toBe("");
  });
});
