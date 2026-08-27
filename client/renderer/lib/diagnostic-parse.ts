/**
 * Reading a job's diagnostic.xml.
 *
 * This is the one place the element names live. The backend's half of the same
 * contract is asserted in
 * `server/ccp4i2/tests/unit/validation/test_diagnostic_xml_contract.py`, and
 * the two tests name the same tags on purpose: change one side alone and one
 * of them fails.
 *
 * They had drifted. The backend has always written `<class>`; this parser read
 * `<className>`, so every error CCP4i2 ever reported drew an empty heading —
 * the traceback was there, in `<details>`, under a title reading " - 993".
 */

export interface ErrorReport {
  /** The task the report came from, e.g. "refmac". */
  className: string;
  code: string;
  /** Which parameter or object it is about, e.g. "XYZIN". */
  name: string;
  /** What this code means in general, from the task's ERROR_CODES. */
  description: string;
  severity: string;
  /** What went wrong on this occasion. */
  details: string;
  /** A traceback, when there is one. */
  stack: string;
}

export interface HeaderInfo {
  function: string;
  userId: string;
  hostName: string;
  creationTime: string;
  pluginName: string;
  projectName: string;
  projectId: string;
  jobId: string;
  jobNumber: string;
}

const HEADER_FIELDS: (keyof HeaderInfo)[] = [
  "function",
  "userId",
  "hostName",
  "creationTime",
  "pluginName",
  "projectName",
  "projectId",
  "jobId",
  "jobNumber",
];

export const parseDiagnosticXml = (
  xmlString: string
): { header: Partial<HeaderInfo>; errorReports: ErrorReport[] } => {
  const doc = new DOMParser().parseFromString(xmlString, "text/xml");

  const header: Partial<HeaderInfo> = {};
  const headerElement = doc.querySelector("ccp4i2_header");
  if (headerElement) {
    HEADER_FIELDS.forEach((field) => {
      const element = headerElement.querySelector(field);
      if (element) header[field] = element.textContent || "";
    });
  }

  const errorReports: ErrorReport[] = [];
  doc.querySelectorAll("errorReport").forEach((errorElement) => {
    const text = (selector: string) => {
      const value = errorElement.querySelector(selector)?.textContent?.trim() || "";
      // A field the backend left unset is serialised through Python's str(),
      // so it arrives as the word "None" and rendered as one --- the heading
      // read "SubstituteLigand \u00b7 None \u2014 210". Newer jobs write an empty
      // element instead; diagnostic.xml files already on disk do not, and
      // they are still opened.
      return value === "None" ? "" : value;
    };

    errorReports.push({
      // `class` is what the backend writes, and what every diagnostic.xml
      // already on disk contains. `className` is accepted too, so a future
      // producer using the other spelling is not silently ignored.
      className: text("className") || text("class"),
      code: text("code"),
      name: text("name"),
      description: text("description"),
      severity: text("severityName") || text("severity"),
      details: text("details"),
      stack: text("stack"),
    });
  });

  return { header, errorReports };
};

/** The heading for one report: the task, the parameter, and the code. */
export const describeReport = (report: ErrorReport): string => {
  const who = [report.className, report.name].filter(Boolean).join(" · ");
  return `${who || "Error"}${report.code ? ` — ${report.code}` : ""}`;
};
