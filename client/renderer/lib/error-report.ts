/**
 * A job's errors as plain text a user can paste into an email or a bug report.
 *
 * The Diagnostics panel draws a traceback folded inside an accordion, which is
 * readable but not copyable: a user who hit a UnicodeEncodeError in freerflag
 * (2026-10-22) had no easy way to send it on, and what they could select
 * lacked the versions and platform needed to make sense of it. This puts the
 * whole report, with the environment first, in one string.
 */
import { ErrorReport, HeaderInfo, describeReport } from "./diagnostic-parse";
import { PROXY_BASE, authHeaders } from "./staged-upload";

/** The fields of GET version/ a report names (VersionInfo in ccp4i2-api). */
export interface ServerEnvironment {
  version?: string;
  ccp4_version?: string | null;
  python_version?: string;
  platform?: string;
}

let serverEnvironment: Promise<ServerEnvironment | null> | null = null;

/** The server's versions and platform, fetched once per session. */
export function fetchServerEnvironment(): Promise<ServerEnvironment | null> {
  if (!serverEnvironment) {
    serverEnvironment = authHeaders()
      .then((headers) => fetch(`${PROXY_BASE}version/`, { headers }))
      .then((r) => (r.ok ? r.json() : null))
      .catch(() => null);
  }
  return serverEnvironment;
}

const HEADER_LABELS: [keyof HeaderInfo, string][] = [
  ["pluginName", "Task"],
  ["projectName", "Project"],
  ["jobNumber", "Job"],
  ["creationTime", "Created"],
  ["hostName", "Host"],
];

export function formatErrorReport(
  header: Partial<HeaderInfo>,
  errorReports: ErrorReport[],
  server: ServerEnvironment | null,
  userAgent: string
): string {
  const lines: string[] = ["CCP4i2 error report", ""];
  const field = (label: string, value?: string | null) => {
    if (value) lines.push(`${label}: ${value}`);
  };

  field("CCP4i2", server?.version);
  field("CCP4", server?.ccp4_version);
  field("Python", server?.python_version);
  field("Server OS", server?.platform);
  if (!server) lines.push("Server: version information unavailable");
  field("Client", userAgent);
  HEADER_LABELS.forEach(([key, label]) => field(label, header[key]));

  errorReports.forEach((report, index) => {
    lines.push("", `--- ${index + 1}. ${report.severity || "?"}: ${describeReport(report)}`);
    if (report.description) lines.push(report.description);
    if (report.details) lines.push("", report.details);
    if (report.stack) lines.push("", report.stack);
  });

  return lines.join("\n") + "\n";
}
