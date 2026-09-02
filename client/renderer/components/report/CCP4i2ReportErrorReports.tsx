import React, { useMemo } from "react";
import Diagnostic from "../diagnostic";
import { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";

/**
 * The `errorReportList` inside a report element, back as a string.
 *
 * `Diagnostic` parses a document, because that is what the Diagnostics tab
 * hands it — a file read off disk. Here the same payload arrives already
 * parsed, as part of the report tree, so it is serialised back out rather than
 * giving the parser a second entry point to keep in step. Exported for the
 * test: the round trip is the only fragile part of this component.
 */
export const errorReportListXml = (item: Element): string => {
  const list = item.getElementsByTagName("errorReportList")[0];
  return list ? new XMLSerializer().serializeToString(list) : "";
};

/**
 * A report that could not be rendered, drawn by the Diagnostics panel's own
 * renderer.
 *
 * When a report class raises, the server sends the failure as an
 * `errorReportList` — the same shape a job's `diagnostic.xml` carries — rather
 * than inventing a second vocabulary for it. So this element does no rendering
 * of its own: it hands the payload to `Diagnostic`, which already knows how to
 * draw a severity icon, a heading naming the task and the line, and a
 * traceback folded away.
 *
 * Two error reports can arrive together and mean different things: the report
 * class failing, and the job's own record of why it failed. Either happens
 * without the other — a report class can raise on a job that ran perfectly —
 * so both are shown, in that order.
 */
export const CCP4i2ReportErrorReports: React.FC<CCP4i2ReportElementProps> = ({
  item,
}) => {
  const xmlDocument = useMemo(() => errorReportListXml(item), [item]);

  if (!xmlDocument) return null;

  return <Diagnostic xmlDocument={xmlDocument} embedded />;
};

export default CCP4i2ReportErrorReports;
