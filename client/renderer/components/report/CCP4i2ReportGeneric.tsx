import { useMemo } from "react";
import $ from "jquery";
import { CCP4i2RVAPITable } from "./CCP4i2RVAPITable";
import { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";
import { CCP4i2ReportSVG } from "./CCP4i2ReportSVG";

export const CCP4i2ReportGeneric: React.FC<CCP4i2ReportElementProps> = (
  props
) => {
  // Check if the item is an RVAPI table
  // An RVAPI table is identified by the class "rvapi-page
  const isRVAPITable = useMemo(() => {
    return $(props.item).children("table.rvapi-page").length > 0;
  }, [props.item, props.job]);

  // Check if the item is an svg element
  const isSVG = useMemo(() => {
    // Look for any child whose tagName is 'svg' or ends with ':svg' (case-insensitive)
    const children = props.item.children;
    for (let i = 0; i < children.length; i++) {
      const tag = children[i].tagName?.toLowerCase();
      if (tag === "svg" || tag.endsWith(":svg")) {
        return true;
      }
    }
    return false;
  }, [props.item]);

  const tableBody = useMemo(() => {
    const tableBody = $(props.item).find("tbody").get(0);
    return tableBody;
  }, [props.item]);

  // Clean up SVG namespaces for rendering
  const cleanedInnerHTML = useMemo(() => {
    let html = props.item.innerHTML;
    // Remove all namespace prefixes like ns0:
    html = html.replace(/<\s*\/?\s*ns\d*:/g, (match) =>
      match.replace(/ns\d*:/, "")
    );
    // Remove xmlns:ns0="..." attributes
    html = html.replace(/xmlns:ns\d+="[^"]*"/g, "");

    return html;
  }, [props.item.innerHTML]);

  return isRVAPITable && tableBody ? (
    <CCP4i2RVAPITable iItem={props.iItem} item={tableBody} job={props.job} />
  ) : isSVG ? (
    <CCP4i2ReportSVG {...props} />
  ) : (
    <div dangerouslySetInnerHTML={{ __html: props.item.innerHTML }} />
  );
};
