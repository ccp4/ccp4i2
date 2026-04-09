import { useMemo } from "react";
import $ from "jquery";
import { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";
import { AlignmentViewer } from "../alignment-viewer";

export const CCP4i2ReportAlignment: React.FC<CCP4i2ReportElementProps> = ({
  item,
}) => {
  const alignmentText = useMemo(() => {
    // The alignment text is stored as the text content of the XML element
    return $(item).text() || "";
  }, [item]);

  if (!alignmentText.trim()) {
    return null;
  }

  return <AlignmentViewer alignment={alignmentText} />;
};
