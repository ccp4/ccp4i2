/*
 * Copyright (C) 2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
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
