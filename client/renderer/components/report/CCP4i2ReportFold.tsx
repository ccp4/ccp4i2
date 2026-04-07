/*
 * Copyright (C) 2025 Newcastle University
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
import { PropsWithChildren } from "react";
import $ from "jquery";
import {
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Typography,
} from "@mui/material";
import { useTheme } from "@mui/material/styles";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import {
  CCP4i2ReportElement,
  CCP4i2ReportElementProps,
} from "./CCP4i2ReportElement";

// Extended interface that includes PropsWithChildren
interface CCP4i2ReportFoldProps
  extends PropsWithChildren<CCP4i2ReportElementProps> {}

export const CCP4i2ReportFold: React.FC<CCP4i2ReportFoldProps> = (props) => {
  const theme = useTheme();

  // Memoize the content processing to avoid recalculation
  const foldContent = useMemo(() => {
    if (!props.item) return [];

    const children = $(props.item).children().toArray();

    // Process floating elements
    children.forEach((child) => {
      const styleString = $(child).attr("style");
      if (styleString && styleString.includes("float:")) {
        const fixedStyle = styleString
          .replace(/float:\s*left;?/g, "")
          .replace(/float:\s*right;?/g, "");
        $(child).attr("style", fixedStyle);
      }
    });

    // Generate content
    return children.map((child, iChild) => (
      <CCP4i2ReportElement
        key={iChild}
        iItem={iChild}
        item={child}
        job={props.job}
      />
    ));
  }, [props.item, props.job]);

  return (
    <Accordion
      disableGutters
      defaultExpanded={$(props.item).attr("initiallyOpen") === "True"}
    >
      <AccordionSummary
        expandIcon={<ExpandMoreIcon />}
        aria-controls="panel-content"
        sx={{
          backgroundColor: theme.palette.action.hover,
          "&:hover": {
            backgroundColor: theme.palette.action.selected,
          },
        }}
      >
        <Typography variant="subtitle1" sx={{ fontWeight: 500 }}>
          {$(props.item).attr("label") || "Untitled Section"}
        </Typography>
      </AccordionSummary>

      <AccordionDetails sx={{ p: 2 }}>
        {foldContent}
        {props.children}
      </AccordionDetails>
    </Accordion>
  );
};
