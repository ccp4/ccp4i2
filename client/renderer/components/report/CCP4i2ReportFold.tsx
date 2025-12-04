import { useEffect, useState, useMemo } from "react";
import { PropsWithChildren } from "react";
import $ from "jquery";
import {
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Grid2,
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
  const [expanded, setExpanded] = useState(
    $(props.item).attr("initiallyOpen") === "True"
  );

  // Memoize the content processing to avoid recalculation
  const { foldContent, nFloatingChildren } = useMemo(() => {
    if (!props.item) return { foldContent: [], nFloatingChildren: 0 };

    let floatingCount = 0;
    const children = $(props.item).children().toArray();

    // Process floating elements
    children.forEach((child) => {
      const styleString = $(child).attr("style");
      if (styleString && styleString.includes("float:")) {
        const fixedStyle = styleString
          .replace(/float:\s*left;?/g, "")
          .replace(/float:\s*right;?/g, "");
        $(child).attr("style", fixedStyle);
        floatingCount++;
      }
    });

    // Generate content
    const content = children.map((child, iChild) => (
      <CCP4i2ReportElement
        key={iChild}
        iItem={iChild}
        item={child}
        job={props.job}
      />
    ));

    return { foldContent: content, nFloatingChildren: floatingCount };
  }, [props.item, props.job]);

  const handleAccordionChange = (
    _event: React.SyntheticEvent,
    isExpanded: boolean
  ) => {
    setExpanded(isExpanded);
  };

  const renderContent = () => {
    if (false && nFloatingChildren > 0) {
      return (
        <Grid2 container spacing={1}>
          {foldContent.map((content, index) => (
            <Grid2 key={index} size={{ xs: 12 / nFloatingChildren }}>
              {content}
            </Grid2>
          ))}
        </Grid2>
      );
    }
    return foldContent;
  };

  return (
    <Accordion
      expanded={expanded}
      onChange={handleAccordionChange}
      disableGutters
      elevation={1}
    >
      <AccordionSummary
        expandIcon={<ExpandMoreIcon />}
        aria-controls="panel-content"
        id="panel-header"
        sx={{
          backgroundColor: theme.palette.action.hover,
          "&:hover": {
            backgroundColor: theme.palette.action.selected,
          },
          minHeight: 48,
          "& .MuiAccordionSummary-content": {
            margin: "8px 0",
          },
          "& .MuiAccordionSummary-content.Mui-expanded": {
            margin: "8px 0",
          },
        }}
      >
        <Typography variant="subtitle1" sx={{ fontWeight: 500 }}>
          {$(props.item).attr("label") || "Untitled Section"}
        </Typography>
      </AccordionSummary>

      <AccordionDetails sx={{ p: 2 }}>
        {renderContent()}
        {props.children}
      </AccordionDetails>
    </Accordion>
  );
};
