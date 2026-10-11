import { useMemo } from "react";
import $ from "jquery";
import { Typography } from "@mui/material";
import type { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";
import { cssToDict } from "../../lib/report-css";

/**
 * A line of report text (`parent.addText(text=..., style=...)`).
 *
 * Rendered as a block with a little vertical space, in the body type, and
 * with the style the report asked for. It was a bare inline span that
 * dropped its style: consecutive lines ran into each other and the bold or
 * coloured text the Python reports write ("Now: ...", a red failure, an
 * orange warning) never showed as such.
 */
export const CCP4i2ReportText: React.FC<CCP4i2ReportElementProps> = (props) => {
  const { style, innerHTML } = useMemo(() => {
    const cssText = $(props.item).attr("style") || "";
    return {
      style: cssText ? cssToDict(cssText) : {},
      innerHTML: props.item?.innerHTML || "",
    };
  }, [props.item, props.job]);

  if (!innerHTML.trim()) return null;
  return (
    <Typography
      component="div"
      variant="body2"
      sx={{ my: 0.5, lineHeight: 1.5 }}
      style={style}
      dangerouslySetInnerHTML={{ __html: innerHTML }}
    />
  );
};
