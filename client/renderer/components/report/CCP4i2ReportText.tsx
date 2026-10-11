import { useMemo } from "react";
import $ from "jquery";
import { Typography } from "@mui/material";
import type { SxProps, Theme } from "@mui/material";
import type { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";
import { cssToDict } from "../../lib/report-css";

/**
 * A line of report text (`parent.addText(text=..., tone=...)`).
 *
 * A report says what a line is for, not how to paint it: `tone` is one of
 * a closed vocabulary and the theme decides the rest (so dark mode and the
 * house style follow). The html-era `style="..."` strings are still
 * honoured for the long tail that has not been converted.
 */
export type ReportTone =
  | "heading"
  | "emphasis"
  | "info"
  | "success"
  | "warning"
  | "error"
  | "muted"
  | "code";

const TONE_SX: Record<ReportTone, SxProps<Theme>> = {
  heading: { fontSize: "1.05rem", fontWeight: 600, mt: 1.5, mb: 0.5 },
  emphasis: { fontWeight: 600 },
  info: { color: "info.main" },
  success: { color: "success.main" },
  warning: { color: "warning.dark" },
  error: { color: "error.main" },
  muted: { color: "text.secondary", fontStyle: "italic" },
  code: { fontFamily: "monospace", fontSize: "0.85em", whiteSpace: "pre-wrap" },
};

export const toneSx = (tone: string | undefined): SxProps<Theme> =>
  tone && tone in TONE_SX ? TONE_SX[tone as ReportTone] : {};

export const CCP4i2ReportText: React.FC<CCP4i2ReportElementProps> = (props) => {
  const { tone, style, innerHTML } = useMemo(() => {
    const el = $(props.item);
    const cssText = el.attr("style") || "";
    return {
      tone: el.attr("tone") || undefined,
      style: cssText ? cssToDict(cssText) : {},
      innerHTML: props.item?.innerHTML || "",
    };
  }, [props.item, props.job]);

  if (!innerHTML.trim()) return null;
  return (
    <Typography
      component="div"
      variant="body2"
      data-tone={tone}
      sx={{ my: 0.5, lineHeight: 1.5, ...(toneSx(tone) as object) }}
      style={style}
      dangerouslySetInnerHTML={{ __html: innerHTML }}
    />
  );
};
