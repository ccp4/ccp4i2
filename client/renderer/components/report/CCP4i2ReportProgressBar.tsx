import { useMemo } from "react";
import $ from "jquery";
import { Box, LinearProgress, Typography } from "@mui/material";
import { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";

/**
 * A styled progress bar report element.
 *
 * Reads `value`, `max`, and `label` attributes from the server-emitted
 * <CCP4i2ReportProgress> element and renders a labeled MUI LinearProgress.
 *
 * Used by the Python report system via `parent.addProgress(value=..., max=..., label=...)`.
 */
export const CCP4i2ReportProgressBar: React.FC<CCP4i2ReportElementProps> = ({
  item,
}) => {
  const { value, max, label, percent, color } = useMemo(() => {
    const el = $(item);
    const v = parseFloat(el.attr("value") || "0");
    const m = parseFloat(el.attr("max") || "100");
    const pct = m > 0 ? (v / m) * 100 : 0;

    // Color: green >= 90%, blue >= 70%, orange >= 50%, red < 50%
    let c: "success" | "primary" | "warning" | "error" = "success";
    if (pct < 50) c = "error";
    else if (pct < 70) c = "warning";
    else if (pct < 90) c = "primary";

    return {
      value: v,
      max: m,
      label: el.attr("label") || "",
      percent: pct,
      color: c,
    };
  }, [item]);

  return (
    <Box sx={{ display: "flex", alignItems: "center", gap: 1, my: 0.5 }}>
      {label && (
        <Typography variant="body2" sx={{ minWidth: 140, flexShrink: 0 }}>
          {label}
        </Typography>
      )}
      <Box sx={{ flexGrow: 1, minWidth: 100 }}>
        <LinearProgress
          variant="determinate"
          value={Math.min(percent, 100)}
          color={color}
          sx={{ height: 8, borderRadius: 4 }}
        />
      </Box>
      <Typography
        variant="body2"
        sx={{ minWidth: 80, textAlign: "right", flexShrink: 0 }}
      >
        {Math.round(value)}/{Math.round(max)} ({percent.toFixed(1)}%)
      </Typography>
    </Box>
  );
};
