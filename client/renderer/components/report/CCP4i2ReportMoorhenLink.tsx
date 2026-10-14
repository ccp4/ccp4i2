import { useMemo } from "react";
import $ from "jquery";
import { Button } from "@mui/material";
import { ViewInAr } from "@mui/icons-material";
import { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";
import { sceneSourceUrl } from "../../lib/scene-source";

/**
 * Report element: a button that opens Moorhen on a scene the job builds
 * when asked (server: report/actions.py MoorhenLink).
 *
 * Expected XML attributes:
 *   label     — display text (e.g. "xtal-0002 event 1 (score 0.93)")
 *   method    — the job's plugin method that builds the scene ("eventScene")
 *   kwargs    — its keyword arguments, as JSON
 *   projectId — project integer PK, for the viewer page's route
 */
export const CCP4i2ReportMoorhenLink: React.FC<CCP4i2ReportElementProps> = (props) => {
  const { label, href } = useMemo(() => {
    const el = $(props.item);
    const projectId = parseInt(el.attr("projectId") || `${props.job.project}`, 10);
    let kwargs: Record<string, unknown> = {};
    try {
      kwargs = JSON.parse(el.attr("kwargs") || "{}");
    } catch {
      kwargs = {};
    }
    return {
      label: el.attr("label") || "Open in Moorhen",
      href: sceneSourceUrl(projectId, {
        jobId: props.job.id,
        method: el.attr("method") || "",
        kwargs,
      }),
    };
  }, [props.item, props.job]);

  return (
    <Button
      size="small"
      variant="outlined"
      startIcon={<ViewInAr fontSize="small" />}
      sx={{ m: 0.5, textTransform: "none" }}
      disabled={!href}
      onClick={() => href && window.open(href, "_blank", "noopener,noreferrer")}
    >
      {label}
    </Button>
  );
};
