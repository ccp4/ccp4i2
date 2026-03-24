import { useMemo, useCallback } from "react";
import $ from "jquery";
import { Box, Button, Typography } from "@mui/material";
import { Description, OpenInNew } from "@mui/icons-material";
import { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";
import { useFilePreviewContext } from "../../providers/file-preview-context";

/**
 * Detect Monaco language from a filename extension.
 */
function languageFromFilename(filename: string): string {
  const lower = filename.toLowerCase();
  if (lower.endsWith(".json")) return "json";
  if (lower.endsWith(".xml")) return "xml";
  if (lower.endsWith(".cif")) return "cif";
  if (lower.endsWith(".csv")) return "csv";
  return "text";
}

/**
 * Report element that renders a clickable link to a file in the job directory.
 *
 * For text/log files: opens the Monaco preview dialog.
 * For HTML files: opens the file in a new browser tab.
 * For images: renders inline in the report.
 *
 * Expected XML attributes:
 *   label        — display text (e.g. "Show Pointless logfile")
 *   relativePath — path relative to job dir (e.g. "job_1/log.txt")
 *   fileType     — "text" (Monaco), "html" (new tab), or "image" (inline)
 *   projectId    — project integer PK for project_file endpoint
 */
export const CCP4i2ReportFileLink: React.FC<CCP4i2ReportElementProps> = (
  props
) => {
  const { setContentSpecification } = useFilePreviewContext();

  const { label, relativePath, fileType, projectId } = useMemo(() => {
    const el = $(props.item);
    return {
      label: el.attr("label") || "View file",
      relativePath: el.attr("relativePath") || "",
      fileType: el.attr("fileType") || "text",
      projectId: el.attr("projectId") || "",
    };
  }, [props.item]);

  const fileUrl = useMemo(() => {
    if (!projectId || !relativePath) return null;
    // The project_file endpoint serves arbitrary files from the project dir.
    // Job files live under CCP4_JOBS/job_N/ inside the project directory.
    // We use the job number from props.job to build the full project-relative path.
    const jobDirSegments = props.job.number
      .split(".")
      .map((n: string) => `job_${n}`);
    const jobDirPath = `CCP4_JOBS/${jobDirSegments.join("/")}`;
    const fullPath = `${jobDirPath}/${relativePath}`;
    return `/api/proxy/ccp4i2/projects/${projectId}/project_file?path=${encodeURIComponent(fullPath)}`;
  }, [projectId, relativePath, props.job]);

  const handleClick = useCallback(() => {
    if (!fileUrl) return;

    if (fileType === "html") {
      // HTML reports open in a new browser tab
      window.open(fileUrl, "_blank");
    } else if (fileType === "image") {
      // Images can also be opened in the preview dialog at full size
      setContentSpecification({
        url: fileUrl,
        title: label,
        language: "image",
      });
    } else {
      // Text/log files open in the Monaco preview dialog
      const filename = relativePath.split("/").pop() || "file";
      setContentSpecification({
        url: fileUrl,
        title: label,
        language: languageFromFilename(filename),
      });
    }
  }, [fileUrl, fileType, relativePath, label, setContentSpecification]);

  if (!fileUrl) return null;

  // Images render inline in the report
  if (fileType === "image") {
    return (
      <Box sx={{ my: 1 }}>
        {label && (
          <Typography variant="subtitle2" sx={{ mb: 0.5 }}>
            {label}
          </Typography>
        )}
        <img
          src={fileUrl}
          alt={label}
          onClick={handleClick}
          style={{
            maxWidth: "100%",
            maxHeight: "500px",
            objectFit: "contain",
            cursor: "pointer",
            border: "1px solid #e0e0e0",
            borderRadius: "4px",
          }}
          title="Click to enlarge"
        />
      </Box>
    );
  }

  return (
    <Button
      variant="outlined"
      size="small"
      startIcon={fileType === "html" ? <OpenInNew /> : <Description />}
      onClick={handleClick}
      sx={{ m: 0.5 }}
    >
      {label}
    </Button>
  );
};
