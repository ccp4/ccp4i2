import { useMemo, useCallback, useEffect, useState } from "react";
import $ from "jquery";
import { Box, Button, Typography } from "@mui/material";
import { Description, OpenInNew } from "@mui/icons-material";
import { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";
import { useFilePreviewContext } from "../../providers/file-preview-context";
import { apiBlob } from "../../api-fetch";

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
    // Job files live under CCP4_JOBS/job_N/ inside the project directory.
    // We use the job number from props.job to build the full project-relative path.
    const jobDirSegments = props.job.number
      .split(".")
      .map((n: string) => `job_${n}`);
    const jobDirPath = `CCP4_JOBS/${jobDirSegments.join("/")}`;
    const fullPath = `${jobDirPath}/${relativePath}`;

    if (fileType === "html") {
      // Path-based URL so relative links in multi-page HTML reports (e.g. ProSMART)
      // resolve correctly against the directory hierarchy.
      return `/api/proxy/ccp4i2/projects/${projectId}/files_by_path/${fullPath}`;
    }
    // Query-parameter approach is fine for programmatic fetches (text, image).
    return `/api/proxy/ccp4i2/projects/${projectId}/project_file?path=${encodeURIComponent(fullPath)}`;
  }, [projectId, relativePath, fileType, props.job]);

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

  // For images, fetch via authenticated apiBlob and create a local blob URL.
  // This avoids <img src> making unauthenticated browser requests.
  const [blobUrl, setBlobUrl] = useState<string | null>(null);
  useEffect(() => {
    if (fileType !== "image" || !fileUrl) return;
    let revoked = false;
    apiBlob(fileUrl).then((blob) => {
      if (!revoked) setBlobUrl(URL.createObjectURL(blob));
    }).catch(() => {});
    return () => {
      revoked = true;
      setBlobUrl((prev) => { if (prev) URL.revokeObjectURL(prev); return null; });
    };
  }, [fileUrl, fileType]);

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
        {blobUrl ? (
          <img
            src={blobUrl}
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
        ) : (
          <Typography variant="body2" color="text.secondary">Loading image…</Typography>
        )}
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
