import { useMemo, useCallback, useEffect, useState } from "react";
import $ from "jquery";
import { Box, Button, Typography } from "@mui/material";
import { Description, OpenInNew } from "@mui/icons-material";
import { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";
import { useFilePreviewContext } from "../../providers/file-preview-context";
import { apiBlob, apiPost } from "../../api-fetch";

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

  // For images, fetch via authenticated apiBlob and create a local blob URL.
  // This avoids <img src> making unauthenticated browser requests.
  const [blobUrl, setBlobUrl] = useState<string | null>(null);

  // Project-relative path of the file. Job files live under CCP4_JOBS/job_N/,
  // with dotted job numbers (sub-jobs) nesting one directory per component.
  const fullPath = useMemo(() => {
    if (!relativePath) return null;
    const jobDirSegments = props.job.number
      .split(".")
      .map((n: string) => `job_${n}`);
    return `CCP4_JOBS/${jobDirSegments.join("/")}/${relativePath}`;
  }, [relativePath, props.job]);

  // The directory the report lives in — the subtree a grant is scoped to,
  // and the one a multi-page report's relative links stay within.
  const reportDirectory = useMemo(
    () => (fullPath ? fullPath.replace(/\/[^/]*$/, "") : null),
    [fullPath]
  );

  const fileUrl = useMemo(() => {
    if (!projectId || !fullPath) return null;

    if (fileType === "html") {
      // Path-based URL so relative links in multi-page HTML reports (e.g. ProSMART)
      // resolve correctly against the directory hierarchy.
      return `/api/proxy/ccp4i2/projects/${projectId}/files_by_path/${fullPath}`;
    }
    // Query-parameter approach is fine for programmatic fetches (text, image).
    return `/api/proxy/ccp4i2/projects/${projectId}/project_file?path=${encodeURIComponent(fullPath)}`;
  }, [projectId, fullPath, fileType]);

  const handleClick = useCallback(async () => {
    if (!fileUrl) return;

    if (fileType === "html") {
      // HTML reports open in a new browser tab, which issues a plain browser
      // navigation with no Authorization header — as do every image, stylesheet
      // and relative fetch the page then makes for itself. Mint a read grant
      // for the report's own directory and seed the tab with it; the proxy
      // turns it into a cookie scoped to that directory, which the page's
      // subsequent requests inherit.
      try {
        const response = await apiPost<{ data?: { grant?: string } }>(
          `projects/${projectId}/file_grant`,
          { path: reportDirectory }
        );
        const grant = response?.data?.grant;
        if (grant) {
          window.open(
            `${fileUrl}?file_grant=${encodeURIComponent(grant)}`,
            "_blank"
          );
          return;
        }
      } catch (error) {
        console.error("Failed to mint file grant:", error);
      }
      // No grant: open anyway rather than doing nothing, so the failure shows
      // up as the backend's 401 rather than a dead button.
      window.open(fileUrl, "_blank");
    } else if (fileType === "image") {
      // Use the already-fetched blob URL for the preview dialog too
      if (!blobUrl) return;
      setContentSpecification({
        url: blobUrl,
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
  }, [
    fileUrl,
    fileType,
    relativePath,
    label,
    blobUrl,
    projectId,
    reportDirectory,
    setContentSpecification,
  ]);

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
