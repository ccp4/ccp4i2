import { useEffect, useCallback } from "react";
import { Job, Project } from "../types/models";
import { useJobDirectory } from "../utils";
import { Box } from "@mui/material";
import DirectoryBrowser, { FileSystemItem } from "./directory-browser";
import { FileSystemFileMenu } from "./file-system-file-menu";
import { useFileSystemFileBrowser } from "../providers/file-system-file-browser-context";
import { useFilePreviewContext } from "../providers/file-preview-context";

interface JobLogViewerProps {
  job: Job;
  project: Project;
}

export const JobLogViewer: React.FC<JobLogViewerProps> = ({ job, project }) => {
  // Lists the job's own folders, re-listing the open ones while it runs
  const { endpoint, pollInterval, revision } = useJobDirectory(job);
  const { closeMenu } = useFileSystemFileBrowser();
  const { setContentSpecification } = useFilePreviewContext();

  // Clean up virtual anchor when component unmounts
  useEffect(() => {
    return () => {
      const existing = document.getElementById("file-menu-anchor");
      if (existing && document.body.contains(existing)) {
        document.body.removeChild(existing);
      }
    };
  }, []);

  // Handle cleanup when menu closes
  const handleMenuClose = () => {
    const existing = document.getElementById("file-menu-anchor");
    if (existing && document.body.contains(existing)) {
      document.body.removeChild(existing);
    }
    closeMenu();
  };

  // Filter to show only log files (folders are always shown)
  const logFileFilter = useCallback((item: FileSystemItem) => {
    if (item.type === "directory") {
      return true; // Always show directories for navigation
    }

    // Show common log file extensions
    const logExtensions = [".log", ".out", ".err", ".txt", ".xml", ".json", ".inp"];

    //Return draft parameter files
    const excludePattern = [
      /^input_params\.previous_\d+\.xml$/i,
      /^params\.previous_\d+\.xml$/i,
    ];
    if (excludePattern.some((pattern) => pattern.test(item.name))) {
      return false;
    }

    return logExtensions.some((ext) => item.name.toLowerCase().endsWith(ext));
  }, []);

  const handleFileSelect = (menuNode: FileSystemItem) => {
    // Handle file selection if needed, e.g., preview or download
    if (menuNode.type !== "directory") {
      // Set the preview node or perform any action needed
      const composite_path = `/api/proxy/ccp4i2/projects/${
        project.id
      }/project_file?path=${encodeURIComponent(
        menuNode?.path.slice(project.directory.length + 1) || ""
      )}`;

      setContentSpecification({
        url: composite_path,
        title: menuNode.name || "Preview",
        language: menuNode.name.endsWith(".json")
          ? "json"
          : menuNode.name.endsWith(".xml")
          ? "xml"
          : "text",
      });
    }
  };

  return (
    <Box sx={{ display: "flex", height: "100%", width: "100%" }}>
      <DirectoryBrowser
        endpoint={endpoint}
        pollInterval={pollInterval}
        revision={revision}
        fileFilter={logFileFilter}
        onItemClick={handleFileSelect}
      />
      <FileSystemFileMenu onClose={handleMenuClose} />
    </Box>
  );
};
