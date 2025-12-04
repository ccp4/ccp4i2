import { useMemo, useEffect, useContext } from "react";
import { Job, Project } from "../types/models";
import { useProject } from "../utils";
import { Box } from "@mui/material";
import DirectoryBrowser from "./directory-browser";
import { FileSystemFileMenu } from "./file-system-file-menu";
import { useFileSystemFileBrowser } from "../providers/file-system-file-browser-context";
import { useFilePreviewContext } from "../providers/file-preview-context";

interface JobLogViewerProps {
  job: Job;
  project: Project;
}

export const JobLogViewer: React.FC<JobLogViewerProps> = ({ job, project }) => {
  const { directory } = useProject(project.id);
  const { closeMenu } = useFileSystemFileBrowser();
  const { setContentSpecification } = useFilePreviewContext();

  const directoryData = useMemo(() => {
    if (!directory || !job || !directory.container) {
      return null;
    }

    let dirNode = directory.container.find(
      (item: any) => item.name === "CCP4_JOBS"
    );

    if (!dirNode) {
      return null;
    }

    const jobNumberElements = job.number.split(".").reverse();
    while (jobNumberElements.length > 0) {
      const jobNumber = jobNumberElements.pop();
      dirNode = dirNode.contents?.find(
        (item: any) => item.name === `job_${jobNumber}`
      );
      if (!dirNode) {
        return null;
      }
      if (jobNumberElements.length === 0) {
        return dirNode.contents || [];
      }
    }

    return null;
  }, [job, directory]);

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

  // Filter to show only log files and relevant directories
  const logFileFilter = (item: any) => {
    if (item.type === "directory") {
      return true; // Always show directories for navigation
    }

    // Show common log file extensions
    const logExtensions = [".log", ".out", ".err", ".txt", ".xml", ".json"];

    //Return draft parameter files
    const excludePattern = [
      /^input_params\.previous_\d+\.xml$/i,
      /^params\.previous_\d+\.xml$/i,
    ];
    if (excludePattern.some((pattern) => pattern.test(item.name))) {
      return false;
    }

    return logExtensions.some((ext) => item.name.toLowerCase().endsWith(ext));
  };

  const handleFileSelect = (menuNode: any) => {
    // Handle file selection if needed, e.g., preview or download
    if (menuNode.type !== "directory") {
      // Set the preview node or perform any action needed
      const composite_path = `/api/proxy/projects/${
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
      {directoryData ? (
        <>
          <DirectoryBrowser
            directoryTree={directoryData || []}
            title={`Job ${job.number} Files`}
            width="100%"
            height="100%"
            fileFilter={logFileFilter}
            showSearch={true}
            showFileSizes={true}
            onItemClick={handleFileSelect}
          />
          <FileSystemFileMenu onClose={handleMenuClose} />
        </>
      ) : (
        <Box sx={{ p: 2, textAlign: "center", width: "100%" }}>
          No job files found for job {job.number}
        </Box>
      )}
    </Box>
  );
};
