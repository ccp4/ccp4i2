import { useEffect, useMemo } from "react";
import { Job, Project } from "../types/models";
import { useJobDirectory } from "../utils";
import { useFileSystemFileBrowser } from "../providers/file-system-file-browser-context";
import DirectoryBrowser, { FileSystemItem } from "./directory-browser";
import { FileSystemFileMenu } from "./file-system-file-menu";
import { LinearProgress } from "@mui/material";

interface JobDirectoryViewProps {
  job: Job;
  project: Project;
}
export const JobDirectoryView: React.FC<JobDirectoryViewProps> = ({
  job,
  project,
}) => {
  // Use job-aware directory hook for adaptive polling based on job status
  const { directory } = useJobDirectory(project.id, job);

  const { closeMenu } = useFileSystemFileBrowser();

  // Clean up virtual anchor when component unmounts or menu closes
  const handleMenuClose = () => {
    const existing = document.getElementById("file-menu-anchor");
    if (existing && document.body.contains(existing)) {
      document.body.removeChild(existing);
    }

    closeMenu();
  };

  // Clean up on unmount
  useEffect(() => {
    return () => {
      const existing = document.getElementById("file-menu-anchor");
      if (existing && document.body.contains(existing)) {
        document.body.removeChild(existing);
      }
    };
  }, []);

  const directoryData = useMemo(() => {
    console.log("directory.container:", directory.container);
    if (!directory || !job || !directory.container) {
      console.log("Missing directory, job, or container");
      return null;
    }
    let dirNode = directory.container.find(
      (item: any) => item.name === "CCP4_JOBS"
    );
    if (!dirNode) {
      console.log("CCP4_JOBS directory not found");
      return [];
    }
    const jobNumberElements = job.number.split(".").reverse();
    console.log("Job number elements:", jobNumberElements);
    let cumulativePath: string = dirNode.path;
    while (jobNumberElements.length > 0) {
      const jobNumber = jobNumberElements.pop();
      console.log(`Looking for job_${jobNumber} in:`, dirNode);
      if (!dirNode.contents || !Array.isArray(dirNode.contents)) {
        console.error(`dirNode.contents is not an array:`, dirNode.contents);
        return [];
      }
      dirNode = dirNode.contents.find(
        (item: any) => item.name === `job_${jobNumber}`
      );
      cumulativePath += `/job_${jobNumber}`;
      console.log({ cumulativePath, dirNode });
      if (!dirNode) {
        console.error(`job_${jobNumber} not found in ${cumulativePath}`);
        return null;
      }
      if (jobNumberElements.length === 0) {
        console.log("Final dirNode.contents:", dirNode.contents);
        return dirNode.contents || [];
      }
    }
  }, [job, project, directory]);

  return directory ? (
    <>
      <DirectoryBrowser directoryTree={directoryData || []} />
      <FileSystemFileMenu onClose={handleMenuClose} />
    </>
  ) : (
    <LinearProgress />
  );
};
