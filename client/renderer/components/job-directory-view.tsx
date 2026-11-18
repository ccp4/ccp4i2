import { useEffect } from "react";
import { Job, Project } from "../types/models";
import { useJobDirectory } from "../utils";
import { useFileSystemFileBrowser } from "../providers/file-system-file-browser-context";
import DirectoryBrowser, { FileSystemItem } from "./directory-browser";
import { FileSystemFileMenu } from "./file-system-file-menu";
import { useJobFilePreview } from "../hooks/use-job-file-preview";
import { isLogFile } from "../lib/file-kinds";

interface JobDirectoryViewProps {
  job: Job;
  project: Project;
}
export const JobDirectoryView: React.FC<JobDirectoryViewProps> = ({
  job,
  project,
}) => {
  // Lists the job's own folders, re-listing the open ones while it runs
  const { endpoint, pollInterval, revision, isJobActive } = useJobDirectory(job);
  const openFile = useJobFilePreview(job, project, isJobActive);

  // A log opens at once, as in the Logs tab, live while the job runs; any
  // other file opens its menu (#788).
  const openLogs = (item: FileSystemItem) => {
    if (!isLogFile(item.name)) return false;
    openFile(item);
    return true;
  };

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

  return (
    <>
      <DirectoryBrowser
        endpoint={endpoint}
        pollInterval={pollInterval}
        revision={revision}
        onFileClick={openLogs}
      />
      <FileSystemFileMenu onClose={handleMenuClose} livePreviewJob={isJobActive ? job : undefined} />
    </>
  );
};
