import { useEffect } from "react";
import { Job, Project } from "../types/models";
import { useJobDirectory } from "../utils";
import { useFileSystemFileBrowser } from "../providers/file-system-file-browser-context";
import DirectoryBrowser from "./directory-browser";
import { FileSystemFileMenu } from "./file-system-file-menu";

interface JobDirectoryViewProps {
  job: Job;
  project: Project;
}
export const JobDirectoryView: React.FC<JobDirectoryViewProps> = ({
  job,
  project,
}) => {
  // Lists the job's own folders, re-listing the open ones while it runs
  const { endpoint, pollInterval, revision } = useJobDirectory(job);

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
      <DirectoryBrowser endpoint={endpoint} pollInterval={pollInterval} revision={revision} />
      <FileSystemFileMenu onClose={handleMenuClose} />
    </>
  );
};
