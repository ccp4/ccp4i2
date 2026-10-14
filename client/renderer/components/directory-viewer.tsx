"use client";
import { useEffect } from "react";
import DirectoryBrowser from "./directory-browser";
import { useFileSystemFileBrowser } from "../providers/file-system-file-browser-context";
import { FileSystemFileMenu } from "./file-system-file-menu";

/** Polling interval for re-listing the open folders (10 seconds) */
const DIRECTORY_POLL_INTERVAL = 10000;

interface CCP4i2DirectoryViewerProps {
  projectId: number;
  /** Enable periodic polling to refresh directory contents (default: true) */
  enablePolling?: boolean;
}

export const CCP4i2DirectoryViewer: React.FC<CCP4i2DirectoryViewerProps> = ({
  projectId,
  enablePolling = true,
}) => {
  const { closeMenu } = useFileSystemFileBrowser();

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

  return (
    <>
      <DirectoryBrowser
        endpoint={`projects/${projectId}/directory`}
        pollInterval={enablePolling ? DIRECTORY_POLL_INTERVAL : 0}
      />
      <FileSystemFileMenu onClose={handleMenuClose} />
    </>
  );
};
