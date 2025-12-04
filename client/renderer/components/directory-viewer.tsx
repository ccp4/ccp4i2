"use client";
import { LinearProgress } from "@mui/material";
import { useProject } from "../utils";
import { useEffect, useMemo, useState } from "react";
import DirectoryBrowser from "./directory-browser";
import { useFileSystemFileBrowser } from "../providers/file-system-file-browser-context";
import { FileSystemFileMenu } from "./file-system-file-menu";

interface CCP4i2DirectoryViewerProps {
  projectId: number;
  searchTerm?: string;
}

// Utility function to find all paths that contain matching items
const findMatchingPaths = (
  directoryTree: any[],
  searchTerm: string,
  currentPath: string = ""
): Set<string> => {
  const matchingPaths = new Set<string>();

  const searchRecursively = (items: any[], path: string) => {
    items.forEach((item) => {
      const itemPath = path ? `${path}/${item.name}` : item.name;

      // Check if current item matches search
      if (item.name.toLowerCase().includes(searchTerm.toLowerCase())) {
        // Add all parent paths to ensure they're expanded
        const pathParts = itemPath.split("/");
        for (let i = 0; i < pathParts.length; i++) {
          const parentPath = pathParts.slice(0, i + 1).join("/");
          matchingPaths.add(parentPath);
        }
      }

      // Recursively search children
      if (item.children && item.children.length > 0) {
        searchRecursively(item.children, itemPath);
      }
    });
  };

  searchRecursively(directoryTree, currentPath);
  return matchingPaths;
};

// Hook for managing expanded states based on search
const useSearchExpansion = (directoryTree: any[], searchTerm?: string) => {
  const [expandedPaths, setExpandedPaths] = useState<Set<string>>(new Set());

  useEffect(() => {
    if (searchTerm && searchTerm.trim()) {
      // Find all paths that should be expanded to show matches
      const pathsToExpand = findMatchingPaths(directoryTree, searchTerm.trim());
      setExpandedPaths(pathsToExpand);
    } else {
      // Clear expansions when search is cleared
      setExpandedPaths(new Set());
    }
  }, [directoryTree, searchTerm]);

  return {
    expandedPaths,
    shouldAutoExpand: Boolean(searchTerm && searchTerm.trim()),
  };
};

export const CCP4i2DirectoryViewer: React.FC<CCP4i2DirectoryViewerProps> = ({
  projectId,
  searchTerm,
}) => {
  const { directory } = useProject(projectId);
  const { closeMenu } = useFileSystemFileBrowser();

  // Use search expansion hook
  const { expandedPaths, shouldAutoExpand } = useSearchExpansion(
    directory?.container || [],
    searchTerm
  );

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

  return directory ? (
    <>
      <DirectoryBrowser
        directoryTree={directory.container}
        searchTerm={searchTerm}
        autoExpandedPaths={expandedPaths}
        isSearchActive={shouldAutoExpand}
      />
      <FileSystemFileMenu onClose={handleMenuClose} />
    </>
  ) : (
    <LinearProgress />
  );
};
