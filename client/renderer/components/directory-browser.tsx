/*
 * Copyright (C) 2025 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
import React, {
  useState,
  useMemo,
  useCallback,
  useEffect,
} from "react";
import {
  Box,
  Typography,
  IconButton,
  Collapse,
  Stack,
} from "@mui/material";
import {
  ExpandMore as ExpandMoreIcon,
  ChevronRight as ChevronRightIcon,
  Folder as FolderIcon,
  FolderOpen as FolderOpenIcon,
  InsertDriveFile as FileIcon,
  MoreVert as MoreVertIcon,
} from "@mui/icons-material";
import { useFileSystemFileBrowser } from "../providers/file-system-file-browser-context";
import SearchField from "./search-field";

export interface FileSystemItem {
  path: string;
  name: string;
  type: "directory" | "file";
  size?: number;
  mode?: number;
  inode?: number;
  device?: number;
  nlink?: number;
  uid?: number;
  gid?: number;
  atime?: number;
  mtime?: number;
  ctime?: number;
  contents?: FileSystemItem[];
}

export interface DirectoryBrowserProps {
  directoryTree: any[];
  fileFilter?: (item: FileSystemItem) => boolean;
  onItemClick?: (item: FileSystemItem, event: React.MouseEvent) => void;
  searchTerm?: string;
  autoExpandedPaths?: Set<string>;
  isSearchActive?: boolean;
}

interface TreeNodeProps {
  item: FileSystemItem;
  level: number;
  expandedNodes: Set<string>;
  onToggleExpand: (path: string) => void;
  searchTerm: string;
  fileFilter?: (item: FileSystemItem) => boolean;
  onItemClick?: (item: FileSystemItem, event: React.MouseEvent) => void;
}

const DirectoryBrowser: React.FC<DirectoryBrowserProps> = ({
  directoryTree,
  fileFilter,
  onItemClick,
  searchTerm,
  autoExpandedPaths = new Set(),
  isSearchActive = false,
}) => {
  const [expandedNodes, setExpandedNodes] = useState<Set<string>>(new Set());
  const [searchTermState, setSearchTerm] = useState<string>("");

  // Use the FileSystemFileBrowser context
  const { openMenu, anchorEl, menuNode } = useFileSystemFileBrowser();

  const sortedDirectoryTree = useMemo(() => {
    if (!directoryTree) return [];

    // Helper to extract trailing integer from name pattern "name_{integer}"
    const extractTrailingInt = (name: string) => {
      const match = name.match(/_(\d+)$/);
      return match ? parseInt(match[1], 10) : null;
    };

    // Recursive sorting function
    const sortItems = (items: FileSystemItem[]): FileSystemItem[] => {
      return items
        .slice()
        .sort((a: FileSystemItem, b: FileSystemItem) => {
          // Directories first
          if (a.type === "directory" && b.type !== "directory") return -1;
          if (a.type !== "directory" && b.type === "directory") return 1;

          // Check for name_{int} pattern
          const aInt = extractTrailingInt(a.name);
          const bInt = extractTrailingInt(b.name);

          if (aInt !== null && bInt !== null) {
            // If both have trailing ints, sort numerically
            if (a.name.replace(/_\d+$/, "") === b.name.replace(/_\d+$/, "")) {
              return aInt - bInt;
            }
          }

          // Fallback to lexicographical
          return a.name.localeCompare(b.name);
        })
        .map((item) => {
          // Recursively sort directory contents
          if (item.type === "directory" && item.contents) {
            return {
              ...item,
              contents: sortItems(item.contents),
            };
          }
          return item;
        });
    };

    return sortItems(directoryTree);
  }, [directoryTree]);

  // Filter directory tree based on provided filter function
  const filteredTree = useMemo(() => {
    if (!fileFilter) return sortedDirectoryTree;

    const filterItems = (items: FileSystemItem[]): FileSystemItem[] => {
      return items
        .map((item) => {
          if (item.type === "directory") {
            const filteredContents = item.contents
              ? filterItems(item.contents)
              : [];
            if (filteredContents.length > 0) {
              return { ...item, contents: filteredContents };
            }
            return null;
          } else if (item.type === "file" && fileFilter(item)) {
            return item;
          }
          return null;
        })
        .filter((item): item is FileSystemItem => item !== null);
    };

    return filterItems(sortedDirectoryTree);
  }, [sortedDirectoryTree, fileFilter]);

  // Search functionality
  const searchFilteredTree = useMemo(() => {
    if (!searchTermState.trim()) return filteredTree;

    const filterBySearch = (items: FileSystemItem[]): FileSystemItem[] => {
      return items
        .map((item) => {
          if (item.type === "directory") {
            const filteredContents = item.contents
              ? filterBySearch(item.contents)
              : [];
            if (
              filteredContents.length > 0 ||
              item.name.toLowerCase().includes(searchTermState.toLowerCase())
            ) {
              return { ...item, contents: filteredContents };
            }
            return null;
          } else if (
            item.name.toLowerCase().includes(searchTermState.toLowerCase())
          ) {
            return item;
          }
          return null;
        })
        .filter((item): item is FileSystemItem => item !== null);
    };

    return filterBySearch(filteredTree);
  }, [filteredTree, searchTermState]);

  const handleToggleExpand = useCallback((path: string) => {
    setExpandedNodes((prev) => {
      const newSet = new Set(prev);
      if (newSet.has(path)) {
        newSet.delete(path);
      } else {
        newSet.add(path);
      }
      return newSet;
    });
  }, []);

  // Handle menu opening using the context
  const handleMenuOpen = useCallback(
    (item: FileSystemItem, element: HTMLElement) => {
      // Capture the position immediately while the element is still valid
      const rect = element.getBoundingClientRect();

      // Create a stable virtual anchor to avoid DOM removal issues
      const virtualAnchor = document.createElement("div");
      virtualAnchor.style.position = "fixed";
      virtualAnchor.style.top = `${rect.bottom}px`;
      virtualAnchor.style.left = `${rect.left}px`;
      virtualAnchor.style.width = "1px";
      virtualAnchor.style.height = "1px";
      virtualAnchor.style.pointerEvents = "none";
      virtualAnchor.style.visibility = "hidden";
      virtualAnchor.style.zIndex = "9999";
      virtualAnchor.id = "file-menu-anchor";

      // Add to DOM immediately
      document.body.appendChild(virtualAnchor);

      // Use the context's openMenu function
      openMenu(virtualAnchor, item);
    },
    [openMenu]
  );

  const TreeNode: React.FC<TreeNodeProps> = ({
    item,
    level,
    expandedNodes,
    onToggleExpand,
    searchTerm,
    fileFilter,
    onItemClick,
  }) => {
    const isExpanded = expandedNodes.has(item.path);
    const hasChildren =
      item.type === "directory" && item.contents && item.contents.length > 0;

    // Highlight search term
    const highlightText = (text: string, term: string) => {
      if (!term.trim()) return text;
      const regex = new RegExp(`(${term})`, "gi");
      const parts = text.split(regex);
      return parts.map((part, index) =>
        regex.test(part) ? (
          <mark key={index} style={{ backgroundColor: "#ffeb3b", padding: 0 }}>
            {part}
          </mark>
        ) : (
          part
        )
      );
    };

    const handleClick = (event: React.MouseEvent) => {
      event.stopPropagation();
      if (item.type === "directory") {
        onToggleExpand(item.path);
      }
      onItemClick?.(item, event);
    };

    const handleDoubleClick = (event: React.MouseEvent) => {
      event.stopPropagation();
    };

    const handleRightClick = (event: React.MouseEvent) => {
      event.preventDefault();
      event.stopPropagation();
      handleMenuOpen(item, event.currentTarget as HTMLElement);
    };

    const handleMenuClick = (event: React.MouseEvent) => {
      event.preventDefault();
      event.stopPropagation();
      handleMenuOpen(item, event.currentTarget as HTMLElement);
    };

    return (
      <Box>
        <Box
          data-tree-item-row={item.path}
          sx={{
            display: "flex",
            alignItems: "center",
            paddingLeft: `${level * 16}px`,
            paddingY: 0.5,
            cursor: "context-menu", // <-- Change cursor to indicate contextual menu
            backgroundColor: "transparent",
            "&:hover": {
              backgroundColor: "action.hover",
              "& .menu-button": {
                opacity: 1,
              },
            },
            borderRadius: 1,
            margin: 0.25,
            position: "relative",
          }}
          onClick={handleClick}
          onDoubleClick={handleDoubleClick}
          onContextMenu={handleRightClick}
        >
          {/* Always show caret area for directories to maintain alignment */}
          {item.type === "directory" ? (
            <IconButton
              size="small"
              sx={{
                padding: 0.25,
                marginRight: 0.5,
                visibility: hasChildren ? "visible" : "hidden",
              }}
              onClick={(e) => {
                e.stopPropagation();
                if (hasChildren) {
                  onToggleExpand(item.path);
                }
              }}
            >
              {hasChildren ? (
                isExpanded ? (
                  <ExpandMoreIcon fontSize="small" />
                ) : (
                  <ChevronRightIcon fontSize="small" />
                )
              ) : (
                <ChevronRightIcon fontSize="small" style={{ opacity: 0 }} />
              )}
            </IconButton>
          ) : (
            <Box sx={{ width: 32, height: 24, marginRight: 0.5 }} />
          )}

          {item.type === "directory" ? (
            isExpanded ? (
              <FolderOpenIcon
                fontSize="small"
                sx={{ marginRight: 1, color: "primary.main" }}
              />
            ) : (
              <FolderIcon
                fontSize="small"
                sx={{ marginRight: 1, color: "primary.main" }}
              />
            )
          ) : (
            <FileIcon
              fontSize="small"
              sx={{ marginRight: 1, color: "text.secondary" }}
            />
          )}

          <Typography
            variant="body2"
            sx={{
              fontFamily: item.type === "directory" ? "inherit" : "monospace",
              fontSize: "0.875rem",
              color: "text.primary",
              fontWeight: "normal",
              flex: 1,
              minWidth: 0,
              overflow: "hidden",
              textOverflow: "ellipsis",
              whiteSpace: "nowrap",
            }}
          >
            {highlightText(item.name, searchTerm)}
          </Typography>

          {item.type === "file" && item.size && (
            <Typography
              variant="caption"
              sx={{
                color: "text.secondary",
                fontSize: "0.75rem",
                marginRight: 1,
                flexShrink: 0,
              }}
            >
              {(item.size / 1024).toFixed(1)} KB
            </Typography>
          )}

          <IconButton
            className="menu-button"
            size="small"
            sx={{
              opacity: 0,
              transition: "opacity 0.2s",
              padding: 0.25,
              marginLeft: 0.5,
              flexShrink: 0,
            }}
            onClick={handleMenuClick}
          >
            <MoreVertIcon fontSize="small" />
          </IconButton>
        </Box>

        {item.type === "directory" && hasChildren && (
          <Collapse in={isExpanded}>
            <Box>
              {item.contents!.map((child) => (
                <TreeNode
                  key={child.path}
                  item={child}
                  level={level + 1}
                  expandedNodes={expandedNodes}
                  onToggleExpand={onToggleExpand}
                  searchTerm={searchTerm ?? ""}
                  fileFilter={fileFilter}
                  onItemClick={onItemClick}
                />
              ))}
            </Box>
          </Collapse>
        )}
      </Box>
    );
  };

  // Local state for manually expanded items (user clicked)
  const [manuallyExpanded, setManuallyExpanded] = useState<Set<string>>(
    new Set()
  );

  // Function to check if an item should be expanded
  const isExpanded = useCallback(
    (itemPath: string) => {
      if (isSearchActive) {
        // During search, use auto-expanded paths
        return autoExpandedPaths.has(itemPath);
      } else {
        // When not searching, use manually expanded paths
        return manuallyExpanded.has(itemPath);
      }
    },
    [isSearchActive, autoExpandedPaths, manuallyExpanded]
  );

  // Handle manual toggle (user clicking expand/collapse)
  const toggleExpansion = useCallback(
    (itemPath: string) => {
      if (isSearchActive) {
        // Don't allow manual toggle during search
        return;
      }

      setManuallyExpanded((prev) => {
        const newSet = new Set(prev);
        if (newSet.has(itemPath)) {
          newSet.delete(itemPath);
        } else {
          newSet.add(itemPath);
        }
        return newSet;
      });
    },
    [isSearchActive]
  );

  // Clear manual expansions when search becomes active
  useEffect(() => {
    if (isSearchActive) {
      setManuallyExpanded(new Set());
    }
  }, [isSearchActive]);

  // Recursive function to render directory items
  const renderDirectoryItem = useCallback(
    (item: any, currentPath: string = "") => {
      const itemPath = currentPath ? `${currentPath}/${item.name}` : item.name;
      const expanded = isExpanded(itemPath);
      const hasChildren = item.children && item.children.length > 0;

      // Check if this item or any children match the search
      const matchesSearch = searchTerm
        ? item.name.toLowerCase().includes(searchTerm.toLowerCase())
        : false;

      return (
        <div key={itemPath}>
          <div
            onClick={() => hasChildren && toggleExpansion(itemPath)}
            style={{
              cursor: hasChildren ? "pointer" : "default",
              padding: "4px 8px",
              backgroundColor: matchesSearch ? "#ffeb3b" : "transparent", // Highlight matches
              fontWeight: matchesSearch ? "bold" : "normal",
            }}
          >
            {hasChildren && (
              <span style={{ marginRight: "8px" }}>{expanded ? "▼" : "▶"}</span>
            )}
            {item.name}
          </div>

          {hasChildren && expanded && (
            <div style={{ marginLeft: "20px" }}>
              {item.children.map((child: any) =>
                renderDirectoryItem(child, itemPath)
              )}
            </div>
          )}
        </div>
      );
    },
    [isExpanded, toggleExpansion, searchTerm]
  );

  return (
    <Stack sx={{ flex: 1, "overflow": "hidden" }}>
      <SearchField
        value={searchTermState}
        onChange={setSearchTerm}
        placeholder="Search files..."
        size="small"
        sx={{ px: 2, py: 1 }}
      />
      <Box
        sx={{
          flex: 1, // This will take remaining space
          overflow: "auto", // This creates the scrollable area
          scrollbarWidth: "thin",
          p: 1,
          pb: 4,
        }}
      >
        {searchFilteredTree.length > 0 ? (
          searchFilteredTree.map((item) => (
            <TreeNode
              key={item.path}
              item={item}
              level={0}
              expandedNodes={expandedNodes}
              onToggleExpand={handleToggleExpand}
              searchTerm={searchTerm || ""}
              fileFilter={fileFilter}
              onItemClick={onItemClick}
            />
          ))
        ) : (
          <Typography
            variant="body2"
            color="text.secondary"
            sx={{ p: 2, textAlign: "center" }}
          >
            {searchTerm ? "No matching files found" : "No files found"}
          </Typography>
        )}
      </Box>
    </Stack>
  );
};

export default DirectoryBrowser;
