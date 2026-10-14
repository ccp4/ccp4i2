"use client";
import { SyntheticEvent, useEffect, useMemo, useState } from "react";
import {
  Box,
  Divider,
  ListItemIcon,
  ListItemText,
  Menu,
  MenuItem,
  Typography,
} from "@mui/material";
import {
  Download,
  Image as ImageIcon,
  OpenInNew,
  Preview,
  TableChart,
  ViewInAr,
} from "@mui/icons-material";
import useSWR from "swr";
import { useCCP4i2Window } from "../app-context";
import { doDownload, useApi } from "../api";
import { apiBlob, apiJson, apiPost } from "../api-fetch";
import { Project } from "../types/models";
import { useFilePreviewContext } from "../providers/file-preview-context";
import { useFileSystemFileBrowser } from "../providers/file-system-file-browser-context";
import { FileSystemItem } from "./directory-browser";
import {
  fileKind,
  formatSize,
  formatTimestamp,
  MoorhenKind,
} from "../lib/file-kinds";

interface FileSystemFileMenuProps {
  onClose?: () => void;
}

/** What a folder holds that Moorhen can open, counted from its listing. */
function countLoadable(entries: FileSystemItem[]) {
  const counts: Record<MoorhenKind, number> = { coordinates: 0, map: 0, cif: 0 };
  for (const entry of entries) {
    if (entry.type !== "file") continue;
    const kind = fileKind(entry.name).moorhen;
    if (kind) counts[kind] += 1;
  }
  return counts;
}

const plural = (n: number, word: string) => `${n} ${word}${n === 1 ? "" : "s"}`;

export const FileSystemFileMenu: React.FC<FileSystemFileMenuProps> = ({
  onClose,
}) => {
  const { anchorEl, menuNode, closeMenu } = useFileSystemFileBrowser();
  const { projectId } = useCCP4i2Window();
  const { get, post } = useApi();
  const { data: project } = get<Project>(`projects/${projectId}`);
  const { setContentSpecification } = useFilePreviewContext();
  const [error, setError] = useState<string | null>(null);

  useEffect(() => setError(null), [menuNode]);

  const isOpen = Boolean(anchorEl && menuNode);
  const isDirectory = menuNode?.type === "directory";
  const kind = useMemo(() => (menuNode ? fileKind(menuNode.name) : null), [menuNode]);

  /** The node's path relative to the project directory. */
  const relative = useMemo(() => {
    if (!menuNode || !project) return null;
    return menuNode.path.slice(project.directory.length + 1);
  }, [menuNode, project]);

  const fileUrl = relative !== null && project
    ? `/api/proxy/ccp4i2/projects/${project.id}/project_file?path=${encodeURIComponent(relative)}`
    : null;

  // A folder: list it (one level, cheap) to say what Moorhen would open.
  const { data: folderListing } = useSWR(
    isOpen && isDirectory && project && relative !== null
      ? `projects/${project.id}/directory?path=${encodeURIComponent(relative)}`
      : null,
    async (url: string) => (await apiJson<{ container?: FileSystemItem[] }>(url))?.container ?? [],
    { revalidateOnFocus: false }
  );
  const folderCounts = folderListing ? countLoadable(folderListing) : null;

  const handleClose = () => {
    closeMenu();
    onClose?.();
  };

  const act = (action: () => void | Promise<void>) => async (ev: SyntheticEvent) => {
    ev.stopPropagation();
    try {
      await action();
      handleClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : String(err));
    }
  };

  const handleDownload = act(() => {
    if (fileUrl && menuNode) doDownload(fileUrl, menuNode.name);
  });

  const handlePreview = act(() => {
    if (!fileUrl || !menuNode || !kind?.preview) return;
    setContentSpecification({
      url: fileUrl,
      title: menuNode.name,
      language: kind.preview,
    });
  });

  const handleMoorhen = act(() => {
    if (!project || relative === null) return;
    window.open(
      `/ccp4i2/moorhen-page/project/${project.id}?path=${encodeURIComponent(relative)}`,
      "_blank",
      "noopener,noreferrer"
    );
  });

  const handleNewTab = act(async () => {
    if (!project || relative === null || !kind?.newTab) return;
    if (kind.newTab === "pdf") {
      const blob = await apiBlob(fileUrl as string);
      window.open(URL.createObjectURL(new Blob([blob], { type: "application/pdf" })), "_blank");
      return;
    }
    // A new tab sends no Authorization header, nor do the images and links
    // the page then loads: seed it with a read grant for the page's folder,
    // as report file links do. The path-based URL keeps relative links working.
    const pageUrl = `/api/proxy/ccp4i2/projects/${project.id}/files_by_path/${relative}`;
    const folder = relative.includes("/") ? relative.replace(/\/[^/]*$/, "") : "";
    try {
      const response = await apiPost<{ data?: { grant?: string } }>(
        `projects/${project.id}/file_grant`,
        { path: folder }
      );
      const grant = response?.data?.grant;
      if (grant) {
        window.open(`${pageUrl}?file_grant=${encodeURIComponent(grant)}`, "_blank");
        return;
      }
    } catch (err) {
      console.error("Failed to mint file grant:", err);
    }
    window.open(pageUrl, "_blank");
  });

  const serverViewer = (viewer: string) =>
    act(async () => {
      if (!project || !menuNode) return;
      await post<any>(`projects/${project.id}/preview_file/`, {
        viewer,
        path: menuNode.path.slice(project.directory.length),
      });
    });

  if (!menuNode || !anchorEl) return null;

  const items: React.ReactNode[] = [];
  const item = (key: string, icon: React.ReactNode, label: string, onClick: (ev: SyntheticEvent) => void, secondary?: string, disabled = false) =>
    items.push(
      <MenuItem key={key} onClick={onClick} disabled={disabled} dense>
        <ListItemIcon>{icon}</ListItemIcon>
        <ListItemText primary={label} secondary={secondary} />
      </MenuItem>
    );

  if (isDirectory) {
    const loadable = folderCounts
      ? folderCounts.coordinates + folderCounts.map + folderCounts.cif
      : 0;
    const what = folderCounts
      ? [
          folderCounts.coordinates && plural(folderCounts.coordinates, "model"),
          folderCounts.map && plural(folderCounts.map, "map"),
          folderCounts.cif && plural(folderCounts.cif, "CIF file"),
        ].filter(Boolean).join(", ")
      : "";
    item(
      "moorhen",
      <ViewInAr fontSize="small" />,
      "Open folder in Moorhen",
      handleMoorhen,
      !folderCounts ? "Looking…" : loadable ? what : "No models or maps directly in this folder",
      !folderCounts || loadable === 0
    );
  } else if (kind) {
    if (kind.preview) {
      const icon = kind.preview === "image" ? <ImageIcon fontSize="small" />
        : kind.preview === "csv" ? <TableChart fontSize="small" />
        : <Preview fontSize="small" />;
      item("preview", icon, kind.previewLabel, handlePreview);
    }
    if (kind.newTab) {
      item("newtab", <OpenInNew fontSize="small" />, "Open in new tab", handleNewTab);
    }
    if (kind.moorhen) {
      item(
        "moorhen",
        <ViewInAr fontSize="small" />,
        "Open in Moorhen",
        handleMoorhen,
        kind.moorhen === "cif" ? "As a model, or as the monomers a dictionary describes" : undefined
      );
    }
    if (kind.coot) item("coot", <OpenInNew fontSize="small" />, "Coot", serverViewer("coot"));
    if (kind.ccp4mg) item("ccp4mg", <OpenInNew fontSize="small" />, "CCP4MG", serverViewer("ccp4mg"));
    if (kind.viewhkl) item("viewhkl", <OpenInNew fontSize="small" />, "ViewHKL", serverViewer("viewhkl"));
    item("download", <Download fontSize="small" />, "Download", handleDownload);
  }

  const details = [
    !isDirectory && formatSize(menuNode.size),
    formatTimestamp(menuNode.mtime) && `modified ${formatTimestamp(menuNode.mtime)}`,
  ].filter(Boolean).join(" · ");

  return (
    <Menu
      anchorEl={anchorEl}
      open={isOpen}
      onClose={handleClose}
      anchorOrigin={{ vertical: "top", horizontal: "left" }}
      transformOrigin={{ vertical: "top", horizontal: "left" }}
      disablePortal={true}
      keepMounted={false}
      disableScrollLock={true}
      slotProps={{
        paper: {
          sx: {
            zIndex: 10000,
            minWidth: 220,
            maxWidth: "min(420px, calc(100vw - 20px))",
            maxHeight: "calc(100vh - 20px)",
          },
        },
      }}
    >
      <Box sx={{ px: 2, pt: 0.5, pb: 1 }}>
        <Typography variant="body2" sx={{ fontWeight: 500, wordBreak: "break-all" }}>
          {menuNode.name}
        </Typography>
        {details && (
          <Typography variant="caption" color="text.secondary">
            {details}
          </Typography>
        )}
        {menuNode.link && (
          <Typography variant="caption" color="text.secondary" display="block">
            A link to a file elsewhere
          </Typography>
        )}
      </Box>
      <Divider />
      {items}
      {error && (
        <Typography variant="caption" color="error" sx={{ px: 2, display: "block" }}>
          {error}
        </Typography>
      )}
    </Menu>
  );
};
