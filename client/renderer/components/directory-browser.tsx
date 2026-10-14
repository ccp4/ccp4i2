import React, { useCallback, useEffect, useMemo, useRef, useState } from "react";
import useSWR from "swr";
import { useVirtualizer } from "@tanstack/react-virtual";
import {
  Box,
  CircularProgress,
  IconButton,
  LinearProgress,
  Stack,
  Tooltip,
  Typography,
} from "@mui/material";
import {
  ExpandMore as ExpandMoreIcon,
  ChevronRight as ChevronRightIcon,
  Folder as FolderIcon,
  FolderOpen as FolderOpenIcon,
  InsertDriveFile as FileIcon,
  MoreVert as MoreVertIcon,
} from "@mui/icons-material";
import { apiJson } from "../api-fetch";
import { useFileSystemFileBrowser } from "../providers/file-system-file-browser-context";
import { DirectoryEntry, directoryRows } from "../lib/directory-rows";
import { formatModified, formatSize, formatTimestamp } from "../lib/file-kinds";
import SearchField from "./search-field";

export interface FileSystemItem {
  path: string;
  name: string;
  type: "directory" | "file";
  size?: number;
  mtime?: number;
  link?: boolean;
}

export interface DirectoryBrowserProps {
  /**
   * The listing endpoint of what is browsed, e.g. ``jobs/12/directory`` or
   * ``projects/3/directory``. It answers ``?path=`` with one folder and
   * ``?search=`` with matching names from the whole tree.
   */
  endpoint: string | null;
  /** Re-list the open folders this often (ms); 0 never. */
  pollInterval?: number;
  /** Lists again whenever this changes (e.g. a job's status). */
  revision?: string | number;
  fileFilter?: (item: FileSystemItem) => boolean;
  onItemClick?: (item: FileSystemItem, event: React.MouseEvent) => void;
}

interface SearchResult {
  items: (FileSystemItem & { relative: string })[];
  truncated: boolean;
}

const ROW_HEIGHT = 30;
const INDENT = 16;

const listingUrl = (endpoint: string, query: Record<string, string>) =>
  `${endpoint}?${new URLSearchParams(query).toString()}`;

/** Every open folder, root included, listed in one SWR entry. */
async function fetchListings(
  [endpoint, , ...folders]: string[]
): Promise<Record<string, DirectoryEntry[]>> {
  const listed = await Promise.all(
    folders.map(async (folder) => {
      const body = await apiJson<{ container?: DirectoryEntry[] }>(
        listingUrl(endpoint, { path: folder })
      );
      return [folder, Array.isArray(body?.container) ? body.container : []] as const;
    })
  );
  return Object.fromEntries(listed);
}

function useDebounced<T>(value: T, delay: number): T {
  const [debounced, setDebounced] = useState(value);
  useEffect(() => {
    const timer = setTimeout(() => setDebounced(value), delay);
    return () => clearTimeout(timer);
  }, [value, delay]);
  return debounced;
}

/**
 * A file browser that lists a folder only when it is opened.
 *
 * It used to receive the whole tree, walked and stat'ed on the server and
 * re-fetched every few seconds while a job ran: for a PanDDA run over
 * hundreds of datasets that was megabytes per poll from a network share,
 * and every node of it mounted in the page. Now the open folders are listed
 * (and re-listed when polling), only the rows on screen are drawn, and
 * search asks the server, which can see folders the browser has not opened.
 */
const DirectoryBrowser: React.FC<DirectoryBrowserProps> = ({
  endpoint,
  pollInterval = 0,
  revision,
  fileFilter,
  onItemClick,
}) => {
  const [open, setOpen] = useState<Set<string>>(new Set());
  const [searchTerm, setSearchTerm] = useState("");
  const query = useDebounced(searchTerm.trim(), 300);
  const { openMenu } = useFileSystemFileBrowser();
  const scrollRef = useRef<HTMLDivElement>(null);

  // "4 min ago" must not stay four minutes old on a job that has stopped polling.
  const [now, setNow] = useState(() => Date.now());
  useEffect(() => {
    const timer = setInterval(() => setNow(Date.now()), 60000);
    return () => clearInterval(timer);
  }, []);

  // A different root starts closed.
  useEffect(() => setOpen(new Set()), [endpoint]);

  const folders = useMemo(() => ["", ...Array.from(open).sort()], [open]);
  const { data: listings, isLoading } = useSWR(
    endpoint ? [endpoint, String(revision ?? ""), ...folders] : null,
    fetchListings,
    { refreshInterval: pollInterval, revalidateOnFocus: false, keepPreviousData: true }
  );

  const { data: found, isLoading: searching } = useSWR<SearchResult>(
    endpoint && query ? listingUrl(endpoint, { search: query }) : null,
    async (url: string) => (await apiJson<{ container: SearchResult }>(url)).container,
    { revalidateOnFocus: false, keepPreviousData: true }
  );

  const rows = useMemo(
    () => (listings ? directoryRows(listings, open, fileFilter) : []),
    [listings, open, fileFilter]
  );

  const searchRows = useMemo(
    () =>
      (found?.items ?? []).filter(
        (item) => item.type === "directory" || !fileFilter || fileFilter(item)
      ),
    [found, fileFilter]
  );

  const inSearch = Boolean(query);
  const count = inSearch ? searchRows.length : rows.length;

  const virtualizer = useVirtualizer({
    count,
    getScrollElement: () => scrollRef.current,
    estimateSize: () => ROW_HEIGHT,
    overscan: 20,
  });

  const toggle = useCallback((relative: string) => {
    setOpen((prev) => {
      const next = new Set(prev);
      if (next.has(relative)) {
        // Closing a folder closes what was open inside it.
        for (const path of prev) {
          if (path === relative || path.startsWith(`${relative}/`)) next.delete(path);
        }
      } else {
        next.add(relative);
      }
      return next;
    });
  }, []);

  /** Show a search hit in the tree: open its folders and leave search. */
  const reveal = useCallback((relative: string, isDirectory: boolean) => {
    const parts = relative.split("/");
    const upTo = isDirectory ? parts.length : parts.length - 1;
    setOpen((prev) => {
      const next = new Set(prev);
      for (let i = 1; i <= upTo; i++) next.add(parts.slice(0, i).join("/"));
      return next;
    });
    setSearchTerm("");
  }, []);

  const handleMenuOpen = useCallback(
    (item: FileSystemItem, element: HTMLElement) => {
      // A fixed anchor that outlives the row: virtualised rows unmount on scroll.
      const rect = element.getBoundingClientRect();
      document.getElementById("file-menu-anchor")?.remove();
      const anchor = document.createElement("div");
      Object.assign(anchor.style, {
        position: "fixed",
        top: `${rect.bottom}px`,
        left: `${rect.left}px`,
        width: "1px",
        height: "1px",
        pointerEvents: "none",
        visibility: "hidden",
        zIndex: "9999",
      });
      anchor.id = "file-menu-anchor";
      document.body.appendChild(anchor);
      openMenu(anchor, item);
    },
    [openMenu]
  );

  const renderEntry = (
    item: FileSystemItem,
    depth: number,
    label: string,
    caret: "open" | "closed" | "none",
    onActivate: (event: React.MouseEvent) => void
  ) => (
    <Box
      data-tree-item-row={item.path}
      sx={{
        display: "flex",
        alignItems: "center",
        height: ROW_HEIGHT,
        pl: `${depth * INDENT}px`,
        cursor: "context-menu",
        borderRadius: 1,
        "&:hover": { backgroundColor: "action.hover", "& .menu-button": { opacity: 1 } },
      }}
      onClick={(event) => {
        event.stopPropagation();
        onActivate(event);
        onItemClick?.(item, event);
      }}
      onContextMenu={(event) => {
        event.preventDefault();
        event.stopPropagation();
        handleMenuOpen(item, event.currentTarget as HTMLElement);
      }}
    >
      <Box sx={{ width: 28, display: "flex", justifyContent: "center", flexShrink: 0 }}>
        {caret === "open" && <ExpandMoreIcon fontSize="small" />}
        {caret === "closed" && <ChevronRightIcon fontSize="small" />}
      </Box>
      {item.type === "directory" ? (
        caret === "open" ? (
          <FolderOpenIcon fontSize="small" sx={{ mr: 1, color: "primary.main" }} />
        ) : (
          <FolderIcon fontSize="small" sx={{ mr: 1, color: "primary.main" }} />
        )
      ) : (
        <FileIcon fontSize="small" sx={{ mr: 1, color: "text.secondary" }} />
      )}
      <Typography
        variant="body2"
        title={label}
        sx={{
          fontFamily: item.type === "directory" ? "inherit" : "monospace",
          fontSize: "0.875rem",
          flex: 1,
          minWidth: 0,
          overflow: "hidden",
          textOverflow: "ellipsis",
          whiteSpace: "nowrap",
          fontStyle: item.link ? "italic" : "normal",
        }}
      >
        {label}
      </Typography>
      {item.type === "file" && (
        <Typography
          variant="caption"
          sx={{ color: "text.secondary", mr: 1, flexShrink: 0, minWidth: 64, textAlign: "right" }}
        >
          {formatSize(item.size)}
        </Typography>
      )}
      {item.mtime !== undefined && (
        <Tooltip title={`Modified ${formatTimestamp(item.mtime)}`} placement="left">
          <Typography
            variant="caption"
            sx={{ color: "text.secondary", mr: 1, flexShrink: 0, minWidth: 72, textAlign: "right" }}
          >
            {formatModified(item.mtime, now)}
          </Typography>
        </Tooltip>
      )}
      <IconButton
        className="menu-button"
        size="small"
        sx={{ opacity: 0, transition: "opacity 0.2s", p: 0.25, flexShrink: 0 }}
        onClick={(event) => {
          event.preventDefault();
          event.stopPropagation();
          handleMenuOpen(item, event.currentTarget as HTMLElement);
        }}
      >
        <MoreVertIcon fontSize="small" />
      </IconButton>
    </Box>
  );

  const renderRow = (index: number) => {
    if (inSearch) {
      const hit = searchRows[index];
      return renderEntry(hit, 0, hit.relative, "none", () =>
        hit.type === "directory" ? reveal(hit.relative, true) : undefined
      );
    }
    const row = rows[index];
    if (row.kind === "loading") {
      return (
        <Box sx={{ display: "flex", alignItems: "center", height: ROW_HEIGHT, pl: `${row.depth * INDENT + 28}px` }}>
          <CircularProgress size={14} sx={{ mr: 1 }} />
          <Typography variant="caption" color="text.secondary">Listing…</Typography>
        </Box>
      );
    }
    const isDirectory = row.entry.type === "directory";
    return renderEntry(
      row.entry,
      row.depth,
      row.entry.name,
      isDirectory ? (row.open ? "open" : "closed") : "none",
      () => isDirectory && toggle(row.relative)
    );
  };

  let message: string | null = null;
  if (inSearch && !searching && count === 0) message = "No matching files found";
  if (!inSearch && listings && count === 0) message = "No files found";

  return (
    <Stack sx={{ flex: 1, overflow: "hidden" }}>
      <SearchField
        value={searchTerm}
        onChange={setSearchTerm}
        placeholder="Search files..."
        size="small"
        sx={{ px: 2, py: 1 }}
      />
      {(isLoading || (inSearch && searching)) && <LinearProgress />}
      {inSearch && found?.truncated && (
        <Typography variant="caption" color="text.secondary" sx={{ px: 2 }}>
          Showing the first {found.items.length} matches; narrow the search to see the rest.
        </Typography>
      )}
      <Box ref={scrollRef} sx={{ flex: 1, overflow: "auto", scrollbarWidth: "thin", p: 1, pb: 4 }}>
        {message ? (
          <Typography variant="body2" color="text.secondary" sx={{ p: 2, textAlign: "center" }}>
            {message}
          </Typography>
        ) : (
          <Box sx={{ height: virtualizer.getTotalSize(), position: "relative" }}>
            {virtualizer.getVirtualItems().map((virtualRow) => (
              <Box
                key={virtualRow.key}
                sx={{
                  position: "absolute",
                  top: 0,
                  left: 0,
                  right: 0,
                  transform: `translateY(${virtualRow.start}px)`,
                }}
              >
                {renderRow(virtualRow.index)}
              </Box>
            ))}
          </Box>
        )}
      </Box>
    </Stack>
  );
};

export default DirectoryBrowser;
