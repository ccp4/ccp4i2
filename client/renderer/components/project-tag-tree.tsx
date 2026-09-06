"use client";
import React, { useCallback, useEffect, useRef, useState } from "react";
import {
  Box,
  Collapse,
  IconButton,
  List,
  ListItemButton,
  ListItemText,
  Menu,
  MenuItem,
  Skeleton,
  TextField,
  Tooltip,
  Typography,
} from "@mui/material";
import {
  Add,
  ChevronRight,
  DriveFileRenameOutline,
  Delete,
  ExpandMore,
  Inbox,
  LabelOutlined,
  Layers,
  MoreVert,
} from "@mui/icons-material";
import { alpha, Theme, useTheme } from "@mui/material/styles";

import { useApi } from "../api";
import { useDeleteDialog } from "../providers/delete-dialog";
import { usePopcorn } from "../providers/popcorn-provider";
import { ProjectTagNode, ProjectTagTree, TagFilter } from "../types/models";

/**
 * The tag forest, as a browsable and editable tree beside the project list.
 *
 * Tags are the mechanism for organising projects — see
 * docs/organising-projects.md for why hierarchy lives here rather than on
 * ProjectGroup. A node answers for everything filed anywhere beneath it, so an
 * ancestor with no projects of its own is still worth clicking.
 *
 * The tree is also where the hierarchy is *built*: add a child to a node,
 * rename it, drag it onto another node to re-parent it. Nesting is expressed
 * by the gesture rather than by typing separators into a tag's name, which is
 * what lets a tag be called anything at all — including "A/B".
 */

export const PROJECT_DRAG_TYPE = "application/x-ccp4i2-projects";
export const TAG_DRAG_TYPE = "application/x-ccp4i2-tag";

/** Matches ProjectTag.PATH_SEPARATOR server-side. */
const PATH_SEPARATOR = "\x1f";

interface DraggedTag {
  id: number;
  path: string;
}

/**
 * What a tag drag carries in its `TAG_DRAG_TYPE` payload. The label is
 * included so a drop target elsewhere (a project row) can name the tag in
 * its confirmation without having the tree to hand.
 */
export interface TagDragPayload {
  id: number;
  label: string;
}

/** Decode a `TAG_DRAG_TYPE` payload, or null if the drag is not a tag. */
export function readTagDragPayload(
  dataTransfer: DataTransfer
): TagDragPayload | null {
  if (!Array.from(dataTransfer.types).includes(TAG_DRAG_TYPE)) return null;
  try {
    const parsed = JSON.parse(dataTransfer.getData(TAG_DRAG_TYPE));
    if (typeof parsed?.id !== "number") return null;
    return { id: parsed.id, label: String(parsed.label ?? "") };
  } catch {
    return null;
  }
}

/**
 * How a hovered drop target announces itself — shared by the tree's nodes and
 * by the project rows and cards a tag can be dropped onto, so the affordance
 * reads the same wherever a drag happens to be.
 */
export const dropTargetSx = {
  bgcolor: (theme: Theme) => alpha(theme.palette.primary.main, 0.16),
  outline: "1px dashed",
  outlineColor: "primary.main",
};

interface ProjectTagTreeProps {
  selected: TagFilter;
  onSelect: (filter: TagFilter) => void;
  /** Called when a project selection is dropped onto a tag node. */
  onDropProjects?: (tagId: number, tagLabel: string) => void;
  /** Called after the tree itself changes, so the project list can refresh. */
  onTagsChanged?: () => void;
  /** Total projects, for the "All projects" node. */
  totalCount: number;
}

function isSameFilter(a: TagFilter, b: TagFilter) {
  if (a.kind !== b.kind) return false;
  if (a.kind === "tag" && b.kind === "tag") return a.id === b.id;
  return true;
}

function countDescendants(node: ProjectTagNode): number {
  return node.children.reduce(
    (total, child) => total + 1 + countDescendants(child),
    0
  );
}

export default function ProjectTagTreePane({
  selected,
  onSelect,
  onDropProjects,
  onTagsChanged,
  totalCount,
}: ProjectTagTreeProps) {
  const api = useApi();
  const { data: tree, mutate: mutateTree } =
    api.get<ProjectTagTree>("projecttags/tree/");
  const { setMessage } = usePopcorn();
  const deleteDialog = useDeleteDialog();

  const [collapsed, setCollapsed] = useState<Set<number>>(new Set());
  const [dropTarget, setDropTarget] = useState<number | "root" | null>(null);
  const [menuFor, setMenuFor] = useState<{
    anchor: HTMLElement;
    node: ProjectTagNode;
  } | null>(null);
  // A menu action is run only once the menu has finished closing. Both "Add
  // tag inside" and "Rename" open a self-focusing field, and MUI's focus trap
  // stays active for the whole exit transition: it pulls focus back out of the
  // new field, which blurs it, which cancelled it a frame after it appeared.
  // Deferring sidesteps the trap rather than arguing with it.
  const [pendingAction, setPendingAction] = useState<{
    kind: "add" | "rename" | "delete";
    node: ProjectTagNode;
  } | null>(null);
  const [renaming, setRenaming] = useState<{ id: number; text: string } | null>(
    null
  );
  const [adding, setAdding] = useState<{
    parentId: number | null;
    text: string;
  } | null>(null);
  // Held in a ref because dragover can only read dataTransfer *types*, not its
  // payload — and refusing a drop onto the dragged node's own subtree has to
  // happen while hovering, not after the drop.
  const draggedTag = useRef<DraggedTag | null>(null);
  // The element the browser paints under the cursor during a tag drag. Made
  // on dragstart and removed on dragend; without it the default ghost is a
  // snapshot of the whole row, which on a wide list reads as "the list".
  const dragGhost = useRef<HTMLElement | null>(null);
  const theme = useTheme();
  const confirmDeleteRef = useRef<((node: ProjectTagNode) => void) | null>(null);

  const toggle = useCallback((id: number) => {
    setCollapsed((previous) => {
      const next = new Set(previous);
      if (next.has(id)) next.delete(id);
      else next.add(id);
      return next;
    });
  }, []);

  const expand = useCallback((id: number) => {
    setCollapsed((previous) => {
      if (!previous.has(id)) return previous;
      const next = new Set(previous);
      next.delete(id);
      return next;
    });
  }, []);

  const runPendingAction = useCallback(() => {
    if (!pendingAction) return;
    const { kind, node } = pendingAction;
    setPendingAction(null);
    if (kind === "add") {
      expand(node.id);
      setAdding({ parentId: node.id, text: "" });
    } else if (kind === "rename") {
      setRenaming({ id: node.id, text: node.text });
    } else {
      confirmDeleteRef.current?.(node);
    }
  }, [expand, pendingAction]);

  const refresh = useCallback(async () => {
    await mutateTree();
    onTagsChanged?.();
  }, [mutateTree, onTagsChanged]);

  /** A 404 means the tag is already gone, which is what the caller wanted. */
  const isAlreadyGone = (error: unknown) =>
    error instanceof Error && error.message.includes("404");

  const reportFailure = useCallback(
    (what: string, error: unknown) => {
      const detail = error instanceof Error ? error.message : String(error);
      setMessage(`Could not ${what}: ${detail}`, "error");
    },
    [setMessage]
  );

  // -- mutations -------------------------------------------------------

  const createTag = useCallback(
    async (parentId: number | null, text: string) => {
      const trimmed = text.trim();
      if (!trimmed) return;
      try {
        await api.post("projecttags/", {
          text: trimmed,
          parent: parentId,
          projects: [],
        });
        if (parentId !== null) expand(parentId);
      } catch (error) {
        reportFailure(`create the tag "${trimmed}"`, error);
      } finally {
        await refresh();
      }
    },
    [api, expand, refresh, reportFailure]
  );

  const renameTag = useCallback(
    async (node: ProjectTagNode, text: string) => {
      const trimmed = text.trim();
      if (!trimmed || trimmed === node.text) return;
      try {
        await api.patch(`projecttags/${node.id}/`, { text: trimmed });
      } catch (error) {
        reportFailure(`rename "${node.text}"`, error);
      } finally {
        await refresh();
      }
    },
    [api, refresh, reportFailure]
  );

  const reparentTag = useCallback(
    async (dragged: DraggedTag, newParent: ProjectTagNode | null) => {
      try {
        await api.patch(`projecttags/${dragged.id}/`, {
          parent: newParent ? newParent.id : null,
        });
        if (newParent) expand(newParent.id);
      } catch (error) {
        reportFailure("move that tag", error);
      } finally {
        await refresh();
      }
    },
    [api, expand, refresh, reportFailure]
  );

  const confirmDelete = useCallback(
    (node: ProjectTagNode) => {
      const descendants = countDescendants(node);
      const notes: React.ReactNode[] = [];
      if (descendants > 0) {
        notes.push(
          <Typography key="subtree" variant="body2" sx={{ mt: 1 }}>
            This also deletes the {descendants} tag
            {descendants === 1 ? "" : "s"} nested below it.
          </Typography>
        );
      }
      notes.push(
        <Typography
          key="projects"
          variant="body2"
          color="text.secondary"
          sx={{ mt: 1 }}
        >
          No projects are deleted — they simply lose{" "}
          {descendants > 0 ? "these labels" : "this label"}.
        </Typography>
      );

      deleteDialog?.({
        type: "show",
        what: `the tag "${node.display_path}"`,
        children: notes,
        onDelete: async () => {
          try {
            await api.delete(`projecttags/${node.id}/`);
          } catch (error) {
            // Deleting something already deleted is not a failure to report:
            // the tree can lag the database (a second click, another window),
            // and the user's intent is satisfied either way.
            if (!isAlreadyGone(error)) {
              reportFailure(`delete "${node.text}"`, error);
            }
          } finally {
            // Always reconcile with the server. Refreshing only on success
            // leaves a phantom node on screen precisely when the tree and the
            // database have diverged.
            if (selected.kind === "tag" && selected.id === node.id) {
              onSelect({ kind: "all" });
            }
            await refresh();
          }
        },
      });
    },
    [api, deleteDialog, onSelect, refresh, reportFailure, selected]
  );

  useEffect(() => {
    confirmDeleteRef.current = confirmDelete;
  }, [confirmDelete]);

  // -- drag and drop ---------------------------------------------------

  /**
   * May `dragged` be filed under `node` (null meaning "make it a root")?
   *
   * The dragged descriptor is a parameter rather than read from the ref here:
   * the drop handler has to clear that ref, and an earlier version cleared it
   * *before* asking this question, so the answer was always "no" and a drop
   * that had been highlighted as acceptable then did nothing.
   */
  const canAcceptTagDrop = useCallback(
    (dragged: DraggedTag | null, node: ProjectTagNode | null) => {
      if (!dragged) return false;
      if (node === null) return dragged.path.includes(PATH_SEPARATOR); // already a root
      if (node.id === dragged.id) return false;
      // Refuse the node's own subtree: the server rejects it too, but a
      // highlighted target that then fails is a worse answer than no highlight.
      return !node.path.startsWith(dragged.path + PATH_SEPARATOR);
    },
    []
  );

  const handleDragOver = useCallback(
    (event: React.DragEvent, node: ProjectTagNode | null) => {
      const types = event.dataTransfer.types;
      if (types.includes(TAG_DRAG_TYPE)) {
        if (!canAcceptTagDrop(draggedTag.current, node)) return;
        // A tag drag allows copy (onto a project) as well as move (here);
        // say which this is so the cursor matches the outcome.
        event.dataTransfer.dropEffect = "move";
      } else if (types.includes(PROJECT_DRAG_TYPE)) {
        if (!onDropProjects || node === null) return;
      } else {
        return;
      }
      event.preventDefault();
      setDropTarget(node ? node.id : "root");
    },
    [canAcceptTagDrop, onDropProjects]
  );

  const handleDrop = useCallback(
    (event: React.DragEvent, node: ProjectTagNode | null) => {
      const types = event.dataTransfer.types;
      setDropTarget(null);
      if (types.includes(TAG_DRAG_TYPE)) {
        const dragged = draggedTag.current;
        draggedTag.current = null;
        if (!dragged || !canAcceptTagDrop(dragged, node)) return;
        event.preventDefault();
        reparentTag(dragged, node);
        return;
      }
      if (types.includes(PROJECT_DRAG_TYPE) && onDropProjects && node) {
        event.preventDefault();
        onDropProjects(node.id, node.display_path);
      }
    },
    [canAcceptTagDrop, onDropProjects, reparentTag]
  );

  const startTagDrag = useCallback(
    (event: React.DragEvent, node: ProjectTagNode) => {
      draggedTag.current = { id: node.id, path: node.path };
      const payload: TagDragPayload = { id: node.id, label: node.display_path };
      event.dataTransfer.setData(TAG_DRAG_TYPE, JSON.stringify(payload));
      // Move within the tree, copy onto a project — the drop target picks.
      event.dataTransfer.effectAllowed = "copyMove";

      // The ghost must be a rendered element, so it goes into the document
      // (off screen) rather than being built detached.
      if (typeof event.dataTransfer.setDragImage === "function") {
        const ghost = document.createElement("span");
        ghost.textContent = node.text;
        Object.assign(ghost.style, {
          position: "fixed",
          top: "-1000px",
          left: "-1000px",
          padding: "2px 10px",
          borderRadius: "4px",
          whiteSpace: "nowrap",
          font: `${theme.typography.body2.fontSize} ${theme.typography.fontFamily}`,
          color: theme.palette.text.primary,
          background: theme.palette.background.paper,
          border: `1px solid ${theme.palette.primary.main}`,
          boxShadow: theme.shadows[2],
        } as Partial<CSSStyleDeclaration>);
        document.body.appendChild(ghost);
        dragGhost.current = ghost;
        event.dataTransfer.setDragImage(ghost, 0, 0);
      }
    },
    [theme]
  );

  const endTagDrag = useCallback(() => {
    draggedTag.current = null;
    dragGhost.current?.remove();
    dragGhost.current = null;
    setDropTarget(null);
  }, []);

  // -- rendering -------------------------------------------------------

  const addRow = useCallback(
    (parentId: number | null, depth: number) => (
      <Box sx={{ pl: 1 + depth * 1.5, pr: 1, py: 0.5 }}>
        <TextField
          autoFocus
          fullWidth
          size="small"
          variant="standard"
          placeholder="New tag name"
          value={adding?.text ?? ""}
          onChange={(event) =>
            setAdding({ parentId, text: event.target.value })
          }
          onBlur={() => {
            // Commit rather than discard: clicking away from a half-typed tag
            // should not silently bin it. Escape is the way to abandon it.
            const text = adding?.text ?? "";
            setAdding(null);
            if (text.trim()) createTag(parentId, text);
          }}
          onKeyDown={(event) => {
            if (event.key === "Enter") {
              const text = adding?.text ?? "";
              setAdding(null);
              createTag(parentId, text);
            } else if (event.key === "Escape") {
              setAdding(null);
            }
          }}
        />
      </Box>
    ),
    [adding, createTag]
  );

  const renderNode = useCallback(
    (node: ProjectTagNode): React.ReactNode => {
      const expanded = !collapsed.has(node.id);
      const isSelected = isSameFilter(selected, {
        kind: "tag",
        id: node.id,
        path: node.path,
        label: node.text,
      });
      const isDropTarget = dropTarget === node.id;
      const isRenaming = renaming?.id === node.id;

      return (
        <React.Fragment key={node.id}>
          <ListItemButton
            dense
            selected={isSelected}
            draggable={!isRenaming}
            onDragStart={(event) => startTagDrag(event, node)}
            onDragEnd={endTagDrag}
            onClick={() =>
              onSelect({
                kind: "tag",
                id: node.id,
                path: node.path,
                label: node.display_path,
              })
            }
            onDragOver={(event) => handleDragOver(event, node)}
            onDragLeave={() =>
              setDropTarget((current) => (current === node.id ? null : current))
            }
            onDrop={(event) => handleDrop(event, node)}
            sx={{
              pl: 1 + node.depth * 1.5,
              borderRadius: 1,
              "&:hover .tag-node-actions": { opacity: 1 },
              ...(isDropTarget && dropTargetSx),
            }}
          >
            {node.children.length > 0 ? (
              <IconButton
                size="small"
                sx={{ mr: 0.5, p: 0.25 }}
                onClick={(event) => {
                  event.stopPropagation();
                  toggle(node.id);
                }}
              >
                {expanded ? (
                  <ExpandMore fontSize="small" />
                ) : (
                  <ChevronRight fontSize="small" />
                )}
              </IconButton>
            ) : (
              <LabelOutlined
                sx={{ mr: 0.5, fontSize: 16, color: "text.disabled" }}
              />
            )}

            {isRenaming ? (
              <TextField
                autoFocus
                fullWidth
                size="small"
                variant="standard"
                value={renaming.text}
                onClick={(event) => event.stopPropagation()}
                onChange={(event) =>
                  setRenaming({ id: node.id, text: event.target.value })
                }
                onBlur={() => {
                  const text = renaming.text;
                  setRenaming(null);
                  if (text.trim() && text.trim() !== node.text) {
                    renameTag(node, text);
                  }
                }}
                onKeyDown={(event) => {
                  event.stopPropagation();
                  if (event.key === "Enter") {
                    const text = renaming.text;
                    setRenaming(null);
                    renameTag(node, text);
                  } else if (event.key === "Escape") {
                    setRenaming(null);
                  }
                }}
              />
            ) : (
              <>
                <ListItemText
                  primary={node.text}
                  primaryTypographyProps={{
                    variant: "body2",
                    noWrap: true,
                    fontWeight: isSelected ? 600 : 400,
                  }}
                />
                <Box
                  className="tag-node-actions"
                  sx={{ opacity: 0, transition: "opacity 120ms" }}
                >
                  <Tooltip title="Add a tag inside this one">
                    <IconButton
                      size="small"
                      aria-label={`Add a tag inside ${node.text}`}
                      sx={{ p: 0.25 }}
                      onClick={(event) => {
                        event.stopPropagation();
                        expand(node.id);
                        setAdding({ parentId: node.id, text: "" });
                      }}
                    >
                      <Add fontSize="small" />
                    </IconButton>
                  </Tooltip>
                  <IconButton
                    size="small"
                    aria-label={`More actions for ${node.text}`}
                    sx={{ p: 0.25 }}
                    onClick={(event) => {
                      event.stopPropagation();
                      setMenuFor({ anchor: event.currentTarget, node });
                    }}
                  >
                    <MoreVert fontSize="small" />
                  </IconButton>
                </Box>
                <Tooltip
                  title={
                    node.total_project_count === node.project_count
                      ? `${node.project_count} project${node.project_count === 1 ? "" : "s"}`
                      : `${node.project_count} here, ${node.total_project_count} including everything below`
                  }
                >
                  <Typography
                    variant="caption"
                    color="text.secondary"
                    sx={{ ml: 1 }}
                  >
                    {node.total_project_count}
                  </Typography>
                </Tooltip>
              </>
            )}
          </ListItemButton>

          {(node.children.length > 0 || adding?.parentId === node.id) && (
            <Collapse in={expanded} timeout="auto" unmountOnExit>
              <List disablePadding>
                {node.children.map(renderNode)}
                {adding?.parentId === node.id &&
                  addRow(node.id, node.depth + 1)}
              </List>
            </Collapse>
          )}
        </React.Fragment>
      );
    },
    [
      addRow,
      adding,
      collapsed,
      dropTarget,
      expand,
      handleDragOver,
      handleDrop,
      onSelect,
      renameTag,
      renaming,
      selected,
      toggle,
    ]
  );

  if (!tree) {
    return (
      <Skeleton variant="rectangular" height={200} sx={{ borderRadius: 1 }} />
    );
  }

  const hasTags = tree.tags.length > 0;

  return (
    <Box sx={{ height: "100%", overflow: "auto", pr: 1 }}>
      <Box
        sx={{
          display: "flex",
          alignItems: "center",
          justifyContent: "space-between",
          pl: 1,
          pr: 0.5,
          pt: 0.5,
        }}
      >
        <Typography variant="overline" color="text.secondary">
          Tags
        </Typography>
        <Tooltip title="New top-level tag">
          <IconButton
            size="small"
            aria-label="New top-level tag"
            sx={{ p: 0.25 }}
            onClick={() => setAdding({ parentId: null, text: "" })}
          >
            <Add fontSize="small" />
          </IconButton>
        </Tooltip>
      </Box>

      <List disablePadding>
        <ListItemButton
          dense
          selected={selected.kind === "all"}
          onClick={() => onSelect({ kind: "all" })}
          onDragOver={(event) => handleDragOver(event, null)}
          onDragLeave={() =>
            setDropTarget((current) => (current === "root" ? null : current))
          }
          onDrop={(event) => handleDrop(event, null)}
          sx={{
            borderRadius: 1,
            ...(dropTarget === "root" && dropTargetSx),
          }}
        >
          <Layers sx={{ mr: 0.5, fontSize: 16, color: "text.disabled" }} />
          <ListItemText
            primary="All projects"
            primaryTypographyProps={{
              variant: "body2",
              fontWeight: selected.kind === "all" ? 600 : 400,
            }}
          />
          <Typography variant="caption" color="text.secondary">
            {totalCount}
          </Typography>
        </ListItemButton>

        {tree.tags.map(renderNode)}
        {adding?.parentId === null && addRow(null, 0)}

        {tree.untagged_project_count > 0 && (
          <ListItemButton
            dense
            selected={selected.kind === "untagged"}
            onClick={() => onSelect({ kind: "untagged" })}
            sx={{ borderRadius: 1 }}
          >
            <Inbox sx={{ mr: 0.5, fontSize: 16, color: "text.disabled" }} />
            <ListItemText
              primary="Untagged"
              primaryTypographyProps={{
                variant: "body2",
                fontWeight: selected.kind === "untagged" ? 600 : 400,
              }}
            />
            <Typography variant="caption" color="text.secondary">
              {tree.untagged_project_count}
            </Typography>
          </ListItemButton>
        )}
      </List>

      <Menu
        open={Boolean(menuFor)}
        anchorEl={menuFor?.anchor ?? null}
        onClose={() => setMenuFor(null)}
        disableRestoreFocus
        TransitionProps={{ onExited: runPendingAction }}
      >
        <MenuItem
          onClick={() => {
            if (menuFor) setPendingAction({ kind: "add", node: menuFor.node });
            setMenuFor(null);
          }}
        >
          <Add fontSize="small" sx={{ mr: 1 }} />
          Add tag inside
        </MenuItem>
        <MenuItem
          onClick={() => {
            if (menuFor) setPendingAction({ kind: "rename", node: menuFor.node });
            setMenuFor(null);
          }}
        >
          <DriveFileRenameOutline fontSize="small" sx={{ mr: 1 }} />
          Rename
        </MenuItem>
        <MenuItem
          onClick={() => {
            if (menuFor) setPendingAction({ kind: "delete", node: menuFor.node });
            setMenuFor(null);
          }}
        >
          <Delete fontSize="small" sx={{ mr: 1 }} />
          Delete
        </MenuItem>
      </Menu>

      {!hasTags && !adding && (
        <Typography
          variant="caption"
          color="text.disabled"
          sx={{ display: "block", px: 1, pt: 1 }}
        >
          No tags yet. Use + above to make one, then drag tags onto each other
          to nest them.
        </Typography>
      )}
    </Box>
  );
}
