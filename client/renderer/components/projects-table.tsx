"use client";
import React, { useEffect, useMemo, useState, useCallback } from "react";
import { useRouter, useSearchParams } from "next/navigation";
import useSWR from "swr";
import {
  Alert,
  AlertTitle,
  Box,
  Button,
  Card,
  CardContent,
  Checkbox,
  Chip,
  IconButton,
  LinearProgress,
  Paper,
  Skeleton,
  Stack,
  Theme,
  Tooltip,
  Typography,
} from "@mui/material";
import {
  Clear,
  Delete,
  Download,
  Edit,
  DriveFileMove,
  FolderOpen,
  LabelOutlined,
  LinkOff,
  Schedule,
  Science,
  StarBorder,
  Science as CampaignIcon,
} from "@mui/icons-material";
import { alpha } from "@mui/material/styles";
import { useApi } from "../api";
import { apiFetch } from "../api-fetch";
import {
  Project,
  ProjectTagTree,
  TagFilter,
} from "../types/models";
import { shortDate } from "../pipes";
import { useDeleteDialog } from "../providers/delete-dialog";
import { useSet } from "../hooks";
import SearchField from "./search-field";
import { usePopcorn } from "../providers/popcorn-provider";
import { DataTable, Column } from "./data-table";
import { VirtualizedCardGrid } from "./virtualized-card-grid";
import { ProjectTagChips } from "./project-tag-chips";
import ProjectTagTreePane, { PROJECT_DRAG_TYPE } from "./project-tag-tree";
import TagSelectionDialog from "./tag-selection-dialog";
import { ViewMode, ViewModeToggle } from "./view-mode-toggle";
import { MoveProjectDialog } from "./move-project-dialog";
import { ReconnectProjectsDialog } from "./reconnect-projects-dialog";
import { isElectron } from "../utils/platform";

// Type for campaign info returned from API
interface CampaignInfo {
  campaign_id: number;
  campaign_name: string;
  member_count?: number;
  membership_type: "parent" | "member";
}

// Hook to fetch campaign info for projects (includes both parent and member campaigns)
// Uses POST to avoid URL length limits with large numbers of projects
function useProjectCampaigns(projectIds: number[]) {
  // Use a stable key for SWR based on sorted IDs
  const cacheKey = projectIds.length > 0
    ? `project_campaigns:${projectIds.slice().sort((a, b) => a - b).join(",")}`
    : null;

  const { data } = useSWR<Record<string, CampaignInfo>>(
    cacheKey,
    async () => {
      // Use POST to send project IDs in body (avoids URL length limits)
      const response = await apiFetch(
        "/api/proxy/ccp4i2/projectgroups/project_campaigns/",
        {
          method: "POST",
          headers: {
            "Content-Type": "application/json",
          },
          body: JSON.stringify({
            project_ids: projectIds,
            include_members: true,
          }),
        }
      );
      if (!response.ok) return {};
      return response.json();
    }
  );
  return data || {};
}

// Memoized ProjectCard component for performance
const ProjectCard = React.memo(
  ({
    project,
    isSelected,
    onToggleSelection,
    onNavigate,
    onEdit,
    onExport,
    onMove,
    onDelete,
  }: {
    project: Project;
    isSelected: boolean;
    onToggleSelection: () => void;
    onNavigate: () => void;
    onEdit: () => void;
    onExport: () => void;
    /** Omitted in the web build, where there is no local filesystem to move to. */
    onMove?: () => void;
    onDelete: () => void;
  }) => {
    return (
      <Card
        sx={isSelected ? sxSelectedCard : sxProjectCard}
        onClick={(event) => {
          if (!(event.target as HTMLElement).closest(".action-button")) {
            onNavigate();
          }
        }}
      >
        <CardContent sx={{ p: 2.5, "&:last-child": { pb: 2.5 } }}>
          <Box
            sx={{
              display: "flex",
              alignItems: "flex-start",
              justifyContent: "space-between",
              mb: 2,
            }}
          >
            <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
              <Checkbox
                className="action-button"
                size="small"
                checked={isSelected}
                onChange={onToggleSelection}
                sx={{ p: 0.5 }}
              />
            </Box>
            <IconButton
              className="action-button"
              size="small"
              sx={{ color: "text.disabled" }}
            >
              <StarBorder fontSize="small" />
            </IconButton>
          </Box>

          <Box sx={{ display: "flex", alignItems: "center", gap: 1.5, mb: 2 }}>
            <Box
              sx={{
                width: 48,
                height: 48,
                borderRadius: 2,
                bgcolor: "primary.50",
                display: "flex",
                alignItems: "center",
                justifyContent: "center",
                border: "1px solid",
                borderColor: "primary.200",
              }}
            >
              <Science sx={{ color: "primary.main", fontSize: 24 }} />
            </Box>
            <Box sx={{ minWidth: 0, flexGrow: 1 }}>
              <Typography
                variant="h6"
                sx={{
                  fontWeight: 600,
                  fontSize: "1.1rem",
                  lineHeight: 1.3,
                  overflow: "hidden",
                  textOverflow: "ellipsis",
                  whiteSpace: "nowrap",
                }}
              >
                {project.name}
              </Typography>
              {project.description && (
                <Typography
                  variant="body2"
                  color="text.secondary"
                  sx={{
                    mt: 0.25,
                    display: "-webkit-box",
                    WebkitLineClamp: 2,
                    WebkitBoxOrient: "vertical",
                    overflow: "hidden",
                  }}
                >
                  {project.description}
                </Typography>
              )}
            </Box>
          </Box>

          <Stack spacing={1} sx={{ mb: 2 }}>
            <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
              <Schedule sx={{ fontSize: 16, color: "text.disabled" }} />
              <Typography variant="caption" color="text.secondary">
                Created {shortDate(project.creation_time)}
              </Typography>
            </Box>
            <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
              <Schedule sx={{ fontSize: 16, color: "text.disabled" }} />
              <Typography variant="caption" color="text.secondary">
                Last accessed {shortDate(project.last_access)}
              </Typography>
            </Box>
          </Stack>

          {/* Project Tags */}
          <ProjectTagChips tags={project.tags} maxVisible={2} size="small" />

          <Stack
            direction="row"
            spacing={1}
            justifyContent="flex-end"
            sx={{ mt: 2 }}
          >
            <Tooltip title="Edit project name, description and tags">
              <IconButton
                className="action-button"
                size="small"
                onClick={(event) => {
                  event.stopPropagation();
                  onEdit();
                }}
                sx={{
                  color: "primary.main",
                  "&:hover": { bgcolor: "primary.50" },
                }}
              >
                <Edit fontSize="small" />
              </IconButton>
            </Tooltip>
            <Tooltip title="Export project">
              <IconButton
                className="action-button"
                size="small"
                onClick={(event) => {
                  event.stopPropagation();
                  onExport();
                }}
                sx={{
                  color: "primary.main",
                  "&:hover": { bgcolor: "primary.50" },
                }}
              >
                <Download fontSize="small" />
              </IconButton>
            </Tooltip>
            {onMove && (
              <Tooltip title="Move project on disk">
                <IconButton
                  className="action-button"
                  size="small"
                  onClick={(event) => {
                    event.stopPropagation();
                    onMove();
                  }}
                  sx={{
                    color: "primary.main",
                    "&:hover": { bgcolor: "primary.50" },
                  }}
                >
                  <DriveFileMove fontSize="small" />
                </IconButton>
              </Tooltip>
            )}
            <Tooltip title="Delete project">
              <IconButton
                className="action-button"
                size="small"
                onClick={(event) => {
                  event.stopPropagation();
                  onDelete();
                }}
                sx={{ color: "error.main", "&:hover": { bgcolor: "error.50" } }}
              >
                <Delete fontSize="small" />
              </IconButton>
            </Tooltip>
          </Stack>
        </CardContent>
      </Card>
    );
  }
);

ProjectCard.displayName = "ProjectCard";

const sxProjectCard = {
  cursor: "pointer",
  transition: "all 0.2s ease-in-out",
  border: "1px solid",
  borderColor: "divider",
  "&:hover": {
    transform: "translateY(-2px)",
    boxShadow: (theme: Theme) => theme.shadows[8],
    borderColor: "primary.main",
  },
};

const sxSelectedCard = {
  ...sxProjectCard,
  borderColor: "primary.main",
  bgcolor: (theme: Theme) => alpha(theme.palette.primary.main, 0.08),
};

export default function ProjectsTable() {
  const api = useApi();
  const router = useRouter();
  const { data: projects, mutate } = api.get<Project[]>("projects");
  const selectedIds = useSet<number>([]);
  const [query, setQuery] = useState("");
  // Which node of the tag tree the list is showing. The tree is fetched here
  // as well as in the pane purely for its mutate handle — SWR shares the key,
  // so it is still one request.
  const [tagFilter, setTagFilter] = useState<TagFilter>({ kind: "all" });
  const { data: tagTree, mutate: mutateTagTree } =
    api.get<ProjectTagTree>("projecttags/tree/");

  // ?tag=<id> arrives from a tag chip in the project menu bar. Resolving it
  // needs the tree, because a filter carries the node's path as well as its
  // id — so this waits for the tree rather than guessing.
  const searchParams = useSearchParams();
  const requestedTagId = searchParams?.get("tag");
  useEffect(() => {
    if (!requestedTagId || !tagTree) return;
    const wanted = Number(requestedTagId);
    const find = (nodes: typeof tagTree.tags): (typeof tagTree.tags)[0] | null => {
      for (const node of nodes) {
        if (node.id === wanted) return node;
        const inChildren = find(node.children);
        if (inChildren) return inChildren;
      }
      return null;
    };
    const node = find(tagTree.tags);
    if (node) {
      setTagFilter({
        kind: "tag",
        id: node.id,
        path: node.path,
        label: node.display_path,
      });
    }
  }, [requestedTagId, tagTree]);
  const [tagDialogOpen, setTagDialogOpen] = useState(false);
  const [viewMode, setViewMode] = useState<ViewMode>("list");
  const deleteDialog = useDeleteDialog();
  const { setMessage } = usePopcorn();
  const [projectToMove, setProjectToMove] = useState<Project | null>(null);
  // Moving a project means moving a directory on the machine running the
  // server, with a native folder picker to choose where. Neither exists in the
  // browser build, so the action is desktop-only. Resolved in an effect rather
  // than inline so the server-rendered markup and the first client render
  // agree.
  const [canMoveProjects, setCanMoveProjects] = useState(false);
  useEffect(() => setCanMoveProjects(isElectron()), []);

  const [reconnectOpen, setReconnectOpen] = useState(false);
  // Projects whose recorded directory is no longer on disk - a renamed drive,
  // a projects folder moved outside CCP4i2. Desktop only: on the web the
  // directories live on the server and the user cannot go and find them.
  const { data: missingResponse, mutate: refreshMissing } = api.get<any>(
    canMoveProjects ? "projects/missing_directories/" : null
  );
  const brokenProjects = useMemo(
    () => (missingResponse?.data ?? missingResponse)?.projects ?? [],
    [missingResponse]
  );
  const brokenIds = useMemo(
    () => new Set<number>(brokenProjects.map((p: any) => p.id)),
    [brokenProjects]
  );

  // Get campaign info for all projects
  const projectIds = (projects || []).map((p) => p.id);
  const campaignInfo = useProjectCampaigns(projectIds);

  // Does a project sit at or below the selected node of the tag tree?
  const matchesTagFilter = (project: Project) => {
    if (tagFilter.kind === "all") return true;
    const tags = Array.isArray(project.tags) ? project.tags : [];
    const paths = tags
      .map((tag) => (typeof tag === "object" ? tag.path : undefined))
      .filter((path): path is string => Boolean(path));
    if (tagFilter.kind === "untagged") return tags.length === 0;
    return paths.some(
      (path) =>
        path === tagFilter.path || path.startsWith(`${tagFilter.path}\x1f`)
    );
  };

  // Filter and sort projects
  const filteredProjects = useMemo(() => {
    if (!Array.isArray(projects)) return [];
    const term = query.toLowerCase();

    return projects
      .filter(matchesTagFilter)
      .filter((p) => {
        if (!term) return true;
        if (p.name.toLowerCase().includes(term)) return true;
        if (p.description?.toLowerCase().includes(term)) return true;
        if (Array.isArray(p.tags)) {
          return p.tags.some(
            (tag) => typeof tag === "object" && tag.text?.toLowerCase().includes(term)
          );
        }
        return false;
      })
      .sort(
        (a, b) =>
          new Date(b.last_access).getTime() - new Date(a.last_access).getTime()
      );
  }, [projects, query, tagFilter]);

  // Handlers
  async function afterTagging(label: string, count: number) {
    await Promise.all([mutate(), mutateTagTree()]);
    selectedIds.clear();
    setMessage(
      `Tagged ${count} project${count === 1 ? "" : "s"} with "${label}"`,
      "success"
    );
  }

  // Dropping the selection onto a node of the tree files projects under that tag.
  async function dropSelectionOnTag(tagId: number, tagLabel: string) {
    const ids = Array.from(selectedIds.values());
    if (ids.length === 0) return;
    try {
      await api.post(`projecttags/${tagId}/add_projects/`, {
        project_ids: ids,
      });
      await afterTagging(tagLabel, ids.length);
    } catch (error) {
      setMessage(
        `Failed to tag projects: ${error instanceof Error ? error.message : String(error)}`,
        "error"
      );
    }
  }

  function deleteSelected() {
    const selectedProjects = projects?.filter((project) =>
      selectedIds.has(project.id)
    );
    if (selectedProjects) deleteProjects(selectedProjects);
  }

  function deleteProjects(projectsToDelete: Project[]) {
    if (deleteDialog)
      deleteDialog({
        type: "show",
        what:
          projectsToDelete.length === 1
            ? projectsToDelete[0].name
            : `${projectsToDelete.length} projects`,
        onDelete: () => {
          const promises = projectsToDelete.map((project) => {
            selectedIds.delete(project.id);
            return api.delete(`projects/${project.id}`);
          });
          Promise.all(promises).then(() => mutate());
        },
      });
  }

  async function exportProject(project: Project) {
    try {
      const exportResult: any = await api.post(
        `projects/${project.id}/export/`,
        {}
      );

      if (exportResult?.success === false) {
        setMessage(
          `Failed to export "${project.name}": ${exportResult?.error || "Unknown error"}`,
          "error"
        );
        return;
      }

      setMessage(
        `Export started for "${project.name}". Available in File/Projects → Exports when complete.`,
        "success"
      );
    } catch (error) {
      console.error(`Failed to export project "${project.name}":`, error);
      setMessage(
        `Failed to start export for "${project.name}": ${error instanceof Error ? error.message : String(error)}`,
        "error"
      );
    }
  }

  async function exportSelected() {
    try {
      const selectedProjects = projects?.filter((project) =>
        selectedIds.has(project.id)
      );

      if (!selectedProjects || selectedProjects.length === 0) {
        setMessage("No projects selected for export.", "warning");
        return;
      }

      const exportPromises = selectedProjects.map(async (project) => {
        try {
          const result: any = await api.post(
            `projects/${project.id}/export/`,
            {}
          );
          if (result?.success === false) {
            return {
              project: project.name,
              success: false,
              error: result?.error,
            };
          }
          return { project: project.name, success: true, result };
        } catch (error) {
          console.error(`Failed to export project "${project.name}":`, error);
          return { project: project.name, success: false, error };
        }
      });

      const results = await Promise.all(exportPromises);
      const successful = results.filter((r) => r.success);
      const failed = results.filter((r) => !r.success);

      if (successful.length > 0 && failed.length === 0) {
        setMessage(
          `Started exports for ${successful.length} project${successful.length > 1 ? "s" : ""}. Available in File/Projects → Exports when complete.`,
          "success"
        );
      } else if (failed.length > 0 && successful.length === 0) {
        setMessage(
          `Failed to export ${failed.length} project${failed.length > 1 ? "s" : ""}: ${failed.map((r) => r.project).join(", ")}`,
          "error"
        );
      } else if (successful.length > 0 && failed.length > 0) {
        setMessage(
          `Exported ${successful.length}, failed ${failed.length}: ${failed.map((r) => r.project).join(", ")}`,
          "warning"
        );
      }
    } catch (error) {
      console.error("Failed to export selected projects:", error);
      setMessage(
        `Failed to start exports: ${error instanceof Error ? error.message : String(error)}`,
        "error"
      );
    }
  }

  function toggleAll() {
    if (filteredProjects) {
      if (selectedIds.size === filteredProjects.length) {
        selectedIds.clear();
      } else {
        filteredProjects.forEach((project) => selectedIds.add(project.id));
      }
    }
  }

  // Table columns definition
  const tableColumns: Column<Project>[] = [
    {
      key: "checkbox",
      label: "",
      width: 50,
      header: (
        <Tooltip
          title={
            selectedIds.size === filteredProjects.length
              ? "Deselect all projects"
              : "Select all projects"
          }
        >
          <Checkbox
            size="small"
            checked={
              selectedIds.size > 0 &&
              selectedIds.size === filteredProjects.length
            }
            indeterminate={
              selectedIds.size > 0 &&
              selectedIds.size < filteredProjects.length
            }
            onChange={toggleAll}
          />
        </Tooltip>
      ),
      render: (_, project) => (
        <Checkbox
          checked={selectedIds.has(project.id)}
          onChange={(event) => {
            event.stopPropagation();
            selectedIds.has(project.id)
              ? selectedIds.delete(project.id)
              : selectedIds.add(project.id);
          }}
          onClick={(e) => e.stopPropagation()}
        />
      ),
    },
    {
      key: "name",
      label: "Project Name",
      sortable: true,
      searchable: true,
      render: (_, project) => {
        const campaign = campaignInfo[String(project.id)];
        const isParent = campaign?.membership_type === "parent";
        const isMember = campaign?.membership_type === "member";
        return (
          <Box sx={{ display: "flex", alignItems: "center", gap: 1.5 }}>
            <Box sx={{ minWidth: 0 }}>
              <Stack direction="row" spacing={0.5} alignItems="center">
                <Typography variant="body2" sx={{ fontWeight: 500 }} noWrap>
                  {project.name}
                </Typography>
                {brokenIds.has(project.id) && (
                  <Tooltip title="This project's folder cannot be found on disk">
                    <LinkOff sx={{ fontSize: 16, color: "warning.main" }} />
                  </Tooltip>
                )}
              </Stack>
              {project.description && (
                <Tooltip title={project.description}>
                  <Typography
                    variant="caption"
                    color="text.secondary"
                    noWrap
                    sx={{ display: "block" }}
                  >
                    {project.description}
                  </Typography>
                </Tooltip>
              )}
              <Stack direction="row" spacing={0.5} flexWrap="wrap" useFlexGap>
                {isParent && (
                  <Tooltip title={`Campaign parent: ${campaign.campaign_name} (${campaign.member_count} datasets)`}>
                    <Chip
                      icon={<CampaignIcon sx={{ fontSize: "12px !important" }} />}
                      label={`${campaign.member_count} datasets`}
                      size="small"
                      color="secondary"
                      sx={{ height: 16, fontSize: "0.65rem", mt: 0.5, cursor: "pointer" }}
                      onClick={(e) => {
                        e.stopPropagation();
                        router.push(`/ccp4i2/campaigns/${campaign.campaign_id}`);
                      }}
                    />
                  </Tooltip>
                )}
                {isMember && (
                  <Tooltip title={`Member of campaign: ${campaign.campaign_name}`}>
                    <Chip
                      icon={<CampaignIcon sx={{ fontSize: "12px !important" }} />}
                      label={campaign.campaign_name}
                      size="small"
                      color="secondary"
                      variant="outlined"
                      sx={{ height: 16, fontSize: "0.65rem", mt: 0.5, cursor: "pointer" }}
                      onClick={(e) => {
                        e.stopPropagation();
                        router.push(`/ccp4i2/campaigns/${campaign.campaign_id}`);
                      }}
                    />
                  </Tooltip>
                )}
              </Stack>
            </Box>
          </Box>
        );
      },
    },
    {
      key: "tags",
      label: "Tags",
      render: (_, project) => (
        <ProjectTagChips tags={project.tags} maxVisible={3} size="small" />
      ),
    },
    {
      key: "creation_time",
      label: "Created",
      sortable: true,
      width: 120,
      render: (value) => (
        <Typography variant="body2" color="text.secondary">
          {shortDate(value)}
        </Typography>
      ),
    },
    {
      key: "last_access",
      label: "Last Accessed",
      sortable: true,
      width: 120,
      render: (value) => (
        <Typography variant="body2" color="text.secondary">
          {shortDate(value)}
        </Typography>
      ),
    },
    {
      key: "actions",
      label: "",
      width: canMoveProjects ? 180 : 140,
      render: (_, project) => (
        <Stack
          direction="row"
          spacing={0.5}
          justifyContent="flex-end"
          onClick={(e) => e.stopPropagation()}
        >
          <Tooltip title="Edit project name, description and tags">
            <IconButton
              size="small"
              onClick={() => router.push(`/ccp4i2/edit-project/${project.id}`)}
              sx={{ color: "primary.main" }}
            >
              <Edit fontSize="small" />
            </IconButton>
          </Tooltip>
          <Tooltip title="Export project">
            <IconButton
              size="small"
              onClick={() => exportProject(project)}
              sx={{ color: "primary.main" }}
            >
              <Download fontSize="small" />
            </IconButton>
          </Tooltip>
          {canMoveProjects && !brokenIds.has(project.id) && (
            <Tooltip title="Move project on disk">
              <IconButton
                size="small"
                onClick={() => setProjectToMove(project)}
                sx={{ color: "primary.main" }}
              >
                <DriveFileMove fontSize="small" />
              </IconButton>
            </Tooltip>
          )}
          <Tooltip title="Delete project">
            <IconButton
              size="small"
              onClick={() => deleteProjects([project])}
              sx={{ color: "error.main" }}
            >
              <Delete fontSize="small" />
            </IconButton>
          </Tooltip>
        </Stack>
      ),
    },
  ];

  // Card renderer for virtualized grid
  const renderProjectCard = useCallback(
    (project: Project) => (
      <ProjectCard
        project={project}
        isSelected={selectedIds.has(project.id)}
        onToggleSelection={() => {
          selectedIds.has(project.id)
            ? selectedIds.delete(project.id)
            : selectedIds.add(project.id);
        }}
        onNavigate={() => router.push(`/ccp4i2/project/${project.id}`)}
        onEdit={() => router.push(`/ccp4i2/edit-project/${project.id}`)}
        onExport={() => exportProject(project)}
        onMove={
          canMoveProjects && !brokenIds.has(project.id)
            ? () => setProjectToMove(project)
            : undefined
        }
        onDelete={() => deleteProjects([project])}
      />
    ),
    [selectedIds, router, canMoveProjects, brokenIds]
  );

  if (projects === undefined) return <LinearProgress />;
  if (projects.length === 0)
    return (
      <Box sx={{ textAlign: "center", py: 8 }}>
        <FolderOpen sx={{ fontSize: 80, color: "text.disabled", mb: 2 }} />
        <Typography variant="h6" color="text.secondary" gutterBottom>
          No Projects Yet
        </Typography>
        <Typography variant="body2" color="text.disabled">
          Create your first crystallography project to get started
        </Typography>
      </Box>
    );

  return Array.isArray(projects) ? (
    <Box
      sx={{
        height: "100%",
        display: "flex",
        flexDirection: "column",
        overflow: "hidden",
      }}
    >
      {brokenProjects.length > 0 && (
        <Alert
          severity="warning"
          icon={<LinkOff />}
          sx={{ mb: 2, flexShrink: 0 }}
          action={
            <Button
              color="inherit"
              size="small"
              onClick={() => setReconnectOpen(true)}
            >
              Reconnect
            </Button>
          }
        >
          <AlertTitle>
            {brokenProjects.length} project
            {brokenProjects.length === 1 ? "" : "s"} cannot be found on disk
          </AlertTitle>
          If a drive was renamed or the projects folder moved, say where they
          are now and they will be reconnected. Nothing is moved or deleted.
        </Alert>
      )}

      {/* Header with search, view toggle, and selection actions */}
      <Box sx={{ mb: 2, flexShrink: 0 }}>
        {selectedIds.size === 0 ? (
          <Box
            sx={{
              display: "flex",
              alignItems: "center",
              justifyContent: "space-between",
              mb: 2,
            }}
          >
            <Stack direction="row" spacing={2} alignItems="center">
              {viewMode === "cards" && (
                <Tooltip
                  title={
                    selectedIds.size === filteredProjects.length
                      ? "Deselect all projects"
                      : "Select all projects"
                  }
                >
                  <Checkbox
                    size="small"
                    checked={
                      selectedIds.size > 0 &&
                      selectedIds.size === filteredProjects.length
                    }
                    indeterminate={
                      selectedIds.size > 0 &&
                      selectedIds.size < filteredProjects.length
                    }
                    onChange={toggleAll}
                  />
                </Tooltip>
              )}
              <Typography
                variant="h5"
                component="h2"
                sx={{ fontWeight: 600, color: "text.primary" }}
              >
                Your Projects
              </Typography>
            </Stack>
            <Stack direction="row" spacing={2} alignItems="center">
              <SearchField
                value={query}
                onChange={setQuery}
                placeholder="Search projects..."
                size="small"
              />
              <ViewModeToggle mode={viewMode} onChange={setViewMode} />
            </Stack>
          </Box>
        ) : (
          <Paper
            elevation={2}
            sx={{
              p: 2,
              mb: 2,
              bgcolor: "primary.50",
              border: "1px solid",
              borderColor: "primary.200",
              borderRadius: 2,
            }}
          >
            <Stack direction="row" alignItems="center" spacing={2}>
              <Tooltip title="Clear selection">
                <IconButton onClick={selectedIds.clear} size="small">
                  <Clear />
                </IconButton>
              </Tooltip>
              <Tooltip title="Drag onto a tag to file these projects there">
                <Typography
                  color="primary.main"
                  variant="subtitle1"
                  draggable
                  onDragStart={(event) => {
                    event.dataTransfer.setData(
                      PROJECT_DRAG_TYPE,
                      JSON.stringify(Array.from(selectedIds.values()))
                    );
                    event.dataTransfer.effectAllowed = "copy";
                  }}
                  sx={{ fontWeight: 600, flexGrow: 1, cursor: "grab" }}
                >
                  {selectedIds.size} project{selectedIds.size !== 1 ? "s" : ""}{" "}
                  selected
                </Typography>
              </Tooltip>
              <Tooltip title="Add a tag to selected projects">
                <IconButton
                  onClick={() => setTagDialogOpen(true)}
                  size="small"
                  color="primary"
                >
                  <LabelOutlined />
                </IconButton>
              </Tooltip>
              <Tooltip title="Export selected projects">
                <IconButton
                  onClick={exportSelected}
                  size="small"
                  color="primary"
                >
                  <Download />
                </IconButton>
              </Tooltip>
              <Tooltip title="Delete selected projects">
                <IconButton onClick={deleteSelected} size="small" color="error">
                  <Delete />
                </IconButton>
              </Tooltip>
            </Stack>
          </Paper>
        )}
      </Box>

      {/* Tag tree beside the list, Gmail-style: the tree filters, it does not
          move anything. */}
      <Box
        sx={{
          flex: 1,
          minHeight: 0,
          display: "flex",
          gap: 2,
          overflow: "hidden",
        }}
      >
        <Box
          sx={{
            width: 220,
            flexShrink: 0,
            display: { xs: "none", md: "block" },
            borderRight: "1px solid",
            borderColor: "divider",
            minHeight: 0,
          }}
        >
          <ProjectTagTreePane
            selected={tagFilter}
            onSelect={setTagFilter}
            onDropProjects={dropSelectionOnTag}
            onTagsChanged={() => mutate()}
            totalCount={projects.length}
          />
        </Box>

      <Box sx={{ flex: 1, overflow: "hidden", minHeight: 0 }}>
        {viewMode === "cards" ? (
          // Virtualized Card Grid View
          <VirtualizedCardGrid
            data={filteredProjects}
            renderCard={renderProjectCard}
            getItemKey={(project) => project.id}
            columns={{ xs: 1, sm: 2, md: 3, lg: 4, xl: 4 }}
            estimateCardHeight={280}
            gap={24}
            fillHeight
            emptyMessage="No projects found"
          />
        ) : (
          // Virtualized Table View
          <DataTable
            data={filteredProjects}
            columns={tableColumns}
            getRowKey={(project) => project.id}
            onRowClick={(project) => router.push(`/ccp4i2/project/${project.id}`)}
            emptyMessage="No projects found"
          />
        )}

        {/* Show message when no projects match search */}
        {filteredProjects.length === 0 && query && (
          <Box sx={{ textAlign: "center", py: 4 }}>
            <Typography variant="h6" color="text.secondary" gutterBottom>
              No projects found
            </Typography>
            <Typography variant="body2" color="text.disabled">
              Try adjusting your search term
            </Typography>
          </Box>
        )}

        {canMoveProjects && (
          <>
            <MoveProjectDialog
              project={projectToMove}
              open={Boolean(projectToMove)}
              onClose={() => setProjectToMove(null)}
              onMoved={() => {
                mutate();
                refreshMissing();
              }}
            />
            <ReconnectProjectsDialog
              open={reconnectOpen}
              onClose={() => setReconnectOpen(false)}
              brokenProjects={brokenProjects}
              onReconnected={() => {
                mutate();
                refreshMissing();
              }}
            />
          </>
        )}
      </Box>
      </Box>

      <TagSelectionDialog
        open={tagDialogOpen}
        projectIds={Array.from(selectedIds.values())}
        onClose={() => setTagDialogOpen(false)}
        onApplied={afterTagging}
      />
    </Box>
  ) : (
    <Skeleton variant="rectangular" height={400} />
  );
}
