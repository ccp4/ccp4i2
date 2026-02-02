"use client";
import React, { useMemo, useState, useCallback } from "react";
import { useRouter } from "next/navigation";
import useSWR from "swr";
import {
  Box,
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
  ToggleButton,
  ToggleButtonGroup,
} from "@mui/material";
import {
  Clear,
  Delete,
  Download,
  FolderOpen,
  Schedule,
  Science,
  StarBorder,
  ViewModule,
  ViewList,
  Science as CampaignIcon,
} from "@mui/icons-material";
import { alpha } from "@mui/material/styles";
import { useApi } from "../api";
import { apiFetch } from "../api-fetch";
import { Project } from "../types/models";
import { shortDate } from "../pipes";
import { useDeleteDialog } from "../providers/delete-dialog";
import { useSet } from "../hooks";
import SearchField from "./search-field";
import { usePopcorn } from "../providers/popcorn-provider";
import { DataTable, Column } from "./data-table";
import { VirtualizedCardGrid } from "./virtualized-card-grid";
import { ProjectTagChips } from "./project-tag-chips";

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
    onExport,
    onDelete,
  }: {
    project: Project;
    isSelected: boolean;
    onToggleSelection: () => void;
    onNavigate: () => void;
    onExport: () => void;
    onDelete: () => void;
  }) => {
    const isRecent =
      new Date(project.last_access).getTime() >
      Date.now() - 7 * 24 * 60 * 60 * 1000;

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
              {isRecent && (
                <Chip
                  label="Recent"
                  size="small"
                  color="success"
                  sx={{ height: 20, fontSize: "0.7rem" }}
                />
              )}
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

type ViewMode = "cards" | "table";

export default function ProjectsTable() {
  const api = useApi();
  const router = useRouter();
  const { data: projects, mutate } = api.get<Project[]>("projects");
  const selectedIds = useSet<number>([]);
  const [query, setQuery] = useState("");
  const [viewMode, setViewMode] = useState<ViewMode>("table");
  const deleteDialog = useDeleteDialog();
  const { setMessage } = usePopcorn();

  // Get campaign info for all projects
  const projectIds = useMemo(
    () => (projects || []).map((p) => p.id),
    [projects]
  );
  const campaignInfo = useProjectCampaigns(projectIds);

  // Filter and sort projects
  const filteredProjects = useMemo(() => {
    if (!Array.isArray(projects)) return [];
    return projects
      ?.filter((project) => {
        const searchTerm = query.toLowerCase();

        // Search in project name
        if (project.name.toLowerCase().includes(searchTerm)) {
          return true;
        }

        // Search in project description
        if (project.description?.toLowerCase().includes(searchTerm)) {
          return true;
        }

        // Search in project tags
        if (Array.isArray(project.tags)) {
          const tagMatches = project.tags.some((tag) => {
            if (typeof tag === "object" && tag.text) {
              return tag.text.toLowerCase().includes(searchTerm);
            }
            return false;
          });
          if (tagMatches) {
            return true;
          }
        }

        return false;
      })
      .sort(
        (a, b) =>
          new Date(b.last_access).getTime() - new Date(a.last_access).getTime()
      );
  }, [projects, query]);

  // Handlers
  const deleteSelected = useCallback(() => {
    const selectedProjects = projects?.filter((project) =>
      selectedIds.has(project.id)
    );
    if (selectedProjects) deleteProjects(selectedProjects);
  }, [projects, selectedIds]);

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
  const tableColumns: Column<Project>[] = useMemo(
    () => [
      {
        key: "checkbox",
        label: "",
        width: 50,
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
              <Box
                sx={{
                  width: 32,
                  height: 32,
                  borderRadius: 1,
                  bgcolor: campaign ? "secondary.50" : "primary.50",
                  display: "flex",
                  alignItems: "center",
                  justifyContent: "center",
                  border: "1px solid",
                  borderColor: campaign ? "secondary.200" : "primary.200",
                  flexShrink: 0,
                }}
              >
                {campaign ? (
                  <CampaignIcon sx={{ color: "secondary.main", fontSize: 16 }} />
                ) : (
                  <Science sx={{ color: "primary.main", fontSize: 16 }} />
                )}
              </Box>
              <Box sx={{ minWidth: 0 }}>
                <Typography variant="body2" sx={{ fontWeight: 500 }} noWrap>
                  {project.name}
                </Typography>
                <Stack direction="row" spacing={0.5} flexWrap="wrap" useFlexGap>
                  {new Date(project.last_access).getTime() >
                    Date.now() - 7 * 24 * 60 * 60 * 1000 && (
                    <Chip
                      label="Recent"
                      size="small"
                      color="success"
                      sx={{ height: 16, fontSize: "0.65rem", mt: 0.5 }}
                    />
                  )}
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
        width: 100,
        render: (_, project) => (
          <Stack
            direction="row"
            spacing={0.5}
            justifyContent="flex-end"
            onClick={(e) => e.stopPropagation()}
          >
            <Tooltip title="Export project">
              <IconButton
                size="small"
                onClick={() => exportProject(project)}
                sx={{ color: "primary.main" }}
              >
                <Download fontSize="small" />
              </IconButton>
            </Tooltip>
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
    ],
    [selectedIds, campaignInfo, router]
  );

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
        onExport={() => exportProject(project)}
        onDelete={() => deleteProjects([project])}
      />
    ),
    [selectedIds, router]
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
      {/* Header with search, view toggle, and selection actions */}
      <Box sx={{ mb: 2, flexShrink: 0 }}>
        {selectedIds.size === 0 ? (
          <Box>
            {/* Title and Controls Row */}
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                justifyContent: "space-between",
                mb: 2,
              }}
            >
              <Typography
                variant="h5"
                component="h2"
                sx={{ fontWeight: 600, color: "text.primary" }}
              >
                Your Projects
              </Typography>
              <Stack direction="row" spacing={2} alignItems="center">
                <SearchField what="projects" onDelay={setQuery} />
                <ToggleButtonGroup
                  value={viewMode}
                  exclusive
                  onChange={(_, newView) => newView && setViewMode(newView)}
                  size="small"
                  sx={{ height: 36 }}
                >
                  <ToggleButton value="cards" aria-label="card view">
                    <Tooltip title="Card view">
                      <ViewModule />
                    </Tooltip>
                  </ToggleButton>
                  <ToggleButton value="table" aria-label="table view">
                    <Tooltip title="Table view">
                      <ViewList />
                    </Tooltip>
                  </ToggleButton>
                </ToggleButtonGroup>
              </Stack>
            </Box>

            {/* Project count and select all */}
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                justifyContent: "space-between",
                mb: 2,
              }}
            >
              <Box sx={{ display: "flex", alignItems: "center", gap: 2 }}>
                <Tooltip
                  title={
                    selectedIds.size == filteredProjects.length
                      ? "Deselect all projects"
                      : "Select all projects"
                  }
                >
                  <Chip
                    icon={
                      <Checkbox
                        size="small"
                        checked={
                          selectedIds.size == filteredProjects.length &&
                          filteredProjects.length > 0
                        }
                        indeterminate={
                          selectedIds.size > 0 &&
                          selectedIds.size < filteredProjects.length
                        }
                        onClick={toggleAll}
                      />
                    }
                    label={`${filteredProjects.length} projects`}
                    variant="outlined"
                    clickable
                    onClick={toggleAll}
                    sx={{ "& .MuiChip-icon": { mr: 0.5 }, borderRadius: 2 }}
                  />
                </Tooltip>
              </Box>
            </Box>
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
              <Typography
                color="primary.main"
                variant="subtitle1"
                sx={{ fontWeight: 600, flexGrow: 1 }}
              >
                {selectedIds.size} project{selectedIds.size !== 1 ? "s" : ""}{" "}
                selected
              </Typography>
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

      {/* Render based on view mode */}
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
            hideHeader
            fillHeight
            estimateRowHeight={60}
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
      </Box>
    </Box>
  ) : (
    <Skeleton variant="rectangular" height={400} />
  );
}
