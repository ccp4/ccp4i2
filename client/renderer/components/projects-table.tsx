"use client";
import React, { useMemo, useState } from "react";
import { useRouter } from "next/navigation";
import {
  Box,
  Card,
  CardContent,
  Checkbox,
  Chip,
  Grid,
  IconButton,
  LinearProgress,
  Paper,
  Skeleton,
  Stack,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  Theme,
  Tooltip,
  Typography,
  ToggleButton,
  ToggleButtonGroup,
  Pagination,
  TableContainer,
} from "@mui/material";
import {
  Clear,
  Delete,
  Download,
  FolderOpen,
  Schedule,
  Science,
  StarBorder,
  Star,
  ViewModule,
  ViewList,
} from "@mui/icons-material";
import { alpha } from "@mui/material/styles";
import { useApi } from "../api";
import { Project, ProjectTag } from "../types/models";
import { shortDate } from "../pipes";
import { useDeleteDialog } from "../providers/delete-dialog";
import { useSet } from "../hooks";
import SearchField from "./search-field";
import { usePopcorn } from "../providers/popcorn-provider";

// Component to display project tag chips
const ProjectTagChips = React.memo(
  ({
    project,
    maxVisible = 3,
    size = "small" as "small" | "medium",
  }: {
    project: Project;
    maxVisible?: number;
    size?: "small" | "medium";
  }) => {
    // Handle both old format (number[]) and new format (ProjectTag[])
    const projectTagsData = Array.isArray(project.tags)
      ? project.tags.filter(
          (tag): tag is ProjectTag => typeof tag === "object" && tag !== null
        )
      : [];

    if (projectTagsData.length === 0) {
      // Return empty container to maintain consistent layout
      return (
        <Box
          sx={{
            minHeight: size === "small" ? 20 : 24,
            display: "flex",
            alignItems: "center",
          }}
        >
          <Typography
            variant="caption"
            color="text.disabled"
            sx={{ fontSize: size === "small" ? "0.65rem" : "0.7rem" }}
          >
            No tags
          </Typography>
        </Box>
      );
    }

    const visibleTags = projectTagsData.slice(0, maxVisible);
    const hiddenCount = projectTagsData.length - maxVisible;

    return (
      <Box
        sx={{
          display: "flex",
          flexWrap: "wrap",
          gap: 0.5,
          mt: size === "small" ? 0.5 : 1,
          minHeight: size === "small" ? 20 : 24, // Maintain consistent height even when no tags
        }}
      >
        {visibleTags.map((tag) => (
          <Chip
            key={tag.id}
            label={tag.text}
            size={size}
            variant="outlined"
            sx={{
              height: size === "small" ? 20 : 24,
              fontSize: size === "small" ? "0.7rem" : "0.75rem",
              bgcolor: "primary.50",
              borderColor: "primary.200",
              color: "primary.700",
              fontWeight: 500,
              "&:hover": {
                bgcolor: "primary.100",
                borderColor: "primary.300",
              },
            }}
          />
        ))}
        {hiddenCount > 0 && (
          <Chip
            label={`+${hiddenCount}`}
            size={size}
            variant="outlined"
            sx={{
              height: size === "small" ? 20 : 24,
              fontSize: size === "small" ? "0.7rem" : "0.75rem",
              bgcolor: "grey.100",
              borderColor: "grey.300",
              color: "grey.600",
              fontWeight: 500,
            }}
          />
        )}
      </Box>
    );
  }
);

ProjectTagChips.displayName = "ProjectTagChips";

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
          <ProjectTagChips project={project} maxVisible={2} size="small" />

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

// Memoized TableRow component for performance
const ProjectTableRow = React.memo(
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
  }) => (
    <TableRow
      hover
      onClick={onNavigate}
      sx={{ cursor: "pointer", ...(isSelected && sxSelected) }}
    >
      <TableCell padding="checkbox">
        <Checkbox
          checked={isSelected}
          onChange={(event) => {
            event.stopPropagation();
            onToggleSelection();
          }}
        />
      </TableCell>
      <TableCell>
        <Box sx={{ display: "flex", alignItems: "center", gap: 1.5 }}>
          <Box
            sx={{
              width: 32,
              height: 32,
              borderRadius: 1,
              bgcolor: "primary.50",
              display: "flex",
              alignItems: "center",
              justifyContent: "center",
              border: "1px solid",
              borderColor: "primary.200",
            }}
          >
            <Science sx={{ color: "primary.main", fontSize: 16 }} />
          </Box>
          <Box>
            <Typography variant="body2" sx={{ fontWeight: 500 }}>
              {project.name}
            </Typography>
            {new Date(project.last_access).getTime() >
              Date.now() - 7 * 24 * 60 * 60 * 1000 && (
              <Chip
                label="Recent"
                size="small"
                color="success"
                sx={{ height: 16, fontSize: "0.65rem", mt: 0.5 }}
              />
            )}
          </Box>
        </Box>
      </TableCell>
      <TableCell>
        <ProjectTagChips project={project} maxVisible={3} size="small" />
      </TableCell>
      <TableCell>
        <Typography variant="body2" color="text.secondary">
          {shortDate(project.creation_time)}
        </Typography>
      </TableCell>
      <TableCell>
        <Typography variant="body2" color="text.secondary">
          {shortDate(project.last_access)}
        </Typography>
      </TableCell>
      <TableCell align="right">
        <Stack direction="row" spacing={0.5} justifyContent="flex-end">
          <Tooltip title="Export project">
            <IconButton
              size="small"
              onClick={(event) => {
                event.stopPropagation();
                onExport();
              }}
              sx={{ color: "primary.main" }}
            >
              <Download fontSize="small" />
            </IconButton>
          </Tooltip>
          <Tooltip title="Delete project">
            <IconButton
              size="small"
              onClick={(event) => {
                event.stopPropagation();
                onDelete();
              }}
              sx={{ color: "error.main" }}
            >
              <Delete fontSize="small" />
            </IconButton>
          </Tooltip>
        </Stack>
      </TableCell>
    </TableRow>
  )
);

const sxSelected = {
  bgcolor: (theme: Theme) =>
    alpha(theme.palette.primary.main, theme.palette.action.activatedOpacity),
};

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

const ITEMS_PER_PAGE = 50; // For performance with large datasets

export default function ProjectsTable() {
  const api = useApi();
  const router = useRouter();
  const { data: projects, mutate } = api.get<Project[]>("projects");
  const selectedIds = useSet<number>([]);
  const [query, setQuery] = useState("");
  const [viewMode, setViewMode] = useState<ViewMode>("cards");
  const [currentPage, setCurrentPage] = useState(1);
  const deleteDialog = useDeleteDialog();
  const { setMessage } = usePopcorn();

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
            // Handle both legacy (number[]) and enhanced (ProjectTag[]) formats
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

  const paginatedProjects = useMemo(() => {
    const startIndex = (currentPage - 1) * ITEMS_PER_PAGE;
    const endIndex = startIndex + ITEMS_PER_PAGE;
    return filteredProjects.slice(startIndex, endIndex);
  }, [filteredProjects, currentPage]);

  const totalPages = Math.ceil(filteredProjects.length / ITEMS_PER_PAGE);

  // Reset to first page when search changes
  useMemo(() => {
    setCurrentPage(1);
  }, [query]);

  function deleteSelected() {
    const selectedProjects = projects?.filter((project) =>
      selectedIds.has(project.id)
    );
    if (selectedProjects) deleteProjects(selectedProjects);
  }

  function deleteProjects(projects: Project[]) {
    if (deleteDialog)
      deleteDialog({
        type: "show",
        what:
          projects.length === 1
            ? projects[0].name
            : `${projects.length} projects`,
        onDelete: () => {
          const promises = projects.map((project) => {
            selectedIds.delete(project.id);
            return api.delete(`projects/${project.id}`);
          });
          Promise.all(promises).then(() => mutate());
        },
      });
  }

  async function exportProject(project: Project) {
    try {
      // Start the export process by calling the API endpoint
      const exportResult: any = await api.post(`projects/${project.id}/export/`, {});

      if (exportResult?.success === false) {
        setMessage(`Failed to export "${project.name}": ${exportResult?.error || "Unknown error"}`, "error");
        return;
      }

      // Show success notification
      console.log(
        `Export started for project "${project.name}":`,
        exportResult
      );
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

      // Start exports for all selected projects
      const exportPromises = selectedProjects.map(async (project) => {
        try {
          const result: any = await api.post(`projects/${project.id}/export/`, {});
          if (result?.success === false) {
            return { project: project.name, success: false, error: result?.error };
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

      // Show results to user
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

      console.log("Bulk export results:", results);
    } catch (error) {
      console.error("Failed to export selected projects:", error);
      setMessage(`Failed to start exports: ${error instanceof Error ? error.message : String(error)}`, "error");
    }
  }

  function toggleAll() {
    if (projects) {
      if (selectedIds.size === projects.length) {
        selectedIds.clear();
      } else {
        projects.forEach((project) => selectedIds.add(project.id));
      }
    }
  }

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
    <Box>
      {/* Header with search, view toggle, and selection actions */}
      <Box sx={{ mb: 3 }}>
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

                {filteredProjects.length > ITEMS_PER_PAGE && (
                  <Typography variant="body2" color="text.secondary">
                    Showing {Math.min(ITEMS_PER_PAGE, filteredProjects.length)}{" "}
                    of {filteredProjects.length}
                  </Typography>
                )}
              </Box>

              {/* Pagination for large datasets */}
              {totalPages > 1 && (
                <Pagination
                  count={totalPages}
                  page={currentPage}
                  onChange={(_, page) => setCurrentPage(page)}
                  size="small"
                  color="primary"
                />
              )}
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
      {viewMode === "cards" ? (
        // Card View - Better for visual browsing, limited by pagination for performance
        <Grid container spacing={3}>
          {paginatedProjects.map((project: Project) => (
            <Grid item xs={12} sm={6} md={4} lg={3} key={project.id}>
              <ProjectCard
                project={project}
                isSelected={selectedIds.has(project.id)}
                onToggleSelection={() => {
                  selectedIds.has(project.id)
                    ? selectedIds.delete(project.id)
                    : selectedIds.add(project.id);
                }}
                onNavigate={() => router.push(`/project/${project.id}`)}
                onExport={() => exportProject(project)}
                onDelete={() => deleteProjects([project])}
              />
            </Grid>
          ))}
        </Grid>
      ) : (
        // Table View - Better for large datasets, more compact
        <TableContainer component={Paper} elevation={1}>
          <Table size="small" stickyHeader>
            <TableHead>
              <TableRow>
                <TableCell padding="checkbox">
                  <Checkbox
                    checked={
                      selectedIds.size == paginatedProjects.length &&
                      paginatedProjects.length > 0
                    }
                    indeterminate={
                      selectedIds.size > 0 &&
                      selectedIds.size < paginatedProjects.length
                    }
                    onChange={toggleAll}
                  />
                </TableCell>
                <TableCell>Project Name</TableCell>
                <TableCell>Tags</TableCell>
                <TableCell>Created</TableCell>
                <TableCell>Last Accessed</TableCell>
                <TableCell align="right">Actions</TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              {paginatedProjects.map((project: Project) => (
                <ProjectTableRow
                  key={project.id}
                  project={project}
                  isSelected={selectedIds.has(project.id)}
                  onToggleSelection={() => {
                    selectedIds.has(project.id)
                      ? selectedIds.delete(project.id)
                      : selectedIds.add(project.id);
                  }}
                  onNavigate={() => router.push(`/project/${project.id}`)}
                  onExport={() => exportProject(project)}
                  onDelete={() => deleteProjects([project])}
                />
              ))}
            </TableBody>
          </Table>
        </TableContainer>
      )}

      {/* Bottom pagination for large datasets */}
      {totalPages > 1 && (
        <Box sx={{ display: "flex", justifyContent: "center", mt: 4 }}>
          <Pagination
            count={totalPages}
            page={currentPage}
            onChange={(_, page) => setCurrentPage(page)}
            color="primary"
            showFirstButton
            showLastButton
          />
        </Box>
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
  ) : (
    <Skeleton variant="rectangular" height={400} />
  );
}
