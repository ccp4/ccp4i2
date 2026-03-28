/*
 * Copyright (C) 2026 Newcastle University
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
import React, { useState, useCallback, useMemo } from "react";
import {
  Box,
  Paper,
  Typography,
  List,
  ListItem,
  ListItemButton,
  ListItemText,
  ListItemIcon,
  Divider,
  Chip,
  Stack,
  CircularProgress,
  Alert,
  IconButton,
  Breadcrumbs,
  Link,
  TextField,
  InputAdornment,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
} from "@mui/material";
import {
  FolderOutlined,
  WorkOutline,
  InsertDriveFileOutlined,
  ChevronRight,
  ArrowBack,
  Home,
  Search,
  Clear,
  Close,
} from "@mui/icons-material";
import { useApi } from "../../../api";
import {
  Project,
  Job,
  File as DjangoFile,
} from "../../../types/models";
import { useTheme } from "../../../theme/theme-provider";
import { useProjectJobs } from "../../../utils";

// ---------------------------------------------------------------------------
// Shared sub-components (compact list items identical to hierarchy browser)
// ---------------------------------------------------------------------------

interface CustomColors {
  ui: Record<string, string>;
}

const ProjectItem: React.FC<{
  project: Project;
  onSelect: (p: Project) => void;
  customColors: CustomColors;
}> = ({ project, onSelect, customColors }) => (
  <ListItem disablePadding sx={{ py: 0 }}>
    <ListItemButton
      onClick={() => onSelect(project)}
      sx={{
        py: 0.5,
        px: 1,
        "&:hover": { backgroundColor: customColors.ui.veryLightGray },
      }}
    >
      <ListItemIcon sx={{ minWidth: 28 }}>
        <FolderOutlined color="primary" fontSize="small" />
      </ListItemIcon>
      <ListItemText
        primary={
          <Typography variant="body2" sx={{ fontWeight: 500, lineHeight: 1.2 }}>
            {project.name || `Project ${project.id}`}
          </Typography>
        }
        secondary={
          <Typography
            variant="caption"
            color="text.secondary"
            sx={{ lineHeight: 1.2 }}
          >
            ID: {project.id} •{" "}
            {new Date(project.creation_time).toLocaleDateString()}
          </Typography>
        }
        sx={{ my: 0 }}
      />
      <ChevronRight color="action" fontSize="small" />
    </ListItemButton>
  </ListItem>
);

const getStatusColor = (status: number) => {
  switch (status) {
    case 1:
      return "warning";
    case 2:
      return "info";
    case 6:
      return "success";
    case 5:
      return "error";
    default:
      return "default";
  }
};

const getStatusText = (status: number) => {
  switch (status) {
    case 1:
      return "Pending";
    case 2:
      return "Running";
    case 6:
      return "Finished";
    case 5:
      return "Failed";
    default:
      return "Unknown";
  }
};

const JobItem: React.FC<{
  job: Job;
  onSelect: (j: Job) => void;
  customColors: CustomColors;
}> = ({ job, onSelect, customColors }) => (
  <ListItem disablePadding sx={{ py: 0 }}>
    <ListItemButton
      onClick={() => onSelect(job)}
      sx={{
        py: 0.5,
        px: 1,
        "&:hover": { backgroundColor: customColors.ui.veryLightGray },
      }}
    >
      <ListItemIcon sx={{ minWidth: 28 }}>
        <WorkOutline color="primary" fontSize="small" />
      </ListItemIcon>
      <ListItemText
        primary={
          <Stack direction="row" alignItems="center" spacing={0.5}>
            <Typography
              variant="body2"
              sx={{ fontWeight: 500, lineHeight: 1.2 }}
            >
              {job.title || `Job ${job.number}`}
            </Typography>
            <Chip
              size="small"
              label={getStatusText(job.status)}
              color={getStatusColor(job.status) as any}
              variant="outlined"
              sx={{
                height: 16,
                fontSize: "0.625rem",
                "& .MuiChip-label": { px: 0.5 },
              }}
            />
          </Stack>
        }
        secondary={
          <Typography
            variant="caption"
            color="text.secondary"
            sx={{ lineHeight: 1.2 }}
          >
            #{job.number} • {job.task_name || "Unknown task"}
          </Typography>
        }
        sx={{ my: 0 }}
      />
      <ChevronRight color="action" fontSize="small" />
    </ListItemButton>
  </ListItem>
);

const getFileTypeColor = (type: string) => {
  if (type.includes("pdb")) return "primary";
  if (type.includes("mtz")) return "secondary";
  if (type.includes("cif")) return "info";
  return "default";
};

const FileItem: React.FC<{
  file: DjangoFile;
  onSelect: (f: DjangoFile) => void;
  customColors: CustomColors;
}> = ({ file, onSelect, customColors }) => (
  <ListItem disablePadding sx={{ py: 0 }}>
    <ListItemButton
      onClick={() => onSelect(file)}
      sx={{
        py: 0.5,
        px: 1,
        "&:hover": { backgroundColor: customColors.ui.veryLightGray },
      }}
    >
      <ListItemIcon sx={{ minWidth: 28 }}>
        <InsertDriveFileOutlined color="primary" fontSize="small" />
      </ListItemIcon>
      <ListItemText
        primary={
          <Stack direction="row" alignItems="center" spacing={0.5}>
            <Typography
              variant="body2"
              sx={{ fontWeight: 500, lineHeight: 1.2 }}
            >
              {file.annotation || file.name}
            </Typography>
            <Chip
              size="small"
              label={file.type.split("/").pop() || "unknown"}
              color={getFileTypeColor(file.type) as any}
              variant="outlined"
              sx={{
                height: 16,
                fontSize: "0.625rem",
                "& .MuiChip-label": { px: 0.5 },
              }}
            />
          </Stack>
        }
        secondary={
          <Typography
            variant="caption"
            color="text.secondary"
            sx={{ lineHeight: 1.2 }}
          >
            {file.job_param_name || file.name}
          </Typography>
        }
        sx={{ my: 0 }}
      />
    </ListItemButton>
  </ListItem>
);

// ---------------------------------------------------------------------------
// Compact panel wrapper with search
// ---------------------------------------------------------------------------

const PanelHeader: React.FC<{
  title: string;
  onBack?: () => void;
  breadcrumbs?: React.ReactNode;
  searchValue: string;
  onSearchChange: (v: string) => void;
  searchPlaceholder?: string;
  customColors: CustomColors;
  children: React.ReactNode;
}> = ({
  title,
  onBack,
  breadcrumbs,
  searchValue,
  onSearchChange,
  searchPlaceholder = "Search...",
  customColors,
  children,
}) => (
  <Paper
    elevation={0}
    sx={{
      height: "100%",
      display: "flex",
      flexDirection: "column",
      border: `1px solid ${customColors.ui.mediumGray}`,
    }}
  >
    {/* Header */}
    <Box
      sx={{
        px: 1.5,
        py: 0.75,
        backgroundColor: customColors.ui.lightBlue,
        color: "white",
        borderBottom: `1px solid ${customColors.ui.mediumGray}`,
        display: "flex",
        alignItems: "center",
        gap: 0.5,
        minHeight: "36px",
      }}
    >
      {onBack && (
        <IconButton
          onClick={onBack}
          sx={{ color: "white", p: 0.25 }}
          size="small"
        >
          <ArrowBack fontSize="small" />
        </IconButton>
      )}
      <Box sx={{ flex: 1, minWidth: 0 }}>
        <Typography
          variant="body2"
          component="h2"
          sx={{ fontWeight: 600, lineHeight: 1.2 }}
        >
          {title}
        </Typography>
        {breadcrumbs && <Box sx={{ mt: 0.25 }}>{breadcrumbs}</Box>}
      </Box>
    </Box>

    {/* Search */}
    <Box sx={{ p: 1, borderBottom: `1px solid ${customColors.ui.mediumGray}` }}>
      <TextField
        fullWidth
        size="small"
        placeholder={searchPlaceholder}
        value={searchValue}
        onChange={(e) => onSearchChange(e.target.value)}
        InputProps={{
          startAdornment: (
            <InputAdornment position="start">
              <Search sx={{ fontSize: "1rem" }} color="action" />
            </InputAdornment>
          ),
          endAdornment: searchValue ? (
            <InputAdornment position="end">
              <IconButton
                size="small"
                onClick={() => onSearchChange("")}
                edge="end"
                sx={{ p: 0.25 }}
              >
                <Clear fontSize="small" />
              </IconButton>
            </InputAdornment>
          ) : null,
          sx: {
            fontSize: "0.875rem",
            "& .MuiInputBase-input": { py: 0.5 },
          },
        }}
        sx={{
          "& .MuiOutlinedInput-root": { backgroundColor: "white" },
        }}
      />
    </Box>

    {/* Content */}
    <Box sx={{ flex: 1, overflow: "auto" }}>{children}</Box>
  </Paper>
);

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

const sortJobsByNumber = (jobs: Job[]): Job[] =>
  [...jobs].sort((a, b) => {
    const first = (n: string) => parseInt(n.split(".")[0], 10) || 0;
    return first(b.number) - first(a.number);
  });

const renderLoading = () => (
  <Box sx={{ display: "flex", justifyContent: "center", p: 1.5 }}>
    <CircularProgress size={24} />
  </Box>
);

const renderError = (error: any) => (
  <Box sx={{ p: 1 }}>
    <Alert severity="error" sx={{ fontSize: "0.875rem", py: 0.5 }}>
      Failed to load data: {error?.message || "Unknown error"}
    </Alert>
  </Box>
);

const renderEmpty = (message: string) => (
  <Box sx={{ p: 1.5, textAlign: "center" }}>
    <Typography variant="caption" color="text.secondary">
      {message}
    </Typography>
  </Box>
);

// ---------------------------------------------------------------------------
// Main dialog
// ---------------------------------------------------------------------------

export interface BrowseProjectFilesDialogProps {
  open: boolean;
  onClose: () => void;
  /** Called when the user picks a file. Receives the full DjangoFile object. */
  onFileSelect: (file: DjangoFile) => void;
  /** MIME type names to include (e.g. ["chemical/x-pdb"]). If null/empty, all files shown. */
  allowedTypes?: string[] | null;
  /** Content flags required (e.g. [1,2] for anomalous). Null = no filter. */
  requiredContentFlag?: number[] | null;
  /** Label shown in the dialog title */
  fileTypeLabel?: string;
}

/** Content-flag conversion table (mirrors CAN_CONVERT_TO from CCP4XtalData). */
const CAN_CONVERT_TO: Record<number, number[]> = {
  1: [1, 2, 3, 4],
  2: [2, 4],
  3: [3, 4],
  4: [4],
};

const canConvertToRequired = (
  fileContent: number | null | undefined,
  requiredFlags: number[]
): boolean => {
  if (fileContent == null) return false;
  const convertibleTo = CAN_CONVERT_TO[fileContent];
  if (!convertibleTo) return false;
  return requiredFlags.some((r) => convertibleTo.includes(r));
};

export const BrowseProjectFilesDialog: React.FC<
  BrowseProjectFilesDialogProps
> = ({
  open,
  onClose,
  onFileSelect,
  allowedTypes,
  requiredContentFlag,
  fileTypeLabel = "file",
}) => {
  const { customColors } = useTheme();
  const api = useApi();

  // Navigation state
  const [selectedProject, setSelectedProject] = useState<Project | null>(null);
  const [selectedJob, setSelectedJob] = useState<Job | null>(null);

  // Search state
  const [projectSearch, setProjectSearch] = useState("");
  const [jobSearch, setJobSearch] = useState("");
  const [fileSearch, setFileSearch] = useState("");

  // Data fetching
  const { data: projects, isLoading: projectsLoading, error: projectsError } =
    api.get<Project[]>(open ? "/projects" : null);

  const { jobs } = useProjectJobs(selectedProject?.id, 10000);
  const jobsLoading = !jobs && !!selectedProject;

  const { data: files, isLoading: filesLoading, error: filesError } =
    api.get_endpoint<DjangoFile[]>(
      selectedJob
        ? { type: "jobs", id: selectedJob.id, endpoint: "files" }
        : null
    );

  // ----- Filtered lists -----

  const filteredProjects = useMemo(() => {
    if (!projects) return projects;
    if (!projectSearch.trim()) return projects;
    const q = projectSearch.toLowerCase();
    return projects.filter((p) => {
      const name = (p.name || `Project ${p.id}`).toLowerCase();
      return name.includes(q) || p.id.toString().includes(q);
    });
  }, [projects, projectSearch]);

  const filteredJobs = useMemo(() => {
    if (!jobs) return jobs;
    const parentJobs = sortJobsByNumber(jobs.filter((j) => j.parent === null));
    if (!jobSearch.trim()) return parentJobs;
    const q = jobSearch.toLowerCase();
    return parentJobs.filter((j) => {
      const title = (j.title || `Job ${j.number}`).toLowerCase();
      return (
        title.includes(q) ||
        j.number.toLowerCase().includes(q) ||
        (j.task_name || "").toLowerCase().includes(q)
      );
    });
  }, [jobs, jobSearch]);

  const filteredFiles = useMemo(() => {
    if (!files) return files;
    let result = files;

    // Filter by allowed MIME types
    if (allowedTypes && allowedTypes.length > 0) {
      result = result.filter(
        (f) =>
          allowedTypes.includes(f.type) || allowedTypes.includes("Unknown")
      );
    }

    // Filter by content flag convertibility
    if (requiredContentFlag && requiredContentFlag.length > 0) {
      result = result.filter((f) =>
        canConvertToRequired(f.content, requiredContentFlag)
      );
    }

    if (!fileSearch.trim()) return result;
    const q = fileSearch.toLowerCase();
    return result.filter(
      (f) =>
        f.name.toLowerCase().includes(q) ||
        f.annotation.toLowerCase().includes(q) ||
        f.type.toLowerCase().includes(q)
    );
  }, [files, allowedTypes, requiredContentFlag, fileSearch]);

  // ----- Navigation handlers -----

  const handleProjectSelect = useCallback((project: Project) => {
    setSelectedProject(project);
    setSelectedJob(null);
    setJobSearch("");
    setFileSearch("");
  }, []);

  const handleJobSelect = useCallback((job: Job) => {
    setSelectedJob(job);
    setFileSearch("");
  }, []);

  const handleFileClick = useCallback(
    (file: DjangoFile) => {
      onFileSelect(file);
      onClose();
    },
    [onFileSelect, onClose]
  );

  const handleBackToProjects = useCallback(() => {
    setSelectedProject(null);
    setSelectedJob(null);
    setJobSearch("");
    setFileSearch("");
  }, []);

  const handleBackToJobs = useCallback(() => {
    setSelectedJob(null);
    setFileSearch("");
  }, []);

  // Reset navigation when dialog reopens
  const handleClose = useCallback(() => {
    setSelectedProject(null);
    setSelectedJob(null);
    setProjectSearch("");
    setJobSearch("");
    setFileSearch("");
    onClose();
  }, [onClose]);

  // ----- Breadcrumbs -----

  const breadcrumbs = useMemo(() => {
    if (selectedJob) {
      return (
        <Breadcrumbs
          separator="›"
          sx={{ color: "rgba(255,255,255,0.8)", fontSize: "0.75rem" }}
        >
          <Link
            component="button"
            variant="caption"
            onClick={handleBackToProjects}
            sx={{
              color: "rgba(255,255,255,0.8)",
              textDecoration: "underline",
              display: "flex",
              alignItems: "center",
            }}
          >
            <Home sx={{ mr: 0.25, fontSize: "0.875rem" }} />
            Projects
          </Link>
          <Link
            component="button"
            variant="caption"
            onClick={handleBackToJobs}
            sx={{
              color: "rgba(255,255,255,0.8)",
              textDecoration: "underline",
            }}
          >
            {selectedProject?.name || `Project ${selectedProject?.id}`}
          </Link>
          <Typography variant="caption" sx={{ color: "white" }}>
            Job {selectedJob.number}
          </Typography>
        </Breadcrumbs>
      );
    }
    if (selectedProject) {
      return (
        <Breadcrumbs
          separator="›"
          sx={{ color: "rgba(255,255,255,0.8)", fontSize: "0.75rem" }}
        >
          <Link
            component="button"
            variant="caption"
            onClick={handleBackToProjects}
            sx={{
              color: "rgba(255,255,255,0.8)",
              textDecoration: "underline",
              display: "flex",
              alignItems: "center",
            }}
          >
            <Home sx={{ mr: 0.25, fontSize: "0.875rem" }} />
            Projects
          </Link>
          <Typography variant="caption" sx={{ color: "white" }}>
            {selectedProject.name || `Project ${selectedProject.id}`}
          </Typography>
        </Breadcrumbs>
      );
    }
    return null;
  }, [selectedProject, selectedJob, handleBackToProjects, handleBackToJobs]);

  // ----- Render -----

  const renderPanel = () => {
    // Files panel
    if (selectedJob) {
      return (
        <PanelHeader
          title={`Files in Job ${selectedJob.number}`}
          onBack={handleBackToJobs}
          breadcrumbs={breadcrumbs}
          searchValue={fileSearch}
          onSearchChange={setFileSearch}
          searchPlaceholder="Search files..."
          customColors={customColors}
        >
          {filesLoading && renderLoading()}
          {filesError && renderError(filesError)}
          {filteredFiles && filteredFiles.length === 0 && renderEmpty(
            fileSearch
              ? `No matching ${fileTypeLabel} files found`
              : `No compatible ${fileTypeLabel} files in this job`
          )}
          {filteredFiles && filteredFiles.length > 0 && (
            <List sx={{ py: 0 }}>
              {filteredFiles.map((file, i) => (
                <React.Fragment key={file.id}>
                  <FileItem
                    file={file}
                    onSelect={handleFileClick}
                    customColors={customColors}
                  />
                  {i < filteredFiles.length - 1 && <Divider sx={{ my: 0 }} />}
                </React.Fragment>
              ))}
            </List>
          )}
        </PanelHeader>
      );
    }

    // Jobs panel
    if (selectedProject) {
      return (
        <PanelHeader
          title={`Jobs in ${selectedProject.name || `Project ${selectedProject.id}`}`}
          onBack={handleBackToProjects}
          breadcrumbs={breadcrumbs}
          searchValue={jobSearch}
          onSearchChange={setJobSearch}
          searchPlaceholder="Search jobs..."
          customColors={customColors}
        >
          {jobsLoading && renderLoading()}
          {filteredJobs && filteredJobs.length === 0 && renderEmpty(
            jobSearch ? `No jobs match "${jobSearch}"` : "No jobs in this project"
          )}
          {filteredJobs && filteredJobs.length > 0 && (
            <List sx={{ py: 0 }}>
              {filteredJobs.map((job, i) => (
                <React.Fragment key={job.id}>
                  <JobItem
                    job={job}
                    onSelect={handleJobSelect}
                    customColors={customColors}
                  />
                  {i < filteredJobs.length - 1 && <Divider sx={{ my: 0 }} />}
                </React.Fragment>
              ))}
            </List>
          )}
        </PanelHeader>
      );
    }

    // Projects panel
    return (
      <PanelHeader
        title="Projects"
        searchValue={projectSearch}
        onSearchChange={setProjectSearch}
        searchPlaceholder="Search projects..."
        customColors={customColors}
      >
        {projectsLoading && renderLoading()}
        {projectsError && renderError(projectsError)}
        {filteredProjects && filteredProjects.length === 0 && renderEmpty(
          projectSearch
            ? `No projects match "${projectSearch}"`
            : "No projects found"
        )}
        {filteredProjects && filteredProjects.length > 0 && (
          <List sx={{ py: 0 }}>
            {filteredProjects.map((project, i) => (
              <React.Fragment key={project.id}>
                <ProjectItem
                  project={project}
                  onSelect={handleProjectSelect}
                  customColors={customColors}
                />
                {i < filteredProjects.length - 1 && (
                  <Divider sx={{ my: 0 }} />
                )}
              </React.Fragment>
            ))}
          </List>
        )}
      </PanelHeader>
    );
  };

  return (
    <Dialog
      open={open}
      onClose={handleClose}
      maxWidth="sm"
      fullWidth
      PaperProps={{ sx: { height: "70vh", maxHeight: 600 } }}
    >
      <DialogTitle
        sx={{
          display: "flex",
          justifyContent: "space-between",
          alignItems: "center",
          py: 1,
          px: 2,
        }}
      >
        <Typography variant="subtitle1" sx={{ fontWeight: 600 }}>
          Browse {fileTypeLabel} files from other projects
        </Typography>
        <IconButton onClick={handleClose} size="small" edge="end">
          <Close fontSize="small" />
        </IconButton>
      </DialogTitle>

      <DialogContent sx={{ p: 0, display: "flex", flexDirection: "column" }}>
        <Box sx={{ flex: 1, minHeight: 0, p: 1 }}>{renderPanel()}</Box>
      </DialogContent>

      <DialogActions sx={{ px: 2, py: 1 }}>
        <Button onClick={handleClose} size="small">
          Cancel
        </Button>
      </DialogActions>
    </Dialog>
  );
};
