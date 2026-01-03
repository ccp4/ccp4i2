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
} from "@mui/icons-material";
import { useApi } from "../../api";
import { Project, Job, File as DjangoFile } from "../../types/models";
import { useTheme } from "../../theme/theme-provider";

interface CCP4i2HierarchyBrowserProps {
  onFileSelect: (fileId: number) => Promise<void>;
}

interface HierarchyPanelProps {
  title: string;
  children: React.ReactNode;
  onBack?: () => void;
  breadcrumbs?: React.ReactNode;
  searchValue: string;
  onSearchChange: (value: string) => void;
  searchPlaceholder?: string;
  customColors: typeof import("../../theme/palette").lightCustomColors;
}

const HierarchyPanel: React.FC<HierarchyPanelProps> = ({
  title,
  children,
  onBack,
  breadcrumbs,
  searchValue,
  onSearchChange,
  searchPlaceholder = "Search...",
  customColors,
}) => (
  <Paper
    elevation={1}
    sx={{
      height: "100%",
      display: "flex",
      flexDirection: "column",
      border: `1px solid ${customColors.ui.mediumGray}`,
    }}
  >
    {/* Compact Header */}
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

    {/* Compact Search Box */}
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
          endAdornment: searchValue && (
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
          ),
          sx: {
            fontSize: "0.875rem",
            "& .MuiInputBase-input": {
              py: 0.5,
            },
          },
        }}
        sx={{
          "& .MuiOutlinedInput-root": {
            backgroundColor: "white",
          },
        }}
      />
    </Box>

    {/* Content */}
    <Box sx={{ flex: 1, overflow: "auto" }}>{children}</Box>
  </Paper>
);

interface ProjectItemProps {
  project: Project;
  onSelect: (project: Project) => void;
  customColors: typeof import("../../theme/palette").lightCustomColors;
}

const ProjectItem: React.FC<ProjectItemProps> = ({
  project,
  onSelect,
  customColors,
}) => (
  <ListItem disablePadding sx={{ py: 0 }}>
    <ListItemButton
      onClick={() => onSelect(project)}
      sx={{
        py: 0.5,
        px: 1,
        "&:hover": {
          backgroundColor: customColors.ui.veryLightGray,
        },
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

interface JobItemProps {
  job: Job;
  onSelect: (job: Job) => void;
  customColors: typeof import("../../theme/palette").lightCustomColors;
}

const JobItem: React.FC<JobItemProps> = ({ job, onSelect, customColors }) => {
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

  return (
    <ListItem disablePadding sx={{ py: 0 }}>
      <ListItemButton
        onClick={() => onSelect(job)}
        sx={{
          py: 0.5,
          px: 1,
          "&:hover": {
            backgroundColor: customColors.ui.veryLightGray,
          },
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
                color={getStatusColor(job.status)}
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
};

interface FileItemProps {
  file: DjangoFile;
  onSelect?: (file: DjangoFile) => void;
  customColors: typeof import("../../theme/palette").lightCustomColors;
}

const FileItem: React.FC<FileItemProps> = ({
  file,
  onSelect,
  customColors,
}) => {
  const getFileTypeColor = (type: string) => {
    if (type.includes("pdb")) return "primary";
    if (type.includes("mtz")) return "secondary";
    if (type.includes("cif")) return "info";
    return "default";
  };

  return (
    <ListItem disablePadding sx={{ py: 0 }}>
      <ListItemButton
        onClick={() => onSelect?.(file)}
        sx={{
          py: 0.5,
          px: 1,
          "&:hover": {
            backgroundColor: customColors.ui.veryLightGray,
          },
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
                {file.name}
              </Typography>
              <Chip
                size="small"
                label={file.type.split("/").pop() || "unknown"}
                color={getFileTypeColor(file.type)}
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
              {file.job_param_name || "N/A"}
            </Typography>
          }
          sx={{ my: 0 }}
        />
      </ListItemButton>
    </ListItem>
  );
};

// Helper function to sort jobs by hierarchical number
const sortJobsByNumber = (jobs: Job[]): Job[] => {
  return [...jobs].sort((a, b) => {
    const getFirstIndex = (jobNumber: string): number => {
      const parts = jobNumber.split(".");
      return parseInt(parts[0], 10) || 0;
    };

    const aFirstIndex = getFirstIndex(a.number);
    const bFirstIndex = getFirstIndex(b.number);

    return bFirstIndex - aFirstIndex;
  });
};

const ALLOWED_FILE_TYPES = [
  "chemical/x-pdb",
  "application/CCP4-map",
  "application/CCP4-mtz-map",
  "application/refmac-dictionary",
] as const;

const isAllowedFileType = (fileType: string): boolean => {
  return ALLOWED_FILE_TYPES.includes(fileType as any);
};

// Constants for refresh intervals
const REFRESH_INTERVALS = {
  PROJECTS: 10000, // 10 seconds - projects may be created by other processes
  JOBS: 10000, // 10 seconds - jobs may be created/updated by other processes
  FILES_STATIC: 0, // No refresh for static files (finished jobs)
  FILES_DYNAMIC: 5000, // 5 seconds for running jobs or when job status might change
} as const;

export const CCP4i2HierarchyBrowser: React.FC<CCP4i2HierarchyBrowserProps> = ({
  onFileSelect,
}) => {
  const { customColors } = useTheme();
  const api = useApi();

  // State management
  const [selectedProject, setSelectedProject] = useState<Project | null>(null);
  const [selectedJob, setSelectedJob] = useState<Job | null>(null);

  // Search state
  const [projectSearchTerm, setProjectSearchTerm] = useState<string>("");
  const [jobSearchTerm, setJobSearchTerm] = useState<string>("");
  const [fileSearchTerm, setFileSearchTerm] = useState<string>("");

  // Calculate refresh interval for files based on job status
  const getFilesRefreshInterval = useCallback((job: Job | null): number => {
    if (!job) return REFRESH_INTERVALS.FILES_STATIC;

    // Refresh files if job is running (status 2) or pending (status 1)
    // Also refresh if job just finished to catch final files
    const isJobActive = job.status === 1 || job.status === 2;

    return isJobActive
      ? REFRESH_INTERVALS.FILES_DYNAMIC
      : REFRESH_INTERVALS.FILES_STATIC;
  }, []);

  // API calls with appropriate refresh intervals
  const {
    data: projects,
    isLoading: projectsLoading,
    error: projectsError,
  } = api.get<Project[]>("/projects", REFRESH_INTERVALS.PROJECTS);

  const {
    data: jobs,
    isLoading: jobsLoading,
    error: jobsError,
  } = api.get_endpoint<Job[]>(
    {
      type: "projects",
      id: selectedProject?.id,
      endpoint: "jobs",
    },
    REFRESH_INTERVALS.JOBS
  );

  const {
    data: files,
    isLoading: filesLoading,
    error: filesError,
  } = api.get_endpoint<DjangoFile[]>(
    {
      type: "jobs",
      id: selectedJob?.id,
      endpoint: "files",
    },
    getFilesRefreshInterval(selectedJob)
  );

  // Filtered and sorted data based on search terms
  const filteredProjects = useMemo(() => {
    if (!projects || !projectSearchTerm.trim()) return projects;

    const searchLower = projectSearchTerm.toLowerCase();
    return projects.filter((project) => {
      const projectName = (
        project.name || `Project ${project.id}`
      ).toLowerCase();
      const projectId = project.id.toString();
      return (
        projectName.includes(searchLower) || projectId.includes(searchLower)
      );
    });
  }, [projects, projectSearchTerm]);

  const filteredJobs = useMemo(() => {
    if (!jobs) return jobs;

    const parentJobs = jobs.filter((job) => job.parent === null);
    const sortedJobs = sortJobsByNumber(parentJobs);

    if (!jobSearchTerm.trim()) return sortedJobs;

    const searchLower = jobSearchTerm.toLowerCase();
    return sortedJobs.filter((job) => {
      const jobTitle = (job.title || `Job ${job.number}`).toLowerCase();
      const jobNumber = job.number.toLowerCase();
      const taskName = (job.task_name || "").toLowerCase();
      return (
        jobTitle.includes(searchLower) ||
        jobNumber.includes(searchLower) ||
        taskName.includes(searchLower)
      );
    });
  }, [jobs, jobSearchTerm]);

  const filteredFiles = useMemo(() => {
    if (!files) return files;

    const allowedFiles = files.filter((file) => isAllowedFileType(file.type));

    if (!fileSearchTerm.trim()) return allowedFiles;

    const searchLower = fileSearchTerm.toLowerCase();
    return allowedFiles.filter((file) => {
      const fileName = file.name.toLowerCase();
      const fileType = file.type.toLowerCase();
      return fileName.includes(searchLower) || fileType.includes(searchLower);
    });
  }, [files, fileSearchTerm]);

  // Event handlers
  const handleProjectSelect = useCallback((project: Project) => {
    setSelectedProject(project);
    setSelectedJob(null);
    setJobSearchTerm("");
    setFileSearchTerm("");
  }, []);

  const handleJobSelect = useCallback((job: Job) => {
    setSelectedJob(job);
    setFileSearchTerm("");
  }, []);

  const handleFileSelect = useCallback(
    (file: DjangoFile) => {
      console.log("File selected:", file);
      onFileSelect(file.id);
    },
    [onFileSelect]
  );

  const handleBackToProjects = useCallback(() => {
    setSelectedProject(null);
    setSelectedJob(null);
    setJobSearchTerm("");
    setFileSearchTerm("");
  }, []);

  const handleBackToJobs = useCallback(() => {
    setSelectedJob(null);
    setFileSearchTerm("");
  }, []);

  // Render loading state
  const renderLoading = () => (
    <Box sx={{ display: "flex", justifyContent: "center", p: 1.5 }}>
      <CircularProgress size={24} />
    </Box>
  );

  // Render error state
  const renderError = (error: any) => (
    <Box sx={{ p: 1 }}>
      <Alert severity="error" sx={{ fontSize: "0.875rem", py: 0.5 }}>
        Failed to load data: {error?.message || "Unknown error"}
      </Alert>
    </Box>
  );

  // Render empty state
  const renderEmpty = (message: string) => (
    <Box sx={{ p: 1.5, textAlign: "center" }}>
      <Typography variant="caption" color="text.secondary">
        {message}
      </Typography>
    </Box>
  );

  // Create breadcrumbs
  const createBreadcrumbs = () => {
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
    } else if (selectedProject) {
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
  };

  return (
    <Box sx={{ height: "100%" }}>
      {/* Show Files Panel if job is selected */}
      {selectedJob && (
        <HierarchyPanel
          title={`Files in Job ${selectedJob.number}`}
          onBack={handleBackToJobs}
          breadcrumbs={createBreadcrumbs()}
          searchValue={fileSearchTerm}
          onSearchChange={setFileSearchTerm}
          searchPlaceholder="Search files..."
          customColors={customColors}
        >
          {filesLoading && renderLoading()}
          {filesError && renderError(filesError)}
          {filteredFiles &&
            filteredFiles.length === 0 &&
            files &&
            files.filter((file) => isAllowedFileType(file.type)).length > 0 &&
            renderEmpty(`No supported files match "${fileSearchTerm}"`)}
          {filteredFiles &&
            filteredFiles.length === 0 &&
            (!files ||
              files.filter((file) => isAllowedFileType(file.type)).length ===
                0) &&
            renderEmpty("No supported files found")}
          {filteredFiles && filteredFiles.length > 0 && (
            <List sx={{ py: 0 }}>
              {filteredFiles.map((file, index) => (
                <React.Fragment key={file.id}>
                  <FileItem
                    file={file}
                    onSelect={handleFileSelect}
                    customColors={customColors}
                  />
                  {index < filteredFiles.length - 1 && (
                    <Divider sx={{ my: 0 }} />
                  )}
                </React.Fragment>
              ))}
            </List>
          )}
        </HierarchyPanel>
      )}

      {/* Show Jobs Panel if project is selected but no job */}
      {selectedProject && !selectedJob && (
        <HierarchyPanel
          title={`Jobs in ${
            selectedProject.name || `Project ${selectedProject.id}`
          }`}
          onBack={handleBackToProjects}
          breadcrumbs={createBreadcrumbs()}
          searchValue={jobSearchTerm}
          onSearchChange={setJobSearchTerm}
          searchPlaceholder="Search jobs..."
          customColors={customColors}
        >
          {jobsLoading && renderLoading()}
          {jobsError && renderError(jobsError)}
          {filteredJobs &&
            filteredJobs.length === 0 &&
            jobs &&
            jobs.length > 0 &&
            renderEmpty(`No jobs match "${jobSearchTerm}"`)}
          {filteredJobs &&
            filteredJobs.length === 0 &&
            (!jobs || jobs.length === 0) &&
            renderEmpty("No jobs found in this project")}
          {filteredJobs && filteredJobs.length > 0 && (
            <List sx={{ py: 0 }}>
              {filteredJobs.map((job, index) => (
                <React.Fragment key={job.id}>
                  <JobItem
                    job={job}
                    onSelect={handleJobSelect}
                    customColors={customColors}
                  />
                  {index < filteredJobs.length - 1 && (
                    <Divider sx={{ my: 0 }} />
                  )}
                </React.Fragment>
              ))}
            </List>
          )}
        </HierarchyPanel>
      )}

      {/* Show Projects Panel if nothing is selected */}
      {!selectedProject && (
        <HierarchyPanel
          title="Projects"
          searchValue={projectSearchTerm}
          onSearchChange={setProjectSearchTerm}
          searchPlaceholder="Search projects..."
          customColors={customColors}
        >
          {projectsLoading && renderLoading()}
          {projectsError && renderError(projectsError)}
          {filteredProjects &&
            filteredProjects.length === 0 &&
            projects &&
            projects.length > 0 &&
            renderEmpty(`No projects match "${projectSearchTerm}"`)}
          {filteredProjects &&
            filteredProjects.length === 0 &&
            (!projects || projects.length === 0) &&
            renderEmpty("No projects found")}
          {filteredProjects && filteredProjects.length > 0 && (
            <List sx={{ py: 0 }}>
              {filteredProjects.map((project, index) => (
                <React.Fragment key={project.id}>
                  <ProjectItem
                    project={project}
                    onSelect={handleProjectSelect}
                    customColors={customColors}
                  />
                  {index < filteredProjects.length - 1 && (
                    <Divider sx={{ my: 0 }} />
                  )}
                </React.Fragment>
              ))}
            </List>
          )}
        </HierarchyPanel>
      )}
    </Box>
  );
};
