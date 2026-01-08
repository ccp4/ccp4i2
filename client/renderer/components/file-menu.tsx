"use client";
import { useState } from "react";
import { useRouter } from "next/navigation";
import {
  Button,
  Divider,
  ListItemIcon,
  ListItemText,
  Menu,
  MenuItem,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Typography,
  Box,
} from "@mui/material";
import { Add, Download, Menu as MenuIcon, Upload } from "@mui/icons-material";
import { useApi } from "../api";
import { Project } from "../types/models";
import { useCCP4i2Window } from "../app-context";
import { apiPost } from "../api-fetch";
import { ProjectExportsDialog } from "./project-exports";

interface ExportResult {
  status: string;
  export_file_name: string;
  log_file_name: string;
  process_id: number;
}

export default function FileMenu() {
  const router = useRouter();
  const api = useApi();
  const { projectId } = useCCP4i2Window();
  const { data: projects, mutate: mutateProjects } =
    api.get<Project[]>("projects");

  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
  const open = Boolean(anchorEl);
  const [exportsDialogOpen, setExportsDialogOpen] = useState(false);
  const [exportSuccessDialogOpen, setExportSuccessDialogOpen] = useState(false);
  const [exportResult, setExportResult] = useState<ExportResult | null>(null);

  const handleClick = (event: React.MouseEvent<HTMLButtonElement>) => {
    setAnchorEl(event.currentTarget);
  };
  const handleClose = () => {
    setAnchorEl(null);
  };

  const handleBrowser = () => {
    const newWindow = window.open("/ccp4i2");
    setAnchorEl(null);
    // Check if the window was successfully opened
    if (newWindow) {
      console.log("New window opened successfully!");
    } else {
      console.log(
        "Failed to open new window. It might be blocked by a popup blocker."
      );
    }
  };

  const handleExportProject = () => {
    apiPost(`projects/${projectId}/export`, {})
      .then((result: ExportResult) => {
        console.log("Export started:", result);
        setExportResult(result);
        setExportSuccessDialogOpen(true);
        handleClose();
      })
      .catch((error) => {
        console.error("Export failed:", error);
        // You might want to show an error dialog here too
      });
  };

  const handleImportProject = () => {
    const newWindow = window.open("/ccp4i2/import-project");
    setAnchorEl(null);
    // Check if the window was successfully opened
    if (newWindow) {
      console.log("New window opened successfully!");
    } else {
      console.log(
        "Failed to open new window. It might be blocked by a popup blocker."
      );
    }
  };

  const handleOpenExportsDialog = () => {
    setExportsDialogOpen(true);
    handleClose();
  };

  return (
    <>
      <Button color="inherit" onClick={handleClick}>
        File/Projects
      </Button>
      <Menu anchorEl={anchorEl} open={open} onClose={handleClose}>
        <MenuItem
          key="Manage"
          onClick={() => {
            handleClose();
            router.push("/ccp4i2");
          }}
        >
          <ListItemIcon>
            <MenuIcon fontSize="small" />
          </ListItemIcon>
          <ListItemText>Manage/open projects</ListItemText>
        </MenuItem>
        <MenuItem
          key="Add"
          onClick={() => {
            handleClose();
            router.push("/ccp4i2/new-project");
          }}
        >
          <ListItemIcon>
            <Add fontSize="small" />
          </ListItemIcon>
          <ListItemText>New project</ListItemText>
        </MenuItem>
        <MenuItem key="Export" onClick={handleExportProject}>
          <ListItemIcon>
            <Download fontSize="small" />
          </ListItemIcon>
          <ListItemText>Export project</ListItemText>
        </MenuItem>
        <MenuItem key="Import" onClick={handleImportProject}>
          <ListItemIcon>
            <Upload fontSize="small" />
          </ListItemIcon>
          <ListItemText>Import project</ListItemText>
        </MenuItem>
        <MenuItem key="Exports" onClick={handleOpenExportsDialog}>
          <ListItemIcon>
            <Download fontSize="small" />
          </ListItemIcon>
          <ListItemText>Exports...</ListItemText>
        </MenuItem>
        <Divider />
        {Array.isArray(projects) &&
          projects
            .sort((a: Project, b: Project) => {
              const dateA = new Date(a.last_access);
              const dateB = new Date(b.last_access);
              return dateA.getTime() - dateB.getTime();
            })
            .slice(-10)
            .reverse()
            .map((project: Project) => (
              <MenuItem
                key={project.id}
                onClick={async () => {
                  setAnchorEl(null);
                  const formData = new FormData();
                  const nowString = new Date().toISOString();
                  formData.set("last_access", nowString);
                  const result = await api
                    .patch(`projects/${project.id}`, formData)
                    .then(() => {
                      mutateProjects();
                    });
                  router.push(`/ccp4i2/project/${project.id}`);
                }}
              >
                {project.name} - {`${new Date(project.last_access)}`}
              </MenuItem>
            ))}
        <MenuItem key="MoreProjects" onClick={handleClose}>
          More Projects
        </MenuItem>
        <Divider />
        <MenuItem key="CCP4i" onClick={handleClose}>
          View old CCP4i projects
        </MenuItem>
        <MenuItem key="Browser" onClick={handleBrowser}>
          Browser
        </MenuItem>
        <MenuItem key="Quit" onClick={handleClose}>
          Quit CCP4i2
        </MenuItem>
      </Menu>

      {/* Export Success Dialog */}
      <Dialog
        open={exportSuccessDialogOpen}
        onClose={() => setExportSuccessDialogOpen(false)}
        maxWidth="sm"
        fullWidth
      >
        <DialogTitle>Project Export Started</DialogTitle>
        <DialogContent>
          <Box sx={{ mb: 2 }}>
            <Typography variant="body1" gutterBottom>
              Your project export has been started successfully.
            </Typography>
          </Box>

          {exportResult && (
            <Box sx={{ mb: 2 }}>
              <Typography variant="body2" color="text.secondary" gutterBottom>
                <strong>Export File:</strong> {exportResult.export_file_name}
              </Typography>
              <Typography variant="body2" color="text.secondary" gutterBottom>
                <strong>Log File:</strong> {exportResult.log_file_name}
              </Typography>
              <Typography variant="body2" color="text.secondary" gutterBottom>
                <strong>Process ID:</strong> {exportResult.process_id}
              </Typography>
            </Box>
          )}

          <Typography variant="body2" color="text.secondary">
            The export is running in the background. You can monitor its
            progress and download the completed export by going to{" "}
            <strong>File/Projects → Exports...</strong>
          </Typography>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setExportSuccessDialogOpen(false)}>
            Close
          </Button>
          <Button
            variant="contained"
            onClick={() => {
              setExportSuccessDialogOpen(false);
              setExportsDialogOpen(true);
            }}
          >
            View Exports
          </Button>
        </DialogActions>
      </Dialog>

      <ProjectExportsDialog
        open={exportsDialogOpen}
        onClose={() => setExportsDialogOpen(false)}
      />
    </>
  );
}
