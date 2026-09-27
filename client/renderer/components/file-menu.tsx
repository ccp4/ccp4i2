"use client";
import { useEffect, useState } from "react";
import {
  Button,
  Divider,
  Menu,
  MenuItem,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Typography,
  Box,
} from "@mui/material";
import { Add, Download, Edit, FolderOpen, Menu as MenuIcon, OpenInNew, Upload, PowerSettingsNew } from "@mui/icons-material";
import { useApi } from "../api";
import { Project } from "../types/models";
import { apiPost } from "../api-fetch";
import { ProjectExportsDialog } from "./project-exports";
import { CCP4i2MenuItem } from "./menu-item";
import { useProjectScope } from "../lib/project-scope";
import { isElectron } from "../utils/platform";
import { shortDate } from "../pipes";

interface ExportResult {
  status: string;
  export_file_name: string;
  log_file_name: string;
  process_id: number;
}

export default function FileMenu() {
  const api = useApi();
  const projectId = useProjectScope();
  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
  const open = Boolean(anchorEl);
  // Ten recent projects, fetched when the menu is opened rather than on every
  // page that mounts this bar: the list is the whole catalogue, 2.45 MB at
  // DDU's 7,204 projects, for ten menu items.
  const { data: projects, mutate: mutateProjects } = api.get<Project[]>(
    open ? "projects" : null
  );
  // Quit is Electron's to offer, and the check is deferred to an effect: this
  // menu is rendered by the app bar on every route, so reading `window` during
  // render broke server rendering of all of them.
  const [canQuit, setCanQuit] = useState(false);
  useEffect(() => setCanQuit(isElectron()), []);
  const [exportsDialogOpen, setExportsDialogOpen] = useState(false);
  const [exportSuccessDialogOpen, setExportSuccessDialogOpen] = useState(false);
  const [exportResult, setExportResult] = useState<ExportResult | null>(null);

  const handleClick = (event: React.MouseEvent<HTMLButtonElement>) => {
    setAnchorEl(event.currentTarget);
  };
  const handleClose = () => {
    setAnchorEl(null);
  };

  const handleNewWindow = () => {
    const newWindow = window.open("/ccp4i2");
    setAnchorEl(null);
    // Check if the window was successfully opened
    if (!newWindow) {
      console.warn("Failed to open new window. It might be blocked by a popup blocker.");
    }
  };

  const handleExportProject = () => {
    apiPost(`projects/${projectId}/export`, {})
      .then((result: ExportResult) => {
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
    if (!newWindow) {
      console.warn("Failed to open new window. It might be blocked by a popup blocker.");
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
        <CCP4i2MenuItem
          text="Manage/open projects"
          icon={MenuIcon}
          onClick={() => {
            handleClose();
            window.open("/ccp4i2");
          }}
        />
        <CCP4i2MenuItem
          text="New project"
          icon={Add}
          onClick={() => {
            handleClose();
            window.open("/ccp4i2/new-project");
          }}
        />
        {projectId && (
          <CCP4i2MenuItem
            text="Edit project..."
            icon={Edit}
            onClick={() => {
              handleClose();
              window.open(`/ccp4i2/edit-project/${projectId}`);
            }}
          />
        )}
        {projectId && (
          <CCP4i2MenuItem
            text="Export Project"
            icon={Download}
            onClick={handleExportProject}
          />
        )}
        <CCP4i2MenuItem
          text="Import Project"
          icon={Upload}
          onClick={handleImportProject}
        />
        <CCP4i2MenuItem
          text="Exports..."
          icon={Download}
          onClick={handleOpenExportsDialog}
        />
        <Divider />
        {Array.isArray(projects) &&
          // The server returns the list ordered by -last_access, so the ten
          // most recent are the first ten. The sort this replaces ran in
          // place, which reordered SWR's cached array under every other
          // consumer of it.
          projects
            .slice(0, 10)
            .map((project: Project) => (
              <CCP4i2MenuItem
                key={project.id}
                text={project.name}
                icon={FolderOpen}
                secondary={shortDate(new Date(project.last_access))}
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
                  window.open(`/ccp4i2/project/${project.id}`);
                }}
              />
            ))}
        <Divider />
        <CCP4i2MenuItem
          text="New Window"
          icon={OpenInNew}
          onClick={handleNewWindow}
        />
        {canQuit && (
          <CCP4i2MenuItem
            text="Quit"
            icon={PowerSettingsNew}
            onClick={() => {
              window.electronAPI?.sendMessage("quit-app", {});
            }}
          />
        )}
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
