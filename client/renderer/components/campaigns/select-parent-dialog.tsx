"use client";

import { useState } from "react";
import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  LinearProgress,
  List,
  ListItemButton,
  ListItemIcon,
  ListItemText,
  TextField,
  Typography,
} from "@mui/material";
import { FolderOpen as FolderIcon } from "@mui/icons-material";
import { useApi } from "../../api";
import { useCampaignsApi } from "../../lib/campaigns-api";
import { Project } from "../../types/models";

interface SelectParentDialogProps {
  open: boolean;
  onClose: () => void;
  campaignId: number;
}

export function SelectParentDialog({
  open,
  onClose,
  campaignId,
}: SelectParentDialogProps) {
  const api = useApi();
  const campaignsApi = useCampaignsApi();
  const { data: projects } = api.get<Project[]>("projects");

  const [selectedProject, setSelectedProject] = useState<Project | null>(null);
  const [searchFilter, setSearchFilter] = useState("");
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const filteredProjects = projects?.filter((p) =>
    p.name.toLowerCase().includes(searchFilter.toLowerCase())
  );

  const handleSave = async () => {
    if (!selectedProject) return;

    setSaving(true);
    setError(null);
    try {
      await campaignsApi.setParentProject(campaignId, selectedProject.id);
      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : "Failed to set parent project");
    } finally {
      setSaving(false);
    }
  };

  const handleClose = () => {
    if (!saving) {
      setSelectedProject(null);
      setSearchFilter("");
      setError(null);
      onClose();
    }
  };

  return (
    <Dialog open={open} onClose={handleClose} maxWidth="sm" fullWidth>
      <DialogTitle>Select Parent Project</DialogTitle>
      <DialogContent>
        <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
          The parent project provides reference coordinates and FreeR flags that
          will be used for all member datasets.
        </Typography>

        <TextField
          fullWidth
          size="small"
          placeholder="Search projects..."
          value={searchFilter}
          onChange={(e) => setSearchFilter(e.target.value)}
          sx={{ mb: 2 }}
        />

        <List sx={{ maxHeight: 300, overflow: "auto" }}>
          {filteredProjects?.map((project) => (
            <ListItemButton
              key={project.id}
              selected={selectedProject?.id === project.id}
              onClick={() => setSelectedProject(project)}
            >
              <ListItemIcon>
                <FolderIcon
                  color={
                    selectedProject?.id === project.id ? "primary" : "inherit"
                  }
                />
              </ListItemIcon>
              <ListItemText
                primary={project.name}
                secondary={project.description || undefined}
              />
            </ListItemButton>
          ))}
          {filteredProjects?.length === 0 && (
            <Typography color="text.secondary" sx={{ p: 2, textAlign: "center" }}>
              No projects found
            </Typography>
          )}
        </List>

        {error && (
          <Typography color="error" sx={{ mt: 1 }}>
            {error}
          </Typography>
        )}

        {saving && <LinearProgress sx={{ mt: 2 }} />}
      </DialogContent>
      <DialogActions>
        <Button onClick={handleClose} disabled={saving}>
          Cancel
        </Button>
        <Button
          onClick={handleSave}
          variant="contained"
          disabled={!selectedProject || saving}
        >
          Set as Parent
        </Button>
      </DialogActions>
    </Dialog>
  );
}
