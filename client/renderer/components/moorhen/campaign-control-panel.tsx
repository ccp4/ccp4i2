"use client";

/**
 * Campaign Control Panel for Moorhen viewer.
 *
 * Provides campaign-specific controls:
 * - Binding site navigation
 * - Member project switching
 * - Copy view link functionality
 * - Push modified structures back to CCP4i2
 */

import React, { useState, useCallback, useMemo } from "react";
import {
  Box,
  Typography,
  Button,
  IconButton,
  List,
  ListItem,
  ListItemText,
  ListItemSecondaryAction,
  Select,
  MenuItem,
  FormControl,
  FormControlLabel,
  Checkbox,
  InputLabel,
  TextField,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Divider,
  Tooltip,
  Paper,
  Switch,
} from "@mui/material";
import {
  Place as PlaceIcon,
  Add as AddIcon,
  Delete as DeleteIcon,
  Edit as EditIcon,
  Home as HomeIcon,
  Upload as UploadIcon,
} from "@mui/icons-material";
import { moorhen } from "moorhen/types/moorhen";
import { CopyViewLinkButton } from "./copy-view-link-button";
import { PushToCCP4i2Panel } from "./push-to-ccp4i2-panel";
import {
  ProjectGroup,
  CampaignSite,
  MemberProjectWithSummary,
} from "../../types/campaigns";
import { Project } from "../../types/models";

interface CampaignControlPanelProps {
  campaign: ProjectGroup;
  sites: CampaignSite[];
  onGoToSite: (site: CampaignSite) => void;
  onSaveCurrentAsSite: (name: string) => Promise<void>;
  onUpdateSite: (index: number, name: string, updatePosition: boolean) => Promise<void>;
  onDeleteSite: (index: number) => Promise<void>;
  memberProjects: MemberProjectWithSummary[];
  selectedMemberProjectId: number | null;
  onSelectMemberProject: (projectId: number | null) => void;
  parentProject: Project | null | undefined;
  getViewUrl: () => string;
  /** Molecules loaded in Moorhen for push-to-CCP4i2 */
  molecules?: moorhen.Molecule[];
}

export const CampaignControlPanel: React.FC<CampaignControlPanelProps> = ({
  campaign,
  sites,
  onGoToSite,
  onSaveCurrentAsSite,
  onUpdateSite,
  onDeleteSite,
  memberProjects,
  selectedMemberProjectId,
  onSelectMemberProject,
  parentProject,
  getViewUrl,
  molecules,
}) => {
  const [showAddSiteDialog, setShowAddSiteDialog] = useState(false);
  const [newSiteName, setNewSiteName] = useState("");
  const [isSaving, setIsSaving] = useState(false);

  // Edit site dialog state
  const [editingSiteIndex, setEditingSiteIndex] = useState<number | null>(null);
  const [editSiteName, setEditSiteName] = useState("");
  const [updatePosition, setUpdatePosition] = useState(false);

  // Push to CCP4i2 dialog state
  const [showPushDialog, setShowPushDialog] = useState(false);
  const [selectedMolecule, setSelectedMolecule] = useState<moorhen.Molecule | null>(null);

  // Display options state
  const [showCBs, setShowCBs] = useState(false);

  // Get the currently viewed project (member project or parent)
  const currentProject = useMemo(() => {
    if (selectedMemberProjectId) {
      const memberProject = memberProjects.find((p) => p.id === selectedMemberProjectId);
      if (memberProject) {
        // Convert to Project type for PushToCCP4i2Panel
        return {
          id: memberProject.id,
          uuid: memberProject.uuid,
          name: memberProject.name,
        } as Project;
      }
    }
    return parentProject || undefined;
  }, [selectedMemberProjectId, memberProjects, parentProject]);

  const handleOpenPushDialog = useCallback((molecule: moorhen.Molecule) => {
    setSelectedMolecule(molecule);
    setShowPushDialog(true);
  }, []);

  const handleClosePushDialog = useCallback(() => {
    setShowPushDialog(false);
    setSelectedMolecule(null);
  }, []);

  const handleToggleCBs = useCallback(async (enabled: boolean) => {
    setShowCBs(enabled);
    if (!molecules) return;

    for (const mol of molecules) {
      try {
        if (enabled) {
          await mol.addRepresentation("CBs", "/*/*/*/*");
        } else {
          mol.clearBuffersOfStyle("CBs");
        }
      } catch (err) {
        console.error(`Failed to toggle CBs for molecule ${mol.name}:`, err);
      }
    }
  }, [molecules]);

  const handleSaveSite = useCallback(async () => {
    if (!newSiteName.trim()) return;
    setIsSaving(true);
    try {
      await onSaveCurrentAsSite(newSiteName.trim());
      setNewSiteName("");
      setShowAddSiteDialog(false);
    } finally {
      setIsSaving(false);
    }
  }, [newSiteName, onSaveCurrentAsSite]);

  const handleDeleteSite = useCallback(
    async (index: number) => {
      await onDeleteSite(index);
    },
    [onDeleteSite]
  );

  const handleOpenEditDialog = useCallback((index: number) => {
    setEditingSiteIndex(index);
    setEditSiteName(sites[index].name);
    setUpdatePosition(false);
  }, [sites]);

  const handleCloseEditDialog = useCallback(() => {
    setEditingSiteIndex(null);
    setEditSiteName("");
    setUpdatePosition(false);
  }, []);

  const handleSaveEdit = useCallback(async () => {
    if (editingSiteIndex === null || !editSiteName.trim()) return;
    setIsSaving(true);
    try {
      await onUpdateSite(editingSiteIndex, editSiteName.trim(), updatePosition);
      handleCloseEditDialog();
    } finally {
      setIsSaving(false);
    }
  }, [editingSiteIndex, editSiteName, updatePosition, onUpdateSite, handleCloseEditDialog]);

  return (
    <Box sx={{ p: 2, height: "100%", display: "flex", flexDirection: "column" }}>
      {/* Header */}
      <Typography variant="h6" sx={{ mb: 2, fontWeight: "bold" }}>
        {campaign.name}
      </Typography>

      {/* Copy View Link */}
      <Box sx={{ mb: 2 }}>
        <CopyViewLinkButton getViewUrl={getViewUrl} />
      </Box>

      {/* Display Options */}
      {molecules && molecules.length > 0 && (
        <Box sx={{ mb: 2 }}>
          <FormControlLabel
            control={
              <Switch
                checked={showCBs}
                onChange={(e) => handleToggleCBs(e.target.checked)}
                size="small"
              />
            }
            label={<Typography variant="body2">Show bonds (CBs)</Typography>}
          />
        </Box>
      )}

      <Divider sx={{ my: 2 }} />

      {/* Member Project Selector */}
      <FormControl fullWidth size="small" sx={{ mb: 2 }}>
        <InputLabel id="member-project-label">View Project</InputLabel>
        <Select
          labelId="member-project-label"
          value={selectedMemberProjectId ?? "parent"}
          label="View Project"
          onChange={(e) => {
            const value = e.target.value;
            onSelectMemberProject(value === "parent" ? null : (value as number));
          }}
        >
          <MenuItem value="parent">
            <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
              <HomeIcon fontSize="small" />
              <span>{parentProject?.name || "Parent Project"}</span>
            </Box>
          </MenuItem>
          <Divider />
          {memberProjects.map((project) => (
            <MenuItem key={project.id} value={project.id}>
              <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
                <Typography variant="body2">{project.name}</Typography>
                {project.job_summary && (
                  <Typography
                    variant="caption"
                    color="text.secondary"
                    sx={{ ml: "auto" }}
                  >
                    ({project.job_summary.finished}/{project.job_summary.total})
                  </Typography>
                )}
              </Box>
            </MenuItem>
          ))}
        </Select>
      </FormControl>

      <Divider sx={{ my: 2 }} />

      {/* Sites Section */}
      <Box sx={{ flex: 1, display: "flex", flexDirection: "column", minHeight: 0 }}>
        <Box
          sx={{
            display: "flex",
            alignItems: "center",
            justifyContent: "space-between",
            mb: 1,
          }}
        >
          <Typography variant="subtitle2" sx={{ fontWeight: "bold" }}>
            Binding Sites
          </Typography>
          <Tooltip title="Save current view as a site">
            <IconButton
              size="small"
              onClick={() => setShowAddSiteDialog(true)}
              color="primary"
            >
              <AddIcon fontSize="small" />
            </IconButton>
          </Tooltip>
        </Box>

        {sites.length === 0 ? (
          <Paper
            variant="outlined"
            sx={{
              p: 2,
              textAlign: "center",
              color: "text.secondary",
              bgcolor: "background.default",
            }}
          >
            <PlaceIcon sx={{ fontSize: 40, mb: 1, opacity: 0.5 }} />
            <Typography variant="body2">
              No sites saved yet.
            </Typography>
            <Typography variant="caption" display="block" sx={{ mt: 1 }}>
              Navigate to a binding site and click + to save it.
            </Typography>
          </Paper>
        ) : (
          <List dense sx={{ flex: 1, overflow: "auto" }}>
            {sites.map((site, index) => (
              <ListItem
                key={index}
                sx={{
                  bgcolor: "background.paper",
                  mb: 0.5,
                  borderRadius: 1,
                  border: "1px solid",
                  borderColor: "divider",
                }}
              >
                <ListItemText
                  primary={site.name}
                  primaryTypographyProps={{ variant: "body2" }}
                />
                <ListItemSecondaryAction>
                  <Tooltip title="Go to site">
                    <IconButton
                      edge="end"
                      size="small"
                      onClick={() => onGoToSite(site)}
                      color="primary"
                    >
                      <PlaceIcon fontSize="small" />
                    </IconButton>
                  </Tooltip>
                  <Tooltip title="Edit site">
                    <IconButton
                      edge="end"
                      size="small"
                      onClick={() => handleOpenEditDialog(index)}
                    >
                      <EditIcon fontSize="small" />
                    </IconButton>
                  </Tooltip>
                  <Tooltip title="Delete site">
                    <IconButton
                      edge="end"
                      size="small"
                      onClick={() => handleDeleteSite(index)}
                      color="error"
                    >
                      <DeleteIcon fontSize="small" />
                    </IconButton>
                  </Tooltip>
                </ListItemSecondaryAction>
              </ListItem>
            ))}
          </List>
        )}
      </Box>

      {/* Push to CCP4i2 Section */}
      {molecules && molecules.length > 0 && (
        <>
          <Divider sx={{ my: 2 }} />
          <Box>
            <Typography variant="subtitle2" sx={{ fontWeight: "bold", mb: 1 }}>
              Push to CCP4i2
            </Typography>
            <List dense>
              {molecules.map((mol) => (
                <ListItem
                  key={mol.molNo}
                  sx={{
                    bgcolor: "background.paper",
                    mb: 0.5,
                    borderRadius: 1,
                    border: "1px solid",
                    borderColor: "divider",
                  }}
                >
                  <ListItemText
                    primary={mol.name}
                    primaryTypographyProps={{ variant: "body2" }}
                  />
                  <ListItemSecondaryAction>
                    <Tooltip title="Push to CCP4i2 project">
                      <IconButton
                        edge="end"
                        size="small"
                        onClick={() => handleOpenPushDialog(mol)}
                        color="primary"
                      >
                        <UploadIcon fontSize="small" />
                      </IconButton>
                    </Tooltip>
                  </ListItemSecondaryAction>
                </ListItem>
              ))}
            </List>
          </Box>
        </>
      )}

      {/* Add Site Dialog */}
      <Dialog
        open={showAddSiteDialog}
        onClose={() => setShowAddSiteDialog(false)}
        maxWidth="xs"
        fullWidth
      >
        <DialogTitle>Save Current View as Site</DialogTitle>
        <DialogContent>
          <TextField
            autoFocus
            margin="dense"
            label="Site Name"
            fullWidth
            variant="outlined"
            value={newSiteName}
            onChange={(e) => setNewSiteName(e.target.value)}
            placeholder="e.g., Active Site, Binding Pocket"
            onKeyDown={(e) => {
              if (e.key === "Enter" && newSiteName.trim()) {
                handleSaveSite();
              }
            }}
          />
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setShowAddSiteDialog(false)}>Cancel</Button>
          <Button
            onClick={handleSaveSite}
            variant="contained"
            disabled={!newSiteName.trim() || isSaving}
          >
            {isSaving ? "Saving..." : "Save"}
          </Button>
        </DialogActions>
      </Dialog>

      {/* Edit Site Dialog */}
      <Dialog
        open={editingSiteIndex !== null}
        onClose={handleCloseEditDialog}
        maxWidth="xs"
        fullWidth
      >
        <DialogTitle>Edit Site</DialogTitle>
        <DialogContent>
          <TextField
            autoFocus
            margin="dense"
            label="Site Name"
            fullWidth
            variant="outlined"
            value={editSiteName}
            onChange={(e) => setEditSiteName(e.target.value)}
            placeholder="e.g., Active Site, Binding Pocket"
            onKeyDown={(e) => {
              if (e.key === "Enter" && editSiteName.trim()) {
                handleSaveEdit();
              }
            }}
          />
          <FormControlLabel
            control={
              <Checkbox
                checked={updatePosition}
                onChange={(e) => setUpdatePosition(e.target.checked)}
              />
            }
            label="Update position to current view"
            sx={{ mt: 1 }}
          />
        </DialogContent>
        <DialogActions>
          <Button onClick={handleCloseEditDialog}>Cancel</Button>
          <Button
            onClick={handleSaveEdit}
            variant="contained"
            disabled={!editSiteName.trim() || isSaving}
          >
            {isSaving ? "Saving..." : "Save"}
          </Button>
        </DialogActions>
      </Dialog>

      {/* Push to CCP4i2 Dialog */}
      <Dialog
        open={showPushDialog}
        onClose={handleClosePushDialog}
        maxWidth="sm"
        fullWidth
      >
        <DialogContent>
          {selectedMolecule && (
            <PushToCCP4i2Panel
              project={currentProject}
              molNo={selectedMolecule.molNo}
              item={selectedMolecule}
              onClose={handleClosePushDialog}
            />
          )}
        </DialogContent>
      </Dialog>
    </Box>
  );
};
