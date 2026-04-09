"use client";

import { useState } from "react";
import { useRouter } from "next/navigation";
import {
  Box,
  Button,
  Container,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  IconButton,
  LinearProgress,
  Paper,
  Skeleton,
  Stack,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  TextField,
  Tooltip,
  Typography,
} from "@mui/material";
import {
  Add as AddIcon,
  Delete as DeleteIcon,
  Refresh as RefreshIcon,
  Science as ScienceIcon,
} from "@mui/icons-material";
import CCP4i2TopBar from "../../../components/ccp4i2-topbar";
import { useCampaignsApi } from "../../../lib/campaigns-api";
import { ProjectGroup } from "../../../types/campaigns";

export default function CampaignsPage() {
  const router = useRouter();
  const campaignsApi = useCampaignsApi();
  const { data: campaigns, isLoading, mutate } = campaignsApi.useCampaigns(5000);

  const [createDialogOpen, setCreateDialogOpen] = useState(false);
  const [newCampaignName, setNewCampaignName] = useState("");
  const [creating, setCreating] = useState(false);
  const [createError, setCreateError] = useState<string | null>(null);
  const [deleteDialogOpen, setDeleteDialogOpen] = useState(false);
  const [campaignToDelete, setCampaignToDelete] = useState<ProjectGroup | null>(
    null
  );
  const [deleting, setDeleting] = useState(false);

  const handleCreateCampaign = async () => {
    if (!newCampaignName.trim()) return;

    setCreating(true);
    setCreateError(null);
    try {
      // Create the campaign with auto-created parent project
      const newCampaign = await campaignsApi.createCampaignWithParent(
        newCampaignName.trim()
      );
      setCreateDialogOpen(false);
      setNewCampaignName("");
      // Navigate to the new campaign
      router.push(`/ccp4i2/campaigns/${newCampaign.id}`);
    } catch (error) {
      console.error("Failed to create campaign:", error);
      // Extract error message from API response
      const message = error instanceof Error ? error.message : "Failed to create campaign";
      setCreateError(message);
    } finally {
      setCreating(false);
    }
  };

  const handleCloseCreateDialog = () => {
    setCreateDialogOpen(false);
    setNewCampaignName("");
    setCreateError(null);
  };

  const handleDeleteCampaign = async () => {
    if (!campaignToDelete) return;

    setDeleting(true);
    try {
      await campaignsApi.deleteCampaign(campaignToDelete.id);
      setDeleteDialogOpen(false);
      setCampaignToDelete(null);
      mutate();
    } catch (error) {
      console.error("Failed to delete campaign:", error);
    } finally {
      setDeleting(false);
    }
  };

  const handleRowClick = (campaign: ProjectGroup) => {
    router.push(`/ccp4i2/campaigns/${campaign.id}`);
  };

  return (
    <>
      <CCP4i2TopBar title="Fragment Screening Campaigns" showBackButton />
      <Container sx={{ my: 3 }}>
        <Stack spacing={3}>
          {/* Toolbar */}
          <Stack direction="row" spacing={2} alignItems="center">
            <Button
              variant="contained"
              startIcon={<AddIcon />}
              onClick={() => setCreateDialogOpen(true)}
            >
              New Campaign
            </Button>
            <Tooltip title="Refresh">
              <IconButton onClick={() => mutate()}>
                <RefreshIcon />
              </IconButton>
            </Tooltip>
          </Stack>

          {/* Campaigns Table */}
          {isLoading ? (
            <Skeleton variant="rectangular" width="100%" height={400} />
          ) : campaigns && campaigns.length > 0 ? (
            <TableContainer component={Paper}>
              <Table>
                <TableHead>
                  <TableRow>
                    <TableCell>Campaign Name</TableCell>
                    <TableCell align="center">Type</TableCell>
                    <TableCell align="right">Actions</TableCell>
                  </TableRow>
                </TableHead>
                <TableBody>
                  {campaigns.map((campaign) => (
                    <TableRow
                      key={campaign.id}
                      hover
                      sx={{ cursor: "pointer" }}
                      onClick={() => handleRowClick(campaign)}
                    >
                      <TableCell>
                        <Stack direction="row" spacing={1} alignItems="center">
                          <ScienceIcon color="primary" />
                          <Typography>{campaign.name}</Typography>
                        </Stack>
                      </TableCell>
                      <TableCell align="center">
                        <Typography variant="caption" color="text.secondary">
                          {campaign.type === "fragment_set"
                            ? "Fragment Set"
                            : "General Set"}
                        </Typography>
                      </TableCell>
                      <TableCell align="right">
                        <Tooltip title="Delete campaign">
                          <IconButton
                            size="small"
                            onClick={(e) => {
                              e.stopPropagation();
                              setCampaignToDelete(campaign);
                              setDeleteDialogOpen(true);
                            }}
                          >
                            <DeleteIcon />
                          </IconButton>
                        </Tooltip>
                      </TableCell>
                    </TableRow>
                  ))}
                </TableBody>
              </Table>
            </TableContainer>
          ) : (
            <Paper sx={{ p: 6, textAlign: "center" }}>
              <Stack spacing={2} alignItems="center">
                <ScienceIcon sx={{ fontSize: 64, color: "text.disabled" }} />
                <Typography variant="h6" color="text.secondary">
                  No campaigns yet
                </Typography>
                <Typography color="text.secondary">
                  Create a new fragment screening campaign to get started.
                </Typography>
                <Button
                  variant="outlined"
                  startIcon={<AddIcon />}
                  onClick={() => setCreateDialogOpen(true)}
                >
                  Create Campaign
                </Button>
              </Stack>
            </Paper>
          )}
        </Stack>
      </Container>

      {/* Create Campaign Dialog */}
      <Dialog
        open={createDialogOpen}
        onClose={handleCloseCreateDialog}
        maxWidth="sm"
        fullWidth
      >
        <DialogTitle>Create New Campaign</DialogTitle>
        <DialogContent>
          <TextField
            autoFocus
            margin="dense"
            label="Campaign Name"
            fullWidth
            variant="outlined"
            value={newCampaignName}
            onChange={(e) => {
              setNewCampaignName(e.target.value);
              setCreateError(null); // Clear error when typing
            }}
            onKeyDown={(e) => {
              if (e.key === "Enter") {
                handleCreateCampaign();
              }
            }}
            error={!!createError}
            helperText={createError || "A parent project will be automatically created with this name"}
          />
          {creating && <LinearProgress sx={{ mt: 2 }} />}
        </DialogContent>
        <DialogActions>
          <Button onClick={handleCloseCreateDialog} disabled={creating}>
            Cancel
          </Button>
          <Button
            onClick={handleCreateCampaign}
            variant="contained"
            disabled={creating || !newCampaignName.trim()}
          >
            Create
          </Button>
        </DialogActions>
      </Dialog>

      {/* Delete Confirmation Dialog */}
      <Dialog
        open={deleteDialogOpen}
        onClose={() => setDeleteDialogOpen(false)}
      >
        <DialogTitle>Delete Campaign?</DialogTitle>
        <DialogContent>
          <Typography>
            Are you sure you want to delete the campaign &quot;
            {campaignToDelete?.name}&quot;?
          </Typography>
          <Typography variant="body2" color="text.secondary" sx={{ mt: 1 }}>
            This will remove the campaign grouping but will not delete the
            individual projects.
          </Typography>
          {deleting && <LinearProgress sx={{ mt: 2 }} />}
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setDeleteDialogOpen(false)} disabled={deleting}>
            Cancel
          </Button>
          <Button
            onClick={handleDeleteCampaign}
            color="error"
            variant="contained"
            disabled={deleting}
          >
            Delete
          </Button>
        </DialogActions>
      </Dialog>
    </>
  );
}
