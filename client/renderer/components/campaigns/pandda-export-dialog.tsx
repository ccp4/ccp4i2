"use client";

import { useState, useEffect } from "react";
import {
  Alert,
  Box,
  Button,
  CircularProgress,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  LinearProgress,
  Link,
  List,
  ListItem,
  ListItemIcon,
  ListItemText,
  Stack,
  Typography,
} from "@mui/material";
import {
  CheckCircle as CheckIcon,
  Error as ErrorIcon,
  Download as DownloadIcon,
} from "@mui/icons-material";
import { useRouter } from "next/navigation";
import useSWR from "swr";
import { apiFetch } from "../../api-fetch";

interface PanddaDataset {
  project_name: string;
  project_id: number;
  dimple_job_id: number | null;
  acedrg_job_id: number | null;
  has_dimple: boolean;
  has_acedrg: boolean;
}

interface PanddaDataResponse {
  datasets: PanddaDataset[];
  total_ready: number;
  total_with_dict: number;
}

interface PanddaExportDialogProps {
  open: boolean;
  onClose: () => void;
  campaignId: number;
  campaignName: string;
  /** The campaign's parent project, where a PanDDA run is configured. */
  parentProject?: { id: number; name: string };
}

/**
 * Dialog for preparing and exporting PANDDA-ready data.
 *
 * Shows which datasets have finished Dimple jobs and optional
 * dictionary files, then allows downloading a ZIP file ready
 * for PANDDA analysis.
 */
export function PanddaExportDialog({
  open,
  onClose,
  campaignId,
  campaignName,
  parentProject,
}: PanddaExportDialogProps) {
  const router = useRouter();
  const [exporting, setExporting] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Fetch PANDDA data when dialog opens
  const { data, isLoading, error: fetchError } = useSWR<PanddaDataResponse>(
    open ? `projectgroups/${campaignId}/pandda_data/` : null,
    async (url: string) => {
      const response = await apiFetch(url);
      return response.json();
    }
  );

  // Reset state when dialog closes
  useEffect(() => {
    if (!open) {
      setError(null);
      setExporting(false);
    }
  }, [open]);

  const handleExport = async () => {
    setExporting(true);
    setError(null);

    try {
      const response = await apiFetch(
        `projectgroups/${campaignId}/export_pandda/`,
        { method: "POST" },
        { timeout: 10 * 60 * 1000 }
      );

      // Download the ZIP file
      const blob = await response.blob();
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = `pandda_${campaignName}.zip`;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(url);

      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : "Export failed");
    } finally {
      setExporting(false);
    }
  };

  return (
    <Dialog open={open} onClose={onClose} maxWidth="md" fullWidth>
      <DialogTitle>
        <Stack direction="row" spacing={1} alignItems="center">
          <DownloadIcon color="primary" />
          <Typography variant="h6">Export for PanDDA</Typography>
        </Stack>
      </DialogTitle>
      <DialogContent>
        <Stack spacing={2}>
          <Typography variant="body2" color="text.secondary">
            Export final.pdb, final.mtz, and dictionary files from finished
            Dimple jobs, to run PanDDA (Pan-Dataset Density Analysis)
            elsewhere.
          </Typography>

          {/* Information, not a shortcut: a PanDDA run is hours of compute
              that later results key off, so it stays a job somebody
              configures and checks before running. */}
          <Alert severity="info">
            To run PanDDA in CCP4i2 instead, configure a{" "}
            <strong>PanDDA over a set of datasets</strong> job in the
            campaign&apos;s parent project
            {parentProject ? (
              <>
                {" "}(
                <Link
                  component="button"
                  variant="body2"
                  sx={{ verticalAlign: "baseline" }}
                  onClick={() => {
                    onClose();
                    router.push(`/ccp4i2/project/${parentProject.id}`);
                  }}
                >
                  {parentProject.name}
                </Link>
                )
              </>
            ) : null}{" "}
            and use <strong>Fill from campaign</strong> to list the datasets.
          </Alert>

          {isLoading && (
            <Box sx={{ display: "flex", justifyContent: "center", py: 4 }}>
              <CircularProgress />
            </Box>
          )}

          {fetchError && (
            <Alert severity="error">Failed to load dataset information</Alert>
          )}

          {error && <Alert severity="error">{error}</Alert>}

          {data && (
            <>
              <Stack direction="row" spacing={4}>
                <Box>
                  <Typography variant="h4" color="success.main">
                    {data.total_ready}
                  </Typography>
                  <Typography variant="caption" color="text.secondary">
                    Datasets ready (with Dimple)
                  </Typography>
                </Box>
                <Box>
                  <Typography variant="h4" color="primary.main">
                    {data.total_with_dict}
                  </Typography>
                  <Typography variant="caption" color="text.secondary">
                    With dictionary files
                  </Typography>
                </Box>
              </Stack>

              {data.datasets.length === 0 ? (
                <Alert severity="warning">
                  No datasets have finished Dimple jobs yet. Run i2Dimple on
                  your member projects first.
                </Alert>
              ) : (
                <Box sx={{ maxHeight: 300, overflow: "auto" }}>
                  <List dense>
                    {data.datasets.map((dataset) => (
                      <ListItem key={dataset.project_id}>
                        <ListItemIcon>
                          {dataset.has_dimple ? (
                            <CheckIcon color="success" fontSize="small" />
                          ) : (
                            <ErrorIcon color="error" fontSize="small" />
                          )}
                        </ListItemIcon>
                        <ListItemText
                          primary={dataset.project_name}
                          secondary={
                            <Stack direction="row" spacing={1} component="span">
                              <Typography
                                variant="caption"
                                component="span"
                                color={
                                  dataset.has_dimple
                                    ? "success.main"
                                    : "text.secondary"
                                }
                              >
                                Dimple: {dataset.has_dimple ? "Yes" : "No"}
                              </Typography>
                              <Typography
                                variant="caption"
                                component="span"
                                color={
                                  dataset.has_acedrg
                                    ? "primary.main"
                                    : "text.secondary"
                                }
                              >
                                Dict: {dataset.has_acedrg ? "Yes" : "No"}
                              </Typography>
                            </Stack>
                          }
                        />
                      </ListItem>
                    ))}
                  </List>
                </Box>
              )}
            </>
          )}

          {exporting && <LinearProgress />}
        </Stack>
      </DialogContent>
      <DialogActions>
        <Button onClick={onClose} disabled={exporting}>
          Cancel
        </Button>
        <Button
          variant="contained"
          onClick={handleExport}
          disabled={!data || data.total_ready === 0 || exporting}
          startIcon={exporting ? <CircularProgress size={16} /> : <DownloadIcon />}
        >
          {exporting ? "Exporting..." : `Export ${data?.total_ready || 0} Datasets`}
        </Button>
      </DialogActions>
    </Dialog>
  );
}
