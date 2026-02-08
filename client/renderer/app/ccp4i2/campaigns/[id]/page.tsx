"use client";

import { use, useState } from "react";
import { useRouter } from "next/navigation";
import {
  Alert,
  Box,
  Button,
  CircularProgress,
  Collapse,
  Container,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  Divider,
  FormControlLabel,
  Grid,
  IconButton,
  Paper,
  Skeleton,
  Stack,
  Switch,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableRow,
  TextField,
  Tooltip,
  Typography,
} from "@mui/material";
import {
  Add as AddIcon,
  Download as DownloadIcon,
  Edit as EditIcon,
  ExpandMore as ExpandMoreIcon,
  ExpandLess as ExpandLessIcon,
  Preview as PreviewIcon,
  Refresh as RefreshIcon,
  Upload as UploadIcon,
  Science as ScienceIcon,
  FolderOpen as FolderIcon,
} from "@mui/icons-material";
import CCP4i2TopBar from "../../../../components/ccp4i2-topbar";
import {
  SelectParentDialog,
  BatchImportDialog,
  CoordsImportDialog,
  FreeRImportDialog,
  PanddaExportDialog,
  VirtualizedMemberProjectsTable,
} from "../../../../components/campaigns";
import { useCampaignsApi, useSmilesLookup } from "../../../../lib/campaigns-api";
import { doDownload } from "../../../../api";
import { apiDelete } from "../../../../api-fetch";
import {
  MemberProjectWithSummary,
  parseDatasetFilename,
} from "../../../../types/campaigns";
import { File as CCP4File } from "../../../../types/models";
import { apiPatch } from "../../../../api-fetch";

/**
 * Get a display label for a file.
 * Priority: annotation > "Job X: PARAM_NAME" > filename
 */
function getFileDisplayLabel(file: CCP4File): string {
  if (file.annotation && file.annotation.trim()) {
    return file.annotation;
  }
  // Fallback to job number and param name if available
  if (file.job && file.job_param_name) {
    return `Job ${file.job}: ${file.job_param_name}`;
  }
  return file.name;
}

interface CampaignDetailPageProps {
  params: Promise<{ id: string }>;
}

export default function CampaignDetailPage({ params }: CampaignDetailPageProps) {
  const { id } = use(params);
  const campaignId = parseInt(id);
  const router = useRouter();
  const campaignsApi = useCampaignsApi();

  // Data fetching
  const { data: campaign, isLoading: campaignLoading } =
    campaignsApi.useCampaign(campaignId, 5000);
  const { data: parentProject } =
    campaignsApi.useParentProject(campaignId, 5000);
  const { data: parentFiles, mutate: mutateParentFiles } =
    campaignsApi.useParentFiles(campaignId, 5000);
  const {
    data: memberProjects,
    isLoading: membersLoading,
    mutate: mutateMemberProjects,
  } = campaignsApi.useMemberProjects(campaignId, 5000);

  // Extract reg_ids from member project names for SMILES lookup
  const regIds =
    memberProjects
      ?.map((p) => {
        const parsed = parseDatasetFilename(p.name);
        return parsed.nclId ? parseInt(parsed.nclId) : null;
      })
      .filter((id): id is number => id !== null && id !== 0) || [];
  const { smilesMap } = useSmilesLookup(regIds);

  // Dialog states
  const [showSelectParent, setShowSelectParent] = useState(false);
  const [showImportCoords, setShowImportCoords] = useState(false);
  const [showImportFreeR, setShowImportFreeR] = useState(false);
  const [showBatchImport, setShowBatchImport] = useState(false);
  const [showPanddaExport, setShowPanddaExport] = useState(false);
  const [showSubJobs, setShowSubJobs] = useState(false);
  const [deleteProject, setDeleteProject] = useState<MemberProjectWithSummary | null>(null);

  // Campaign info collapse state - expand by default when configuration is needed
  const needsConfiguration = !parentProject || !parentFiles?.coordinates?.length || !parentFiles?.freer?.length;
  const [infoExpanded, setInfoExpanded] = useState(false);
  const showExpanded = infoExpanded || needsConfiguration;

  // File annotation editing state
  const [editingFile, setEditingFile] = useState<CCP4File | null>(null);
  const [editAnnotation, setEditAnnotation] = useState("");
  const [savingAnnotation, setSavingAnnotation] = useState(false);

  const handleEditAnnotation = (file: CCP4File) => {
    setEditingFile(file);
    setEditAnnotation(file.annotation || "");
  };

  const handleSaveAnnotation = async () => {
    if (!editingFile) return;
    setSavingAnnotation(true);
    try {
      await apiPatch(`files/${editingFile.id}/`, { annotation: editAnnotation });
      mutateParentFiles();
      setEditingFile(null);
    } catch (err) {
      console.error("Failed to save annotation:", err);
    } finally {
      setSavingAnnotation(false);
    }
  };

  const handleRefresh = () => {
    mutateParentFiles();
    mutateMemberProjects();
  };

  const handleProjectClick = (project: MemberProjectWithSummary) => {
    router.push(`/ccp4i2/project/${project.id}`);
  };

  const handleDeleteProject = async () => {
    if (!deleteProject) return;
    try {
      await apiDelete(`projects/${deleteProject.id}/`);
      // Also remove from campaign
      await campaignsApi.removeMember(campaignId, deleteProject.id);
      mutateMemberProjects();
    } catch (err) {
      console.error("Failed to delete project:", err);
    }
    setDeleteProject(null);
  };

  const handleDownloadFile = (file: CCP4File) => {
    doDownload(
      `/api/proxy/ccp4i2/files/${file.uuid}/download_by_uuid/`,
      file.name
    );
  };

  if (campaignLoading) {
    return (
      <>
        <CCP4i2TopBar title="Loading..." showBackButton />
        <Container sx={{ my: 3 }}>
          <Skeleton variant="rectangular" width="100%" height={600} />
        </Container>
      </>
    );
  }

  if (!campaign) {
    return (
      <>
        <CCP4i2TopBar title="Campaign Not Found" showBackButton />
        <Container sx={{ my: 3 }}>
          <Alert severity="error">Campaign not found</Alert>
        </Container>
      </>
    );
  }

  return (
    <>
      <CCP4i2TopBar title={campaign.name} showBackButton />
      <Container sx={{ my: 3 }}>
        <Stack spacing={3}>
          {/* Header with actions and collapsible campaign info */}
          <Paper sx={{ p: 2 }}>
            <Stack direction="row" spacing={2} alignItems="center">
              <ScienceIcon color="primary" fontSize="large" />
              <Box sx={{ flexGrow: 1, minWidth: 0 }}>
                <Typography variant="h5">
                  {campaign.name}
                </Typography>
                {/* Collapsed summary */}
                {!showExpanded && (
                  <Stack direction="row" spacing={2} sx={{ mt: 0.5 }}>
                    <Typography variant="body2" color="text.secondary">
                      <strong>Parent:</strong>{" "}
                      {parentProject?.name || "None"}
                    </Typography>
                    <Typography variant="body2" color="text.secondary">
                      <strong>Coords:</strong>{" "}
                      {parentFiles?.coordinates?.length || 0} file{parentFiles?.coordinates?.length !== 1 ? "s" : ""}
                    </Typography>
                    <Typography variant="body2" color="text.secondary">
                      <strong>FreeR:</strong>{" "}
                      {parentFiles?.freer?.length || 0} file{parentFiles?.freer?.length !== 1 ? "s" : ""}
                    </Typography>
                  </Stack>
                )}
              </Box>
              <Tooltip title="Refresh">
                <IconButton onClick={handleRefresh}>
                  <RefreshIcon />
                </IconButton>
              </Tooltip>
              <Tooltip title={showExpanded ? "Hide details" : "Show details"}>
                <IconButton onClick={() => setInfoExpanded(!infoExpanded)}>
                  {showExpanded ? <ExpandLessIcon /> : <ExpandMoreIcon />}
                </IconButton>
              </Tooltip>
            </Stack>

            <Collapse in={showExpanded} timeout={250}>
              <Divider sx={{ my: 2 }} />

              {/* Parent Project Section */}
              <Stack spacing={2}>
                <Box>
                  <Typography variant="h6" gutterBottom>
                    Parent Project
                  </Typography>
                  {parentProject ? (
                    <Stack direction="row" spacing={1} alignItems="center">
                      <FolderIcon color="primary" />
                      <Typography
                        sx={{ cursor: "pointer", "&:hover": { textDecoration: "underline" } }}
                        onClick={() => router.push(`/ccp4i2/project/${parentProject.id}`)}
                      >
                        {parentProject.name}
                      </Typography>
                    </Stack>
                  ) : (
                    <Alert severity="warning">
                      Parent project not found. This campaign may have been created
                      with an older version. Use the button below to associate an
                      existing project.
                      <Button
                        size="small"
                        sx={{ ml: 2 }}
                        onClick={() => setShowSelectParent(true)}
                      >
                        Select Parent
                      </Button>
                    </Alert>
                  )}
                </Box>

                {/* Reference Files Section */}
                <Grid container spacing={2}>
                  {/* Coordinates */}
                  <Grid item xs={12} md={6}>
                    <Paper variant="outlined" sx={{ p: 2, height: "100%" }}>
                      <Stack
                        direction="row"
                        justifyContent="space-between"
                        alignItems="center"
                        mb={2}
                      >
                        <Typography variant="h6">Starting Coordinates</Typography>
                        {parentProject && (
                          <Button
                            size="small"
                            startIcon={<UploadIcon />}
                            onClick={() => setShowImportCoords(true)}
                            disabled={!parentProject}
                          >
                            Import
                          </Button>
                        )}
                      </Stack>
                      {parentFiles?.coordinates && parentFiles.coordinates.length > 0 ? (
                        <TableContainer>
                          <Table size="small">
                            <TableBody>
                              {parentFiles.coordinates.map((file) => (
                                <TableRow key={file.id}>
                                  <TableCell>
                                    <Tooltip title={file.name}>
                                      <Typography variant="body2" noWrap>
                                        {getFileDisplayLabel(file)}
                                      </Typography>
                                    </Tooltip>
                                  </TableCell>
                                  <TableCell align="right" sx={{ whiteSpace: "nowrap" }}>
                                    <Tooltip title="Edit annotation">
                                      <IconButton
                                        size="small"
                                        onClick={() => handleEditAnnotation(file)}
                                      >
                                        <EditIcon fontSize="small" />
                                      </IconButton>
                                    </Tooltip>
                                    <Tooltip title="Download">
                                      <IconButton
                                        size="small"
                                        onClick={() => handleDownloadFile(file)}
                                      >
                                        <DownloadIcon fontSize="small" />
                                      </IconButton>
                                    </Tooltip>
                                  </TableCell>
                                </TableRow>
                              ))}
                            </TableBody>
                          </Table>
                        </TableContainer>
                      ) : (
                        <Typography color="text.secondary" variant="body2">
                          No coordinates imported yet
                        </Typography>
                      )}
                    </Paper>
                  </Grid>

                  {/* FreeR */}
                  <Grid item xs={12} md={6}>
                    <Paper variant="outlined" sx={{ p: 2, height: "100%" }}>
                      <Stack
                        direction="row"
                        justifyContent="space-between"
                        alignItems="center"
                        mb={2}
                      >
                        <Typography variant="h6">Collective FreeR Set</Typography>
                        {parentProject && (
                          <Button
                            size="small"
                            startIcon={<UploadIcon />}
                            onClick={() => setShowImportFreeR(true)}
                            disabled={!parentProject}
                          >
                            Import
                          </Button>
                        )}
                      </Stack>
                      {parentFiles?.freer && parentFiles.freer.length > 0 ? (
                        <TableContainer>
                          <Table size="small">
                            <TableBody>
                              {parentFiles.freer.map((file) => (
                                <TableRow key={file.id}>
                                  <TableCell>
                                    <Tooltip title={file.name}>
                                      <Typography variant="body2" noWrap>
                                        {getFileDisplayLabel(file)}
                                      </Typography>
                                    </Tooltip>
                                  </TableCell>
                                  <TableCell align="right" sx={{ whiteSpace: "nowrap" }}>
                                    <Tooltip title="Edit annotation">
                                      <IconButton
                                        size="small"
                                        onClick={() => handleEditAnnotation(file)}
                                      >
                                        <EditIcon fontSize="small" />
                                      </IconButton>
                                    </Tooltip>
                                    <Tooltip title="Download">
                                      <IconButton
                                        size="small"
                                        onClick={() => handleDownloadFile(file)}
                                      >
                                        <DownloadIcon fontSize="small" />
                                      </IconButton>
                                    </Tooltip>
                                  </TableCell>
                                </TableRow>
                              ))}
                            </TableBody>
                          </Table>
                        </TableContainer>
                      ) : (
                        <Typography color="text.secondary" variant="body2">
                          No FreeR set imported yet
                        </Typography>
                      )}
                    </Paper>
                  </Grid>
                </Grid>
              </Stack>
            </Collapse>
          </Paper>

          {/* Member Projects Section */}
          <Paper sx={{ p: 2 }}>
            <Stack
              direction="row"
              justifyContent="space-between"
              alignItems="center"
              mb={2}
            >
              <Stack direction="row" spacing={2} alignItems="center">
                <Typography variant="h6">
                  Child Projects ({memberProjects?.length || 0})
                </Typography>
                <FormControlLabel
                  control={
                    <Switch
                      checked={showSubJobs}
                      onChange={(e) => setShowSubJobs(e.target.checked)}
                      size="small"
                    />
                  }
                  label="Show subjobs"
                />
              </Stack>
              <Stack direction="row" spacing={1}>
                <Button
                  variant="contained"
                  startIcon={<AddIcon />}
                  onClick={() => setShowBatchImport(true)}
                  disabled={
                    !parentFiles?.coordinates?.length ||
                    !parentFiles?.freer?.length
                  }
                >
                  Batch Import
                </Button>
                <Button
                  variant="outlined"
                  startIcon={<ScienceIcon />}
                  onClick={() => setShowPanddaExport(true)}
                  disabled={!memberProjects?.length}
                >
                  PANDDA
                </Button>
                <Button
                  variant="outlined"
                  startIcon={<PreviewIcon />}
                  onClick={() => {
                    // Open Summary View - Moorhen with all finished refmac jobs overlaid
                    // For now, navigate to parent project's Moorhen view if there's a reference structure
                    if (parentFiles?.coordinates?.[0]) {
                      // Find first member with finished refmac
                      const memberWithRefmac = memberProjects?.find(p =>
                        p.jobs?.some(j => j.task_name === "refmac" && j.status === 6)
                      );
                      if (memberWithRefmac) {
                        const refmacJob = memberWithRefmac.jobs?.find(
                          j => j.task_name === "refmac" && j.status === 6
                        );
                        if (refmacJob) {
                          window.open(`/ccp4i2/moorhen-page/job-by-id/${refmacJob.id}`, '_blank');
                        }
                      }
                    }
                  }}
                  disabled={!memberProjects?.some(p =>
                    p.jobs?.some(j => j.task_name === "refmac" && j.status === 6)
                  )}
                >
                  Summary View
                </Button>
              </Stack>
            </Stack>

            {!parentFiles?.coordinates?.length || !parentFiles?.freer?.length ? (
              <Alert severity="warning" sx={{ mb: 2 }}>
                Import coordinates and FreeR set before adding member projects.
              </Alert>
            ) : null}

            {membersLoading ? (
              <Skeleton variant="rectangular" width="100%" height={300} />
            ) : (
              <VirtualizedMemberProjectsTable
                projects={memberProjects || []}
                smilesMap={smilesMap}
                showSubJobs={showSubJobs}
                latestCoordsFileId={parentFiles?.coordinates?.[parentFiles.coordinates.length - 1]?.id}
                campaignId={campaignId}
                onRefresh={mutateMemberProjects}
                onDelete={setDeleteProject}
                onProjectClick={handleProjectClick}
                maxHeight={500}
              />
            )}
          </Paper>
        </Stack>
      </Container>

      {/* Select Parent Dialog */}
      <SelectParentDialog
        open={showSelectParent}
        onClose={() => setShowSelectParent(false)}
        campaignId={campaignId}
      />

      {/* Batch Import Dialog */}
      {parentProject && parentFiles?.coordinates?.[0] && (
        <BatchImportDialog
          open={showBatchImport}
          onClose={() => setShowBatchImport(false)}
          onSuccess={() => mutateMemberProjects()}
          campaignId={campaignId}
          parentProjectId={parentProject.id}
          latestCoordsFileId={parentFiles.coordinates[parentFiles.coordinates.length - 1].id}
        />
      )}

      {/* Coords Import Dialog */}
      {parentProject && (
        <CoordsImportDialog
          open={showImportCoords}
          onClose={() => {
            setShowImportCoords(false);
            mutateParentFiles();
          }}
          parentProjectId={parentProject.id}
          parentProjectName={parentProject.name}
        />
      )}

      {/* FreeR Import Dialog */}
      {parentProject && (
        <FreeRImportDialog
          open={showImportFreeR}
          onClose={() => {
            setShowImportFreeR(false);
            mutateParentFiles();
          }}
          parentProjectId={parentProject.id}
          parentProjectName={parentProject.name}
        />
      )}

      {/* Delete Project Dialog */}
      <Dialog open={!!deleteProject} onClose={() => setDeleteProject(null)}>
        <DialogTitle>Delete Project</DialogTitle>
        <DialogContent>
          <Typography>
            Delete project &quot;{deleteProject?.name}&quot; and all associated jobs?
          </Typography>
          <Typography variant="body2" color="error" sx={{ mt: 1 }}>
            This action cannot be undone.
          </Typography>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setDeleteProject(null)}>Cancel</Button>
          <Button onClick={handleDeleteProject} color="error" variant="contained">
            Delete
          </Button>
        </DialogActions>
      </Dialog>

      {/* PANDDA Export Dialog */}
      <PanddaExportDialog
        open={showPanddaExport}
        onClose={() => setShowPanddaExport(false)}
        campaignId={campaignId}
        campaignName={campaign.name}
      />

      {/* Edit File Annotation Dialog */}
      <Dialog
        open={!!editingFile}
        onClose={() => setEditingFile(null)}
        maxWidth="sm"
        fullWidth
      >
        <DialogTitle>Edit File Annotation</DialogTitle>
        <DialogContent>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
            File: {editingFile?.name}
          </Typography>
          <TextField
            autoFocus
            fullWidth
            label="Annotation"
            placeholder="Enter a descriptive label for this file"
            value={editAnnotation}
            onChange={(e) => setEditAnnotation(e.target.value)}
            helperText="This label will be displayed instead of the filename"
          />
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setEditingFile(null)} disabled={savingAnnotation}>
            Cancel
          </Button>
          <Button
            onClick={handleSaveAnnotation}
            variant="contained"
            disabled={savingAnnotation}
            startIcon={savingAnnotation ? <CircularProgress size={16} /> : null}
          >
            {savingAnnotation ? "Saving..." : "Save"}
          </Button>
        </DialogActions>
      </Dialog>
    </>
  );
}

