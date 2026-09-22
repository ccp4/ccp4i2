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

import React, { useState, useCallback, useMemo, useEffect } from "react";
import { useDispatch, useSelector } from "react-redux";
import {
  hideMap,
  showMap,
  hideMolecule,
  showMolecule,
  removeMolecule,
  setRequestDrawScene,
} from "moorhen/react-lib";
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
  ToggleButton,
  ToggleButtonGroup,
  Slider,
  Stack,
} from "@mui/material";
import {
  Place as PlaceIcon,
  Preview as PreviewIcon,
  Add as AddIcon,
  Delete as DeleteIcon,
  Edit as EditIcon,
  Home as HomeIcon,
  Upload as UploadIcon,
  CheckCircle as HitIcon,
  HelpOutline as UnclearIcon,
  DoDisturbOn as NothingThereIcon,
  FolderOpen as FolderOpenIcon,
  Science as ScienceIcon,
  VisibilityOutlined,
  VisibilityOffOutlined,
} from "@mui/icons-material";
import { moorhen } from "moorhen/types/moorhen";
import { CopyViewLinkButton } from "./copy-view-link-button";
import { PasteViewLinkField } from "./paste-view-link-field";
import { PushToCCP4i2Panel } from "./push-to-ccp4i2-panel";
import { CCP4i2HierarchyBrowser } from "./ccp4i2-hierarchy-browser";
import { Ligand2DView } from "../campaigns/ligand-2d-view";
import { AddLigandButton } from "./add-ligand-button";
import {
  ProjectGroup,
  CampaignSite,
  MemberProjectWithSummary,
  SiteEvaluation,
  SiteVerdict,
} from "../../types/campaigns";
import { Project } from "../../types/models";

/**
 * The floor for the two lists that are navigated (sites) and acted on
 * (loaded molecules). Both grow with the campaign, so left to themselves one
 * squeezes the other out of the panel: with a dozen projects loaded the sites
 * list collapsed to a sliver and the campaign could no longer be navigated.
 * Each keeps at least this much and scrolls within it; if the panel is too
 * short for both floors, the panel itself scrolls.
 */
const MIN_LIST_SECTION_HEIGHT = 160;

/** The contour sliders are set once and left alone, so they yield space to the
 *  lists that are worked in, and scroll past this. */
const MAX_CONTOUR_SECTION_HEIGHT = 180;

interface CampaignControlPanelProps {
  campaign: ProjectGroup;
  sites: CampaignSite[];
  onGoToSite: (site: CampaignSite) => void;
  /** Open the site's scene: every hit recorded there, fitted on the pocket.
   *  Going TO a site moves the camera in what is loaded; VIEWING a site
   *  loads what was found there. Different acts, so different controls. */
  onViewSite?: (site: CampaignSite) => void;
  onSaveCurrentAsSite: (name: string) => Promise<void>;
  onUpdateSite: (siteId: number, name: string, updatePosition: boolean) => Promise<void>;
  onDeleteSite: (siteId: number) => Promise<void>;
  memberProjects: MemberProjectWithSummary[];
  selectedMemberProjectId: number | null;
  onSelectMemberProject: (projectId: number | null) => void;
  parentProject: Project | null | undefined;
  getViewUrl: () => string;
  /** Molecules loaded in Moorhen for push-to-CCP4i2 */
  molecules?: moorhen.Molecule[];
  /** Controlled representation state */
  visibleRepresentations: string[];
  onRepresentationsChange: (representations: string[]) => void;
  /** Ligand dictionary file ID for 2D structure display */
  ligandDictFileId?: number | null;
  /** Ligand name for display */
  ligandName?: string | null;
  /** The codes the member's dictionaries define, for "Add ligand here". */
  ligandCodes?: string[];
  /** The file the member's coordinates were loaded from: the molecule the
   *  ligand is added to is found by this, never by which one is active. */
  memberCoordFileId?: number | null;
  /** Place one ligand code at the view centre in the member's molecule. */
  onAddLigand?: (code: string) => Promise<void>;
  /** Maps loaded in Moorhen for contour control */
  maps?: moorhen.Map[];
  /** Callback to change map contour level */
  onMapContourLevelChange?: (molNo: number, level: number) => void;
  /**
   * What has been found at each site in the selected dataset, empties
   * included. An absent entry means nobody has looked there yet, which is a
   * different thing from having looked and found nothing.
   */
  evaluations?: SiteEvaluation[];
  /** Record or change what was found at a site in the selected dataset. */
  onSetVerdict?: (siteId: number, verdict: SiteVerdict) => Promise<void>;
  /** Withdraw a verdict, returning the site to "nobody has looked". */
  onClearVerdict?: (siteId: number) => Promise<void>;
  /** Callback to load a file into the Moorhen session */
  onFileSelect?: (fileId: number) => Promise<void>;
  /** Callback to load all of a job's outputs into the session (the browser's
   *  per-job load button, shown only when this is given). */
  onJobLoad?: (jobId: number) => Promise<void>;
  /** Callback to run servalcat_pipe refinement on a molecule */
  onRunServalcat?: (mol: moorhen.Molecule) => Promise<void>;
}

export const CampaignControlPanel: React.FC<CampaignControlPanelProps> = ({
  campaign,
  sites,
  onGoToSite,
  onViewSite,
  onSaveCurrentAsSite,
  onUpdateSite,
  onDeleteSite,
  memberProjects,
  selectedMemberProjectId,
  onSelectMemberProject,
  parentProject,
  getViewUrl,
  molecules,
  visibleRepresentations,
  onRepresentationsChange,
  ligandDictFileId,
  ligandName,
  ligandCodes,
  memberCoordFileId,
  onAddLigand,
  maps,
  onMapContourLevelChange,
  evaluations,
  onSetVerdict,
  onClearVerdict,
  onFileSelect,
  onJobLoad,
  onRunServalcat,
}) => {
  const dispatch = useDispatch();

  // Get contour levels and visibility from Redux state
  const contourLevels = useSelector(
    (state: moorhen.State) => state.mapContourSettings?.contourLevels || []
  );
  const visibleMaps = useSelector(
    (state: moorhen.State) => state.mapContourSettings?.visibleMaps || []
  );
  const visibleMolecules = useSelector(
    (state: moorhen.State) => state.molecules.visibleMolecules || []
  );

  // Helper to get contour level for a specific map
  // Note: Moorhen internally uses { molNo, contourLevel } despite TypeScript defs saying 'level'
  const getContourLevel = useCallback((molNo: number): number => {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const entry = contourLevels.find((c: any) => c.molNo === molNo);
    const level = entry?.contourLevel;
    // Guard against NaN/undefined - return sensible default
    if (level === undefined || level === null || Number.isNaN(level)) {
      return 1.0;
    }
    return level;
  }, [contourLevels]);

  const [showAddSiteDialog, setShowAddSiteDialog] = useState(false);
  const [newSiteName, setNewSiteName] = useState("");
  const [isSaving, setIsSaving] = useState(false);

  // Edit site dialog state
  const [editingSiteId, setEditingSiteId] = useState<number | null>(null);
  // The site awaiting a delete confirmation. Held as the site itself, not an
  // id, so the dialog can name it and say what it is about to destroy.
  const [sitePendingDelete, setSitePendingDelete] = useState<CampaignSite | null>(
    null
  );
  const [editSiteName, setEditSiteName] = useState("");
  const [updatePosition, setUpdatePosition] = useState(false);

  // Push to CCP4i2 dialog state
  const [showPushDialog, setShowPushDialog] = useState(false);
  const [selectedMolecule, setSelectedMolecule] = useState<moorhen.Molecule | null>(null);

  // Ligand code loaded from CIF file (overrides ligandName prop)
  const [loadedLigandCode, setLoadedLigandCode] = useState<string | null>(null);

  // Import from projects modal state
  const [showImportModal, setShowImportModal] = useState(false);

  // Track which molecule (if any) has a servalcat job running
  const [runningServalcatMolNo, setRunningServalcatMolNo] = useState<number | null>(null);

  // Reset loaded ligand code when file ID changes
  useEffect(() => {
    setLoadedLigandCode(null);
  }, [ligandDictFileId]);

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

  // Handle file selection from the hierarchy browser
  const handleImportFileSelect = useCallback(
    async (fileId: number) => {
      if (onFileSelect) {
        await onFileSelect(fileId);
      }
      setShowImportModal(false);
    },
    [onFileSelect]
  );

  // Handle "load all job outputs" from the hierarchy browser
  const handleImportJobLoad = useCallback(
    async (jobId: number) => {
      if (onJobLoad) {
        await onJobLoad(jobId);
      }
      setShowImportModal(false);
    },
    [onJobLoad]
  );

  const handleOpenPushDialog = useCallback((molecule: moorhen.Molecule) => {
    setSelectedMolecule(molecule);
    setShowPushDialog(true);
  }, []);

  const handleClosePushDialog = useCallback(() => {
    setShowPushDialog(false);
    setSelectedMolecule(null);
  }, []);

  const handleToggleRepresentation = useCallback(
    async (_event: React.MouseEvent<HTMLElement>, newRepresentations: string[]) => {
      if (!molecules) return;

      const previousReps = visibleRepresentations;
      onRepresentationsChange(newRepresentations);

      // Determine which representations were added or removed
      const added = newRepresentations.filter((r) => !previousReps.includes(r));
      const removed = previousReps.filter((r) => !newRepresentations.includes(r));

      for (const mol of molecules) {
        // Add new representations
        for (const rep of added) {
          try {
            // `added` are ToggleButton representation keys typed loosely as
            // string[] through MUI onChange; narrow to the Moorhen union.
            await mol.addRepresentation(rep as moorhen.RepresentationStyles, "/*/*/*/*");
          } catch (err) {
            console.error(`Failed to add ${rep} for molecule ${mol.name}:`, err);
          }
        }
        // Remove old representations
        for (const rep of removed) {
          try {
            mol.clearBuffersOfStyle(rep);
          } catch (err) {
            console.error(`Failed to remove ${rep} for molecule ${mol.name}:`, err);
          }
        }
      }
    },
    [molecules, visibleRepresentations, onRepresentationsChange]
  );

  // What has been recorded at each site for the selected dataset. Absent
  // means nobody has looked there; "empty" means somebody looked and found
  // nothing. The control has to show those differently or it will write the
  // wrong one back.
  const verdictBySite = useMemo(
    () => new Map((evaluations ?? []).map((e) => [e.site_id, e.verdict])),
    [evaluations]
  );

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

  // Deleting a site is not undoable and takes every verdict recorded there
  // with it, in every dataset. It also used to happen on one click of an icon
  // sitting between "go to site" and "edit site" in a dense row, inside a
  // panel whose whole purpose is clicking through sites one after another --
  // so it was deleted by people who meant to navigate. It asks now.
  const handleDeleteSite = useCallback((site: CampaignSite) => {
    setSitePendingDelete(site);
  }, []);

  const handleConfirmDeleteSite = useCallback(async () => {
    if (!sitePendingDelete) return;
    const siteId = sitePendingDelete.id;
    setSitePendingDelete(null);
    await onDeleteSite(siteId);
  }, [sitePendingDelete, onDeleteSite]);

  const handleOpenEditDialog = useCallback((site: CampaignSite) => {
    setEditingSiteId(site.id);
    setEditSiteName(site.name);
    setUpdatePosition(false);
  }, []);

  const handleCloseEditDialog = useCallback(() => {
    setEditingSiteId(null);
    setEditSiteName("");
    setUpdatePosition(false);
  }, []);

  const handleSaveEdit = useCallback(async () => {
    if (editingSiteId === null || !editSiteName.trim()) return;
    setIsSaving(true);
    try {
      await onUpdateSite(editingSiteId, editSiteName.trim(), updatePosition);
      handleCloseEditDialog();
    } finally {
      setIsSaving(false);
    }
  }, [editingSiteId, editSiteName, updatePosition, onUpdateSite, handleCloseEditDialog]);

  return (
    <Box sx={{ p: 1.5, height: "100%", display: "flex", flexDirection: "column", overflowY: "auto", overflowX: "hidden" }}>
      {/* Header */}
      <Typography variant="subtitle1" sx={{ mb: 1, fontWeight: "bold" }}>
        {campaign.name}
      </Typography>

      {/* Copy / Paste View Link */}
      <Stack direction="row" spacing={0.5} sx={{ mb: 1, alignItems: "center" }}>
        <CopyViewLinkButton getViewUrl={getViewUrl} />
        <Box sx={{ flex: 1 }}>
          <PasteViewLinkField />
        </Box>
      </Stack>

      {/* Import from Projects */}
      {onFileSelect && (
        <Box sx={{ mb: 1 }}>
          <Button
            variant="outlined"
            size="small"
            startIcon={<FolderOpenIcon />}
            onClick={() => setShowImportModal(true)}
            fullWidth
          >
            Import from Projects
          </Button>
        </Box>
      )}

      {/* Display Options */}
      {molecules && molecules.length > 0 && (
        <Box sx={{ mb: 1 }}>
          <Typography variant="caption" color="text.secondary" sx={{ mb: 0.5, display: "block" }}>
            Representations
          </Typography>
          <ToggleButtonGroup
            value={visibleRepresentations}
            onChange={handleToggleRepresentation}
            size="small"
            aria-label="molecule representations"
          >
            <ToggleButton value="CBs" aria-label="bonds">
              <Tooltip title="Bonds">
                <span>Bonds</span>
              </Tooltip>
            </ToggleButton>
            <ToggleButton value="CRs" aria-label="ribbons">
              <Tooltip title="Ribbons">
                <span>Ribbons</span>
              </Tooltip>
            </ToggleButton>
            <ToggleButton value="MolecularSurface" aria-label="surface">
              <Tooltip title="Molecular Surface">
                <span>Surface</span>
              </Tooltip>
            </ToggleButton>
          </ToggleButtonGroup>
        </Box>
      )}

      {/* Map Contour Controls */}
      {maps && maps.length > 0 && onMapContourLevelChange && (
        // Two maps per loaded project, so this grows with the campaign as well;
        // left uncapped it pushes the lists below it off the panel entirely.
        <Box sx={{ mb: 1, maxHeight: MAX_CONTOUR_SECTION_HEIGHT, overflowY: "auto", flexShrink: 0 }}>
          {maps.map((map) => {
            const level = getContourLevel(map.molNo!);
            const isVisible = visibleMaps.includes(map.molNo!);
            // mapSubType: 1=normal, 2=difference, 3=anomalous
            // Difference and anomalous maps need 2x higher values for contour slider
            const mapSubType = (map as any).mapSubType as number | undefined;
            const isDiff = map.isDifference;
            const multiplier = isDiff ? 2 : 1;
            // Logarithmic scale: 10^-3 to 10^1 for regular maps (~0.001 to ~10)
            // For diff maps: ~0.002 to ~20
            const minLog = -3;  // 10^-3 = 0.001
            const maxLog = 1.0; // 10^1 = 10
            const minValue = Math.pow(10, minLog) * multiplier;
            const maxValue = Math.pow(10, maxLog) * multiplier;

            // Convert actual value to slider position (0-100)
            const valueToSlider = (v: number) => {
              const clampedV = Math.max(minValue, Math.min(maxValue, v));
              return ((Math.log10(clampedV / multiplier) - minLog) / (maxLog - minLog)) * 100;
            };

            // Convert slider position to actual value
            const sliderToValue = (s: number) => {
              const logValue = minLog + (s / 100) * (maxLog - minLog);
              return Math.pow(10, logValue) * multiplier;
            };

            const sliderPosition = valueToSlider(level);
            // Label based on map sub_type: 1=normal (2Fo-Fc), 2=difference (Fo-Fc), 3=anomalous (Anom), 4=mask (Mask)
            const shortName = mapSubType === 4 ? "Mask" : mapSubType === 3 ? "Anom" : mapSubType === 2 ? "Fo-Fc" : isDiff ? "Fo-Fc" : "2Fo-Fc";
            // The row is labelled by map type, which is what you want while
            // scanning contour sliders and useless for telling apart the two
            // identically-labelled rows a second loaded dataset brings. The
            // hover says which map this actually is: the type spelled out,
            // and the file's own annotation underneath.
            const fullType =
              mapSubType === 4
                ? "Mask"
                : mapSubType === 3
                ? "Anomalous difference map"
                : mapSubType === 2 || isDiff
                ? "Fo-Fc difference map"
                : "2Fo-Fc weighted map";
            const description =
              ((map as any).ccp4i2Description as string | undefined) ||
              map.name ||
              undefined;

            return (
              <Stack key={map.molNo ?? map.uniqueId} direction="row" alignItems="center" spacing={1} sx={{ mb: 0.5 }}>
                <Tooltip
                  title={
                    description && description !== fullType ? (
                      <Box>
                        <Typography variant="body2">{fullType}</Typography>
                        <Typography variant="caption" sx={{ opacity: 0.85 }}>
                          {description}
                        </Typography>
                      </Box>
                    ) : (
                      fullType
                    )
                  }
                >
                  <Typography
                    variant="caption"
                    sx={{
                      minWidth: 42,
                      flexShrink: 0,
                      opacity: isVisible ? 1 : 0.4,
                      // Hints that the abbreviation has more behind it.
                      cursor: "help",
                    }}
                  >
                    {shortName}
                  </Typography>
                </Tooltip>
                <Slider
                  size="small"
                  disabled={!isVisible}
                  value={sliderPosition}
                  onChange={(_e, value) => onMapContourLevelChange(map.molNo!, sliderToValue(value as number))}
                  min={0}
                  max={100}
                  step={1}
                  valueLabelDisplay="auto"
                  valueLabelFormat={(v) => {
                    const actualValue = sliderToValue(v);
                    return actualValue < 0.01 ? actualValue.toExponential(1) : actualValue.toFixed(3);
                  }}
                  sx={{ flex: 1, py: 0, mx: 0.5 }}
                />
                <Typography
                  variant="caption"
                  color="text.secondary"
                  sx={{ minWidth: 45, textAlign: "right", flexShrink: 0, opacity: isVisible ? 1 : 0.4 }}
                >
                  {level < 0.01 ? level.toExponential(1) : level.toFixed(2)}
                </Typography>
                <Tooltip title={isVisible ? "Hide map" : "Show map"}>
                  <IconButton
                    size="small"
                    onClick={() => {
                      if (isVisible) {
                        dispatch(hideMap(map as any));
                      } else {
                        dispatch(showMap(map as any));
                      }
                      dispatch(setRequestDrawScene(true));
                    }}
                    sx={{ p: 0.25, flexShrink: 0 }}
                  >
                    {isVisible
                      ? <VisibilityOutlined sx={{ fontSize: 16 }} />
                      : <VisibilityOffOutlined sx={{ fontSize: 16, opacity: 0.4 }} />
                    }
                  </IconButton>
                </Tooltip>
              </Stack>
            );
          })}
        </Box>
      )}

      <Divider sx={{ my: 1 }} />

      {/* Member Project Selector */}
      <FormControl fullWidth size="small" sx={{ mb: 1 }}>
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

      {/* Ligand 2D Structure */}
      {ligandDictFileId && (
        <>
          <Divider sx={{ my: 1 }} />
          <Stack direction="row" alignItems="center" spacing={1} sx={{ mb: 1 }}>
            <Typography variant="caption" color="text.secondary" sx={{ flexShrink: 0 }}>
              {loadedLigandCode || ligandName || "Ligand"}
            </Typography>
            <Box sx={{ flex: 1, display: "flex", justifyContent: "center" }}>
              <Ligand2DView
                fileId={ligandDictFileId}
                name={ligandName || undefined}
                width={140}
                height={90}
                onLigandCodeLoaded={setLoadedLigandCode}
              />
            </Box>
          </Stack>
        </>
      )}

      {/* Shown for any selected dataset, dictionary or not, so the disabled
          state can say why rather than the button silently not existing. */}
      {selectedMemberProjectId && onAddLigand && (
        <Box sx={{ mb: 1 }}>
          <AddLigandButton
            ligandCodes={ligandCodes ?? []}
            molecules={molecules ?? []}
            memberCoordFileId={memberCoordFileId ?? null}
            onAddLigand={onAddLigand}
          />
        </Box>
      )}

      <Divider sx={{ my: 1 }} />

      {/* Sites Section */}
      <Box
        sx={{
          flex: "1 1 0",
          display: "flex",
          flexDirection: "column",
          minHeight: MIN_LIST_SECTION_HEIGHT,
        }}
      >
        <Box
          sx={{
            display: "flex",
            alignItems: "center",
            justifyContent: "space-between",
            mb: 0.5,
          }}
        >
          <Typography variant="caption" sx={{ fontWeight: "bold" }}>
            Binding Sites
          </Typography>
          <Tooltip title="Save current view as a site">
            <IconButton
              size="small"
              onClick={() => setShowAddSiteDialog(true)}
              color="primary"
              sx={{ p: 0.5 }}
            >
              <AddIcon fontSize="small" />
            </IconButton>
          </Tooltip>
        </Box>

        {sites.length === 0 ? (
          <Paper
            variant="outlined"
            sx={{
              p: 1,
              textAlign: "center",
              color: "text.secondary",
              bgcolor: "background.default",
            }}
          >
            <PlaceIcon sx={{ fontSize: 24, opacity: 0.5 }} />
            <Typography variant="caption" display="block">
              No sites saved. Navigate and click + to save.
            </Typography>
          </Paper>
        ) : (
          <List dense sx={{ flex: 1, overflow: "auto" }}>
            {sites.map((site) => (
              <ListItem
                key={site.id}
                component="div"
                onClick={() => onGoToSite(site)}
                sx={{
                  bgcolor: "background.paper",
                  mb: 0.5,
                  borderRadius: 1,
                  border: "1px solid",
                  borderColor: "divider",
                  cursor: "pointer",
                  "&:hover": {
                    bgcolor: "action.hover",
                  },
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
                  {onViewSite && (
                    <Tooltip title="View this site: every hit found here, overlaid">
                      <IconButton
                        edge="end"
                        size="small"
                        onClick={(event) => {
                          // The row itself goes to the site; this does not.
                          event.stopPropagation();
                          onViewSite(site);
                        }}
                        color="primary"
                      >
                        <PreviewIcon fontSize="small" />
                      </IconButton>
                    </Tooltip>
                  )}
                  {selectedMemberProjectId && onSetVerdict && (
                    <SiteVerdictControl
                      verdict={verdictBySite.get(site.id)}
                      onSet={(verdict) => onSetVerdict(site.id, verdict)}
                      onClear={
                        onClearVerdict
                          ? () => onClearVerdict(site.id)
                          : undefined
                      }
                    />
                  )}
                  <Tooltip title="Edit site">
                    <IconButton
                      edge="end"
                      size="small"
                      onClick={() => handleOpenEditDialog(site)}
                    >
                      <EditIcon fontSize="small" />
                    </IconButton>
                  </Tooltip>
                  <Tooltip title="Delete site (and any verdicts recorded there)">
                    <IconButton
                      edge="end"
                      size="small"
                      onClick={(event) => {
                        // The row itself navigates; without this, asking to
                        // delete also moves the view.
                        event.stopPropagation();
                        handleDeleteSite(site);
                      }}
                      color="error"
                      // Set apart from the controls next to it. The row packs
                      // six targets into a panel built for tapping through
                      // sites in sequence, and on a phone they are all well
                      // under the ~44pt a thumb needs -- which is how this
                      // one got hit by someone who meant to navigate.
                      sx={{ ml: 1 }}
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
          <Divider sx={{ my: 1 }} />
          <Box
            sx={{
              flex: "1 1 0",
              display: "flex",
              flexDirection: "column",
              minHeight: MIN_LIST_SECTION_HEIGHT,
            }}
          >
            <Typography variant="caption" sx={{ fontWeight: "bold", mb: 0.5, display: "block" }}>
              Push to CCP4i2
            </Typography>
            <List dense sx={{ flex: 1, overflow: "auto", minHeight: 0 }}>
              {molecules.map((mol) => {
                const isVisible = visibleMolecules.includes(mol.molNo!);
                return (
                  <ListItem
                    key={mol.molNo ?? mol.uniqueId}
                    sx={{
                      bgcolor: "background.paper",
                      mb: 0.5,
                      borderRadius: 1,
                      border: "1px solid",
                      borderColor: "divider",
                      opacity: isVisible ? 1 : 0.5,
                    }}
                  >
                    <ListItemText
                      primary={mol.name}
                      primaryTypographyProps={{ variant: "body2" }}
                    />
                    <ListItemSecondaryAction>
                      <Tooltip title={isVisible ? "Hide molecule" : "Show molecule"}>
                        <IconButton
                          edge="end"
                          size="small"
                          onClick={() => {
                            if (isVisible) {
                              (mol as any).representations?.forEach((r: any) => r.hide());
                              dispatch(hideMolecule(mol as any));
                            } else {
                              (mol as any).representations?.forEach((r: any) => {
                                if (r.interfaceOption?.visible) r.show();
                              });
                              dispatch(showMolecule(mol as any));
                            }
                            dispatch(setRequestDrawScene(true));
                          }}
                        >
                          {isVisible
                            ? <VisibilityOutlined sx={{ fontSize: 18 }} />
                            : <VisibilityOffOutlined sx={{ fontSize: 18 }} />
                          }
                        </IconButton>
                      </Tooltip>
                      <Tooltip title="Remove from Coot">
                        <IconButton
                          edge="end"
                          size="small"
                          onClick={() => {
                            dispatch(removeMolecule(mol));
                            mol.delete();
                            dispatch(setRequestDrawScene(true));
                          }}
                          color="error"
                        >
                          <DeleteIcon fontSize="small" />
                        </IconButton>
                      </Tooltip>
                      {onRunServalcat && (
                        <Tooltip title="Run servalcat refinement">
                          <span>
                            <IconButton
                              edge="end"
                              size="small"
                              disabled={runningServalcatMolNo !== null}
                              onClick={async () => {
                                setRunningServalcatMolNo(mol.molNo);
                                try {
                                  await onRunServalcat(mol);
                                } finally {
                                  setRunningServalcatMolNo(null);
                                }
                              }}
                              color="secondary"
                            >
                              <ScienceIcon fontSize="small" />
                            </IconButton>
                          </span>
                        </Tooltip>
                      )}
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
                );
              })}
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

      {/* Delete Site Confirmation */}
      <Dialog
        open={sitePendingDelete !== null}
        onClose={() => setSitePendingDelete(null)}
        maxWidth="xs"
        fullWidth
      >
        <DialogTitle>Delete Site</DialogTitle>
        <DialogContent>
          <Typography>
            Delete &quot;{sitePendingDelete?.name}&quot; from this campaign?
          </Typography>
          {/* What it costs, not just that it costs something: a site nobody
              has looked at yet is a different proposition from one carrying a
              campaign's worth of verdicts, and a warning that cannot tell
              them apart is one people learn to click through. */}
          {sitePendingDelete?.evaluation_count ? (
            <Typography variant="body2" color="error" sx={{ mt: 1 }}>
              {sitePendingDelete.evaluation_count === 1
                ? "The one verdict recorded here, in any dataset, is deleted with it."
                : `All ${sitePendingDelete.evaluation_count} verdicts recorded here, across every dataset, are deleted with it.`}
            </Typography>
          ) : (
            <Typography variant="body2" color="text.secondary" sx={{ mt: 1 }}>
              No verdicts have been recorded at this site.
            </Typography>
          )}
          <Typography variant="body2" color="error" sx={{ mt: 1 }}>
            This action cannot be undone.
          </Typography>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setSitePendingDelete(null)}>Cancel</Button>
          <Button
            onClick={handleConfirmDeleteSite}
            color="error"
            variant="contained"
          >
            Delete
          </Button>
        </DialogActions>
      </Dialog>

      {/* Edit Site Dialog */}
      <Dialog
        open={editingSiteId !== null}
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
              molNo={selectedMolecule.molNo ?? undefined}
              item={selectedMolecule}
              onClose={handleClosePushDialog}
            />
          )}
        </DialogContent>
      </Dialog>

      {/* Import from Projects Dialog */}
      <Dialog
        open={showImportModal}
        onClose={() => setShowImportModal(false)}
        maxWidth="sm"
        fullWidth
        PaperProps={{
          sx: { height: "70vh", maxHeight: "600px" },
        }}
      >
        <DialogTitle>
          Import from Projects
          <Typography variant="body2" color="text.secondary">
            Navigate to a project and job to load a file, or all of a job&apos;s
            outputs, into this session
          </Typography>
        </DialogTitle>
        <DialogContent sx={{ p: 0, display: "flex", flexDirection: "column" }}>
          <Box sx={{ flex: 1, minHeight: 0 }}>
            <CCP4i2HierarchyBrowser
              onFileSelect={handleImportFileSelect}
              onJobLoad={onJobLoad ? handleImportJobLoad : undefined}
            />
          </Box>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setShowImportModal(false)}>Cancel</Button>
        </DialogActions>
      </Dialog>
    </Box>
  );
};


/**
 * Record what was found at one site in one dataset.
 *
 * This replaces tagging the project with the site's name. A tag carried one
 * bit where three are needed: an untagged project was indistinguishable from
 * one somebody had examined and found nothing in, and the association broke
 * whenever a site was renamed.
 *
 * Clicking the verdict already recorded withdraws it, which is how "nobody has
 * looked" is reachable again without a fourth button. Withdrawing is not the
 * same as recording "nothing there".
 */
const VERDICT_OPTIONS: {
  verdict: SiteVerdict;
  label: string;
  colour: "success" | "warning" | "standard";
  Icon: typeof HitIcon;
}[] = [
  { verdict: "hit", label: "Ligand present", colour: "success", Icon: HitIcon },
  {
    verdict: "unclear",
    label: "Unclear",
    colour: "warning",
    Icon: UnclearIcon,
  },
  {
    verdict: "empty",
    label: "Looked, nothing there",
    colour: "standard",
    Icon: NothingThereIcon,
  },
];

const SiteVerdictControl: React.FC<{
  verdict: SiteVerdict | undefined;
  onSet: (verdict: SiteVerdict) => Promise<void>;
  onClear?: () => Promise<void>;
}> = ({ verdict, onSet, onClear }) => {
  const [busy, setBusy] = useState(false);

  const handle = useCallback(
    async (choice: SiteVerdict) => {
      setBusy(true);
      try {
        if (choice === verdict && onClear) {
          await onClear();
        } else {
          await onSet(choice);
        }
      } finally {
        setBusy(false);
      }
    },
    [verdict, onSet, onClear]
  );

  return (
    <>
      {VERDICT_OPTIONS.map(({ verdict: choice, label, colour, Icon }) => {
        const active = verdict === choice;
        return (
          <Tooltip
            key={choice}
            title={active ? `${label} — click to withdraw` : label}
          >
            <span>
              <IconButton
                edge="end"
                size="small"
                disabled={busy}
                onClick={() => handle(choice)}
                color={active && colour !== "standard" ? colour : "default"}
                sx={{
                  opacity: active ? 1 : 0.35,
                  p: 0.25,
                }}
              >
                <Icon fontSize="small" />
              </IconButton>
            </span>
          </Tooltip>
        );
      })}
    </>
  );
};
