/*
 * Copyright (C) 2025-2026 Newcastle University
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
import { moorhen } from "moorhen/types/moorhen";
import {
  hideMap,
  showMap,
  hideMolecule,
  showMolecule,
  removeMolecule,
  setRequestDrawScene,
} from "moorhen";
import { useDispatch, useSelector } from "react-redux";
import { useCallback, useEffect, useMemo, useState } from "react";
import {
  Box,
  IconButton,
  Menu,
  MenuItem,
  Slider,
  Stack,
  Tooltip,
  Typography,
  Dialog,
} from "@mui/material";
import {
  VisibilityOutlined,
  VisibilityOffOutlined,
  Delete as DeleteIcon,
  MoreVert,
  Visibility,
  VisibilityOff,
  Science as ScienceIcon,
} from "@mui/icons-material";
import { CCP4i2HierarchyBrowser } from "./ccp4i2-hierarchy-browser";
import { CopyViewLinkButton } from "./copy-view-link-button";
import { PasteViewLinkField } from "./paste-view-link-field";
import { PushToCCP4i2Panel } from "./push-to-ccp4i2-panel";
import { useTheme } from "../../theme/theme-provider";

// Representation types and their compact labels/colours
const REP_DEFS = [
  { key: "CBs", label: "B", title: "Bonds", color: "#1976d2" },
  { key: "CRs", label: "R", title: "Ribbons", color: "#388e3c" },
  { key: "MolecularSurface", label: "S", title: "Surface", color: "#7b1fa2" },
] as const;

interface MoorhenControlPanelProps {
  onFileSelect: (fileId: number) => Promise<void>;
  onJobLoad?: (jobId: number) => Promise<void>;
  getViewUrl?: () => string;
  molecules: moorhen.Molecule[];
  maps: moorhen.Map[];
  onMapContourLevelChange: (molNo: number, level: number) => void;
  onRunServalcat?: (mol: moorhen.Molecule) => Promise<void>;
  servalcatStatus?: string | null;
}

export const MoorhenControlPanel: React.FC<MoorhenControlPanelProps> = ({
  onFileSelect,
  onJobLoad,
  getViewUrl,
  molecules,
  maps,
  onMapContourLevelChange,
  onRunServalcat,
  servalcatStatus,
}) => {
  const { customColors } = useTheme();
  const dispatch = useDispatch();

  const cootInitialized = useSelector(
    (state: moorhen.State) => state.generalStates.cootInitialized
  );
  const visibleMolecules = useSelector(
    (state: moorhen.State) => state.molecules.visibleMolecules || []
  );
  const visibleMaps = useSelector(
    (state: moorhen.State) => state.mapContourSettings?.visibleMaps || []
  );
  const contourLevels = useSelector(
    (state: moorhen.State) => state.mapContourSettings?.contourLevels || []
  );

  // Per-molecule representation state: Map<molNo, Set<repKey>>
  const [molReps, setMolReps] = useState<Map<number, Set<string>>>(new Map());

  // Context menu state for molecules
  const [menuState, setMenuState] = useState<{
    anchorEl: HTMLElement | null;
    mol: moorhen.Molecule | null;
  }>({ anchorEl: null, mol: null });

  // Push to CCP4i2 dialog
  const [pushMol, setPushMol] = useState<moorhen.Molecule | null>(null);

  // Initialize reps for newly loaded molecules (default is CRs from fetchMolecule)
  useEffect(() => {
    if (!molecules || molecules.length === 0) return;
    setMolReps((prev) => {
      const next = new Map(prev);
      let changed = false;
      for (const mol of molecules) {
        if (mol.molNo != null && !next.has(mol.molNo)) {
          next.set(mol.molNo, new Set(["CRs"]));
          changed = true;
        }
      }
      // Clean up entries for molecules that no longer exist
      const activeMolNos = new Set(molecules.map((m) => m.molNo).filter((n) => n != null));
      for (const molNo of next.keys()) {
        if (!activeMolNos.has(molNo)) {
          next.delete(molNo);
          changed = true;
        }
      }
      return changed ? next : prev;
    });
  }, [molecules]);

  const toggleRep = useCallback(
    async (mol: moorhen.Molecule, repKey: string) => {
      if (mol.molNo == null) return;
      const current = molReps.get(mol.molNo) || new Set<string>();
      const isActive = current.has(repKey);

      if (isActive) {
        try {
          mol.clearBuffersOfStyle(repKey);
        } catch (err) {
          console.error(`Failed to remove ${repKey} from ${mol.name}:`, err);
        }
        const next = new Set(current);
        next.delete(repKey);
        setMolReps((prev) => new Map(prev).set(mol.molNo!, next));
      } else {
        try {
          await mol.addRepresentation(repKey, "/*/*/*/*");
        } catch (err) {
          console.error(`Failed to add ${repKey} to ${mol.name}:`, err);
        }
        const next = new Set(current);
        next.add(repKey);
        setMolReps((prev) => new Map(prev).set(mol.molNo!, next));
      }
    },
    [molReps]
  );

  const getContourLevel = useCallback(
    (molNo: number): number => {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      const entry = contourLevels.find((c: any) => c.molNo === molNo);
      const level = entry?.contourLevel;
      if (level === undefined || level === null || Number.isNaN(level)) return 1.0;
      return level;
    },
    [contourLevels]
  );

  // Molecule context menu handlers
  const handleMenuOpen = useCallback(
    (event: React.MouseEvent<HTMLElement>, mol: moorhen.Molecule) => {
      event.stopPropagation();
      setMenuState({ anchorEl: event.currentTarget, mol });
    },
    []
  );
  const handleMenuClose = useCallback(() => {
    setMenuState({ anchorEl: null, mol: null });
  }, []);
  const handleCentre = useCallback(() => {
    if (menuState.mol) {
      menuState.mol.centreOn("/*/*/*/*", false, true);
    }
    handleMenuClose();
  }, [menuState.mol, handleMenuClose]);
  const handlePush = useCallback(() => {
    if (menuState.mol) setPushMol(menuState.mol);
    handleMenuClose();
  }, [menuState.mol, handleMenuClose]);
  const handleServalcat = useCallback(() => {
    if (menuState.mol && onRunServalcat) {
      onRunServalcat(menuState.mol);
    }
    handleMenuClose();
  }, [menuState.mol, onRunServalcat, handleMenuClose]);

  if (!cootInitialized) {
    return (
      <Box
        sx={{
          height: "calc(100vh - 75px)",
          overflowY: "auto",
          display: "flex",
          alignItems: "center",
          justifyContent: "center",
          maxWidth: "100%",
        }}
      >
        Loading...
      </Box>
    );
  }

  return (
    <Stack
      direction="column"
      sx={{
        height: "calc(100vh - 100px)",
        width: "100%",
      }}
    >
      {/* Copy / Paste View Link */}
      {getViewUrl && (
        <Box
          sx={{
            p: 1,
            borderBottom: `1px solid ${customColors.ui.mediumGray}`,
            backgroundColor: customColors.ui.lightGray,
          }}
        >
          <Stack direction="row" spacing={0.5} sx={{ alignItems: "center" }}>
            <CopyViewLinkButton getViewUrl={getViewUrl} disabled={!cootInitialized} />
            <Box sx={{ flex: 1 }}>
              <PasteViewLinkField />
            </Box>
          </Stack>
        </Box>
      )}

      {/* Upper section - CCP4i2HierarchyBrowser */}
      <Box
        sx={{
          flex: 1,
          minHeight: 0,
          overflow: "hidden",
          display: "flex",
          flexDirection: "column",
        }}
      >
        <CCP4i2HierarchyBrowser onFileSelect={onFileSelect} onJobLoad={onJobLoad} />
      </Box>

      {/* Molecules section */}
      {molecules && molecules.length > 0 && (
        <Box
          sx={{
            minHeight: 0,
            overflow: "hidden",
            backgroundColor: customColors.ui.lightGray,
            border: `1px solid ${customColors.ui.mediumGray}`,
            borderTop: `2px solid ${customColors.ui.lightBlue}`,
            position: "relative",
            display: "flex",
            flexDirection: "column",
          }}
        >
          <Typography
            variant="caption"
            sx={{
              position: "absolute",
              top: 2,
              left: 6,
              fontSize: "0.7rem",
              fontWeight: 500,
              color: customColors.ui.lightBlue,
              backgroundColor: "rgba(255, 255, 255, 0.8)",
              px: 0.5,
              py: 0.25,
              borderRadius: "2px",
              zIndex: 1,
            }}
          >
            Molecules
          </Typography>
          <Box sx={{ overflow: "auto", pt: 2.5, px: 0.5, pb: 0.5 }}>
            {molecules.map((mol) => {
              const isVisible = mol.molNo != null && visibleMolecules.includes(mol.molNo);
              const reps = mol.molNo != null ? molReps.get(mol.molNo) || new Set<string>() : new Set<string>();
              return (
                <Stack
                  key={mol.molNo ?? mol.uniqueId}
                  direction="row"
                  alignItems="center"
                  spacing={0.25}
                  sx={{
                    mb: 0.25,
                    opacity: isVisible ? 1 : 0.45,
                    px: 0.5,
                    py: 0.25,
                    borderRadius: 0.5,
                    "&:hover": { backgroundColor: "rgba(0,0,0,0.04)" },
                  }}
                >
                  {/* Molecule name */}
                  <Typography
                    variant="caption"
                    noWrap
                    sx={{
                      flex: 1,
                      minWidth: 0,
                      fontSize: "0.7rem",
                      cursor: "default",
                    }}
                    title={mol.name}
                  >
                    {mol.name}
                  </Typography>

                  {/* Per-molecule representation toggles */}
                  {REP_DEFS.map(({ key, label, title, color }) => {
                    const active = reps.has(key);
                    return (
                      <Tooltip key={key} title={title} enterDelay={400}>
                        <Box
                          component="span"
                          onClick={() => toggleRep(mol, key)}
                          sx={{
                            display: "inline-flex",
                            alignItems: "center",
                            justifyContent: "center",
                            width: 18,
                            height: 18,
                            borderRadius: "3px",
                            fontSize: "0.65rem",
                            fontWeight: 700,
                            cursor: "pointer",
                            userSelect: "none",
                            color: active ? "#fff" : "text.disabled",
                            backgroundColor: active ? color : "transparent",
                            border: `1px solid ${active ? color : "rgba(0,0,0,0.15)"}`,
                            transition: "all 0.15s",
                            "&:hover": {
                              borderColor: color,
                              backgroundColor: active ? color : `${color}20`,
                            },
                          }}
                        >
                          {label}
                        </Box>
                      </Tooltip>
                    );
                  })}

                  {/* Visibility toggle */}
                  <Tooltip title={isVisible ? "Hide" : "Show"} enterDelay={400}>
                    <IconButton
                      size="small"
                      onClick={() => {
                        if (isVisible) {
                          mol.representations?.forEach((r: any) => r.hide());
                          dispatch(hideMolecule(mol as any));
                        } else {
                          mol.representations?.forEach((r: any) => {
                            if (r.interfaceOption?.visible) r.show();
                          });
                          dispatch(showMolecule(mol as any));
                        }
                        dispatch(setRequestDrawScene(true));
                      }}
                      sx={{ p: 0.25 }}
                    >
                      {isVisible ? (
                        <VisibilityOutlined sx={{ fontSize: 14 }} />
                      ) : (
                        <VisibilityOffOutlined sx={{ fontSize: 14, opacity: 0.4 }} />
                      )}
                    </IconButton>
                  </Tooltip>

                  {/* Delete */}
                  <Tooltip title="Remove" enterDelay={400}>
                    <IconButton
                      size="small"
                      onClick={() => {
                        dispatch(removeMolecule(mol));
                        mol.delete();
                        dispatch(setRequestDrawScene(true));
                      }}
                      sx={{ p: 0.25 }}
                      color="error"
                    >
                      <DeleteIcon sx={{ fontSize: 14 }} />
                    </IconButton>
                  </Tooltip>

                  {/* More menu (Centre, Push to CCP4i2) */}
                  <IconButton
                    size="small"
                    onClick={(e) => handleMenuOpen(e, mol)}
                    sx={{ p: 0.25, opacity: 0.6, "&:hover": { opacity: 1 } }}
                  >
                    <MoreVert sx={{ fontSize: 14 }} />
                  </IconButton>
                </Stack>
              );
            })}
          </Box>

          {/* Molecule context menu */}
          <Menu
            anchorEl={menuState.anchorEl}
            open={Boolean(menuState.anchorEl)}
            onClose={handleMenuClose}
            slotProps={{
              paper: { sx: { minWidth: 160, boxShadow: "0 4px 12px rgba(0,0,0,0.15)" } },
            }}
          >
            <MenuItem onClick={handleCentre}>
              <Typography sx={{ mr: 1 }} variant="body2">
                Centre on molecule
              </Typography>
            </MenuItem>
            <MenuItem onClick={handlePush}>
              <Typography variant="body2">Push to CCP4i2</Typography>
            </MenuItem>
            {onRunServalcat && (
              <MenuItem onClick={handleServalcat} disabled={!!servalcatStatus}>
                <ScienceIcon sx={{ fontSize: 16, mr: 1 }} />
                <Typography variant="body2">Run Servalcat</Typography>
              </MenuItem>
            )}
          </Menu>

          {/* Push to CCP4i2 Dialog */}
          <Dialog
            open={pushMol !== null}
            onClose={() => setPushMol(null)}
            maxWidth="sm"
            fullWidth
          >
            {pushMol && (
              <PushToCCP4i2Panel
                molNo={pushMol.molNo ?? undefined}
                item={pushMol}
                onClose={() => setPushMol(null)}
              />
            )}
          </Dialog>

          {/* Servalcat status */}
          {servalcatStatus && (
            <Typography
              variant="caption"
              sx={{
                px: 1,
                py: 0.5,
                fontSize: "0.65rem",
                color: servalcatStatus.startsWith("Failed")
                  ? "error.main"
                  : "text.secondary",
                borderTop: `1px solid ${customColors.ui.mediumGray}`,
                backgroundColor: "rgba(0,0,0,0.02)",
              }}
            >
              {servalcatStatus}
            </Typography>
          )}
        </Box>
      )}

      {/* Maps section */}
      {maps && maps.length > 0 && (
        <Box
          sx={{
            minHeight: 0,
            overflow: "hidden",
            backgroundColor: customColors.ui.lightGray,
            border: `1px solid ${customColors.ui.mediumGray}`,
            borderTop: `2px solid ${customColors.ui.lightBlue}`,
            position: "relative",
            display: "flex",
            flexDirection: "column",
          }}
        >
          <Typography
            variant="caption"
            sx={{
              position: "absolute",
              top: 2,
              left: 6,
              fontSize: "0.7rem",
              fontWeight: 500,
              color: customColors.ui.lightBlue,
              backgroundColor: "rgba(255, 255, 255, 0.8)",
              px: 0.5,
              py: 0.25,
              borderRadius: "2px",
              zIndex: 1,
            }}
          >
            Maps
          </Typography>
          <Box sx={{ overflow: "auto", pt: 2.5, px: 0.5, pb: 0.5 }}>
            {maps.map((map) => {
              const level = getContourLevel(map.molNo!);
              const isVisible = visibleMaps.includes(map.molNo!);
              const mapSubType = (map as any).mapSubType as number | undefined;
              const isDiff = map.isDifference;
              const multiplier = isDiff ? 2 : 1;
              const minLog = -3;
              const maxLog = 1.0;
              const minValue = Math.pow(10, minLog) * multiplier;
              const maxValue = Math.pow(10, maxLog) * multiplier;

              const valueToSlider = (v: number) => {
                const clampedV = Math.max(minValue, Math.min(maxValue, v));
                return ((Math.log10(clampedV / multiplier) - minLog) / (maxLog - minLog)) * 100;
              };
              const sliderToValue = (s: number) => {
                const logValue = minLog + (s / 100) * (maxLog - minLog);
                return Math.pow(10, logValue) * multiplier;
              };

              const sliderPosition = valueToSlider(level);
              const shortName =
                mapSubType === 3
                  ? "Anom"
                  : mapSubType === 2
                    ? "Fo-Fc"
                    : isDiff
                      ? "Fo-Fc"
                      : "2Fo-Fc";

              return (
                <Stack
                  key={map.molNo ?? map.uniqueId}
                  direction="row"
                  alignItems="center"
                  spacing={0.5}
                  sx={{ mb: 0.25, px: 0.5 }}
                >
                  <Typography
                    variant="caption"
                    sx={{
                      minWidth: 38,
                      flexShrink: 0,
                      opacity: isVisible ? 1 : 0.4,
                      fontSize: "0.7rem",
                    }}
                  >
                    {shortName}
                  </Typography>
                  <Slider
                    size="small"
                    disabled={!isVisible}
                    value={sliderPosition}
                    onChange={(_e, value) =>
                      onMapContourLevelChange(map.molNo!, sliderToValue(value as number))
                    }
                    min={0}
                    max={100}
                    step={1}
                    valueLabelDisplay="auto"
                    valueLabelFormat={(v) => {
                      const actualValue = sliderToValue(v);
                      return actualValue < 0.01
                        ? actualValue.toExponential(1)
                        : actualValue.toFixed(3);
                    }}
                    sx={{ flex: 1, py: 0, mx: 0.5 }}
                  />
                  <Typography
                    variant="caption"
                    color="text.secondary"
                    sx={{
                      minWidth: 40,
                      textAlign: "right",
                      flexShrink: 0,
                      opacity: isVisible ? 1 : 0.4,
                      fontSize: "0.7rem",
                    }}
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
                      {isVisible ? (
                        <VisibilityOutlined sx={{ fontSize: 14 }} />
                      ) : (
                        <VisibilityOffOutlined sx={{ fontSize: 14, opacity: 0.4 }} />
                      )}
                    </IconButton>
                  </Tooltip>
                </Stack>
              );
            })}
          </Box>
        </Box>
      )}
    </Stack>
  );
};
