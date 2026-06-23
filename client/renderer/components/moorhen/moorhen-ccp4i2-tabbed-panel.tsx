/**
 * CCP4i2 tabbed side-panel — hosts both Controls and Scenes under one
 * Moorhen side-panel registration.
 *
 * Why: registering two custom extraSidePanels with Moorhen surfaces only
 * one tab (an upstream Moorhen behaviour we couldn't work around in the
 * panel-switcher's local state). MUI Tabs in our own panel content side-
 * steps the issue entirely — one Moorhen side-panel tab, two sub-tabs
 * under our control.
 *
 * The two sub-panels are kept mounted (Tabs hide the inactive one with
 * display:none) so Monaco's editor state survives switching back and
 * forth, and so live-apply continues working while the user is on the
 * Controls tab.
 */
import React, { useState, useMemo } from "react";
import { Box, Tab, Tabs } from "@mui/material";
import { moorhen } from "moorhen/types/moorhen";

import { MoorhenControlPanel } from "./moorhen-control-panel";
import { MoorhenScenesPanel, SceneBundleAssets } from "./moorhen-scenes-panel";
import type { SceneResolveResult } from "../../lib/moorhen-scene-resolver";
import type { SceneLiftHints } from "../../lib/moorhen-scene-lifter";
import type { MoorhenScene } from "../../types/moorhen-scene";

export interface MoorhenCcp4i2TabbedPanelProps {
  onFileSelect: (fileId: number) => Promise<void>;
  onJobLoad?: (jobId: number) => Promise<void>;
  getViewUrl?: () => string;
  molecules: moorhen.Molecule[];
  maps: moorhen.Map[];
  onMapContourLevelChange: (molNo: number, level: number) => void;
  onRunServalcat?: (mol: moorhen.Molecule) => Promise<void>;
  servalcatStatus?: string | null;
  onApplyScene: (
    yamlText: string,
    assets: SceneBundleAssets,
  ) => Promise<SceneResolveResult>;
  onCaptureScene: () => Promise<{
    scene: MoorhenScene;
    hints: SceneLiftHints;
    assets: SceneBundleAssets;
  }>;
  onPromoteSceneToPortable: (
    yamlText: string,
    currentAssets: SceneBundleAssets,
  ) => Promise<{
    yamlText: string;
    assets: SceneBundleAssets;
    warnings: string[];
  }>;
  cootInitialized: boolean;
}

export const MoorhenCcp4i2TabbedPanel: React.FC<MoorhenCcp4i2TabbedPanelProps> = ({
  onApplyScene,
  onCaptureScene,
  onPromoteSceneToPortable,
  cootInitialized,
  ...controlsProps
}) => {
  const [activeTab, setActiveTab] = useState<"controls" | "scenes">("controls");

  // Memoise the two sub-panels so switching tabs doesn't rebuild them
  // (Monaco re-init is expensive; preserving the Controls subtree keeps
  // the existing slider/list state too).
  const controlsContent = useMemo(
    () => <MoorhenControlPanel {...controlsProps} />,
    // Spread deps must match props the panel actually reads. Re-render
    // whenever any of these change, as the original site did.
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [
      controlsProps.onFileSelect,
      controlsProps.onJobLoad,
      controlsProps.getViewUrl,
      controlsProps.molecules,
      controlsProps.maps,
      controlsProps.onMapContourLevelChange,
      controlsProps.onRunServalcat,
      controlsProps.servalcatStatus,
    ],
  );

  const scenesContent = useMemo(
    () => (
      <MoorhenScenesPanel
        onApplyScene={onApplyScene}
        onCaptureScene={onCaptureScene}
        onPromoteSceneToPortable={onPromoteSceneToPortable}
        enabled={cootInitialized}
      />
    ),
    [onApplyScene, onCaptureScene, onPromoteSceneToPortable, cootInitialized],
  );

  return (
    <Box
      sx={{
        display: "flex",
        flexDirection: "column",
        height: "100%",
        width: "100%",
        minWidth: 0,
      }}
    >
      <Tabs
        value={activeTab}
        onChange={(_, v) => setActiveTab(v)}
        variant="fullWidth"
        sx={{ borderBottom: 1, borderColor: "divider", minHeight: 36 }}
      >
        <Tab
          label="Controls"
          value="controls"
          sx={{ minHeight: 36, py: 0.5, fontSize: "0.75rem", textTransform: "none" }}
        />
        <Tab
          label="Scenes"
          value="scenes"
          sx={{ minHeight: 36, py: 0.5, fontSize: "0.75rem", textTransform: "none" }}
        />
      </Tabs>
      <Box
        sx={{
          flex: 1,
          minHeight: 0,
          width: "100%",
          overflow: "hidden",
          display: activeTab === "controls" ? "block" : "none",
        }}
      >
        {controlsContent}
      </Box>
      <Box
        sx={{
          flex: 1,
          minHeight: 0,
          width: "100%",
          overflow: "hidden",
          display: activeTab === "scenes" ? "block" : "none",
        }}
      >
        {scenesContent}
      </Box>
    </Box>
  );
};
