/**
 * Campaign tabbed side-panel — hosts the campaign Controls and the Scenes
 * editor under one Moorhen side-panel registration, mirroring the generic
 * job/file viewer's CCP4i2 tabbed panel.
 *
 * Registering two custom extraSidePanels with Moorhen only ever surfaces one
 * tab (an upstream quirk), so we put our own MUI Tabs inside one panel. Both
 * sub-panels stay mounted (the inactive one is display:none) so Monaco's
 * editor state in the Scenes tab survives switching back to Controls.
 *
 * The campaign viewer applies scenes but does not lift them, so the Scenes
 * panel here is given only onApplyScene — its Capture / Save-self-contained
 * actions are hidden (those props are optional on MoorhenScenesPanel).
 */
import React, { useMemo, useState } from "react";
import { Box, Tab, Tabs } from "@mui/material";

import { CampaignControlPanel } from "./campaign-control-panel";
import { MoorhenScenesPanel, SceneBundleAssets } from "./moorhen-scenes-panel";
import type { SceneResolveResult } from "../../lib/moorhen-scene-resolver";

export interface CampaignMoorhenTabbedPanelProps {
  /** Everything the campaign control panel needs, forwarded verbatim. */
  controlPanelProps: React.ComponentProps<typeof CampaignControlPanel>;
  /** Apply the editor's YAML (+ bundle assets) to the current view. */
  onApplyScene: (
    yamlText: string,
    assets: SceneBundleAssets,
  ) => Promise<SceneResolveResult>;
  /** True once Coot is ready (gates Apply / live-apply). */
  cootInitialized: boolean;
  /** Optional YAML to seed the Scenes editor with (the summary scene). */
  initialSceneYaml?: string;
  /** Apply `initialSceneYaml` automatically through the Scenes panel. */
  autoApplyInitialScene?: boolean;
}

export const CampaignMoorhenTabbedPanel: React.FC<CampaignMoorhenTabbedPanelProps> = ({
  controlPanelProps,
  onApplyScene,
  cootInitialized,
  initialSceneYaml,
  autoApplyInitialScene,
}) => {
  const [activeTab, setActiveTab] = useState<"campaign" | "scenes">("campaign");

  const controlsContent = useMemo(
    () => <CampaignControlPanel {...controlPanelProps} />,
    [controlPanelProps],
  );

  // Stable across control-panel re-renders so Monaco editor state survives.
  const scenesContent = useMemo(
    () => (
      <MoorhenScenesPanel
        onApplyScene={onApplyScene}
        enabled={cootInitialized}
        initialYaml={initialSceneYaml}
        autoApplyInitial={autoApplyInitialScene}
      />
    ),
    [onApplyScene, cootInitialized, initialSceneYaml, autoApplyInitialScene],
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
          label="Campaign"
          value="campaign"
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
          display: activeTab === "campaign" ? "block" : "none",
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
