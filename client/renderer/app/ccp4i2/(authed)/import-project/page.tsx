"use client";
import { ImportProjectContent } from "@/components/import-project-content";
import { Box, Stack } from "@mui/material";
import CCP4i2TopBar from "@/components/ccp4i2-topbar";

export default function ImportProjectPage() {
  // Viewport-high flex column: the top bar is an ordinary flex child, so the
  // content area is whatever is left, however short the window and whatever
  // the bar's height at this breakpoint. (A hardcoded `calc(100vh - 80px)`
  // used to push the card up under the bar on short windows.) The content
  // area scrolls; the card centres itself with `margin: auto`, which — unlike
  // alignItems: center — collapses to 0 when the card is taller than the
  // area, so its top is never clipped.
  return (
    <Stack
      sx={{
        height: "100vh",
        "@supports (height: 100dvh)": { height: "100dvh" },
        overflow: "hidden",
      }}
    >
      <CCP4i2TopBar title="Import Project" showBackButton backPath="/ccp4i2" />
      <Box sx={{ flex: 1, minHeight: 0, overflow: "auto", display: "flex" }}>
        <ImportProjectContent />
      </Box>
    </Stack>
  );
}
