"use client";
import { ProgramLocations } from "@/components/program-locations";
import { CredentialsPanel } from "@/components/credentials-panel";
import { Divider, Paper, Stack } from "@mui/material";
import CCP4i2TopBar from "@/components/ccp4i2-topbar";

/**
 * Preferences — running-app settings (distinct from the launch/get-ready
 * screen at /ccp4i2/config). Sections: Program locations (binary discovery)
 * and Credentials (tokens/passwords for external services). More preference
 * sections can be added here over time.
 */
export default function PreferencesPage() {
  return (
    <Stack
      sx={{
        height: "100vh",
        "@supports (height: 100dvh)": { height: "100dvh" },
        overflow: "hidden",
      }}
    >
      <CCP4i2TopBar title="Preferences" showBackButton backPath="/ccp4i2" />
      <Paper sx={{ flex: 1, minHeight: 0, overflowY: "auto", py: 3 }}>
        <ProgramLocations />
        <Divider sx={{ my: 2 }} />
        <CredentialsPanel />
      </Paper>
    </Stack>
  );
}
