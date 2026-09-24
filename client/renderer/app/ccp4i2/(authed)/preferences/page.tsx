"use client";
import { ProgramLocations } from "@/components/program-locations";
import { ProjectsDirectory } from "@/components/projects-directory";
import { CredentialsPanel } from "@/components/credentials-panel";
import { GeneralPreferencesPanel } from "@/components/general-preferences-panel";
import { Divider, Paper } from "@mui/material";
import { useTopBar } from "@/providers/top-bar-context";

/**
 * Preferences — running-app settings (distinct from the launch/get-ready
 * screen at /ccp4i2/config). Sections: Program locations (binary discovery)
 * and Credentials (tokens/passwords for external services). More preference
 * sections can be added here over time.
 */
export default function PreferencesPage() {
  useTopBar({ title: "Preferences" });
  return (
    <Paper sx={{ flex: 1, minHeight: 0, overflowY: "auto", py: 3 }}>
      <GeneralPreferencesPanel />
      <Divider sx={{ my: 2 }} />
      <ProjectsDirectory />
      <Divider sx={{ my: 2 }} />
      <ProgramLocations />
      <Divider sx={{ my: 2 }} />
      <CredentialsPanel />
    </Paper>
  );
}
