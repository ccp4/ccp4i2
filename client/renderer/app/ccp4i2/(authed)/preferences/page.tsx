"use client";
import { ProgramLocations } from "@/components/program-locations";
import { Paper } from "@mui/material";
import { NavigationShortcutsProvider } from "@/providers/navigation-shortcuts-provider";
import CCP4i2TopBar from "@/components/ccp4i2-topbar";

/**
 * Preferences — running-app settings (distinct from the launch/get-ready
 * screen at /ccp4i2/config). First section is Program locations (binary
 * discovery); more preference sections can be added here over time.
 */
export default function PreferencesPage() {
  return (
    <NavigationShortcutsProvider>
      <CCP4i2TopBar title="Preferences" showBackButton backPath="/ccp4i2" />
      <Paper
        sx={{
          minHeight: "calc(100vh - 80px)",
          overflowY: "auto",
          py: 3,
        }}
      >
        <ProgramLocations />
      </Paper>
    </NavigationShortcutsProvider>
  );
}
