"use client";
import { Paper } from "@mui/material";
import { NavigationShortcutsProvider } from "@/providers/navigation-shortcuts-provider";
import CCP4i2TopBar from "@/components/ccp4i2-topbar";
import { MigrateLegacyContent } from "@/components/migrate-legacy-content";

export default function MigrateLegacyPage() {
  return (
    <NavigationShortcutsProvider>
      <CCP4i2TopBar
        title="Migrate legacy CCP4i2"
        showBackButton
        backPath="/ccp4i2"
      />
      <Paper
        sx={{
          height: "calc(100vh - 80px)",
          overflow: "auto",
          p: 3,
        }}
      >
        <MigrateLegacyContent />
      </Paper>
    </NavigationShortcutsProvider>
  );
}
