"use client";
import { ImportProjectContent } from "@/components/import-project-content";
import { Box} from "@mui/material";
import { NavigationShortcutsProvider } from "@/providers/navigation-shortcuts-provider";
import CCP4i2TopBar from "@/components/ccp4i2-topbar";

export default function ImportProjectPage() {
  return (
    <NavigationShortcutsProvider>
      <CCP4i2TopBar title="Import Project" showBackButton backPath="/ccp4i2" />
      <Box
        sx={{
          display: "flex",
          justifyContent: "center",
          alignItems: "center",
          height: "calc(100vh - 80px)",
        }}
      >
        <ImportProjectContent />
      </Box>
    </NavigationShortcutsProvider>
  );
}
