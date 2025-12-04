import { PropsWithChildren } from "react";
import { ImportProjectContent } from "../../components/import-project-content";
import { Paper } from "@mui/material";
import { NavigationShortcutsProvider } from "../../providers/navigation-shortcuts-provider";

export default function ImportProjectPage() {
  return (
    <NavigationShortcutsProvider>
      <Paper
        sx={{
          display: "flex",
          justifyContent: "center",
          alignItems: "center",
          height: "100vh",
        }}
      >
        <ImportProjectContent />
      </Paper>
    </NavigationShortcutsProvider>
  );
}
