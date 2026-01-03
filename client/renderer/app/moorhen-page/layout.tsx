"use client";
import { PropsWithChildren } from "react";
import { CootProvider } from "../../providers/coot-provider";
import { AppBar, Toolbar, Typography } from "@mui/material";
import { useParams, useRouter } from "next/navigation";
import { NavigationShortcutsProvider } from "../../providers/navigation-shortcuts-provider";
import { ThemeToggle } from "../../components/theme-toggle";

export default function MoorhenPageLayout(props: PropsWithChildren) {
  const router = useRouter();
  const { id } = useParams();
  const fileIds = id ? [parseInt(id as string)] : [];
  return (
    <CootProvider>
      <NavigationShortcutsProvider>
        <AppBar position="static">
          <Toolbar sx={{ gap: 2 }}>
            <ThemeToggle />
            <Typography
              variant="h6"
              style={{ textAlign: "center", margin: "10px" }}
            >
              Moorhen Viewer for file with Id{" "}
              {fileIds.length > 0 ? fileIds[0] : "unknown"}
            </Typography>
          </Toolbar>
        </AppBar>
        {props.children}
      </NavigationShortcutsProvider>
    </CootProvider>
  );
}
