"use client";
import React from "react";
import { ConfigContent } from "../../../components/config-content";
import { NavigationShortcutsProvider } from "../../../providers/navigation-shortcuts-provider";
import { Box, Container, Stack, Typography } from "@mui/material";
import { CCP4Icon } from "../../../components/General/CCP4i2Icons";

/**
 * The launch / "get ready" screen — the first thing an Electron user sees.
 *
 * Deliberately NOT framed as an installer or a "system diagnostics" console.
 * On the desktop, CCP4i2 is a single-user tool, so there is intentionally no
 * teams / authentication / accounts chrome here — those belong only to a web
 * backend with real multi-user auth. This page's whole job is to get the user
 * into their work: confirm the setup is ready (or fix what isn't) and launch.
 */
export default function ConfigPage() {
  return (
    <NavigationShortcutsProvider>
      <Container maxWidth="sm" sx={{ py: { xs: 4, sm: 8 } }}>
        <Stack spacing={4}>
          <Stack spacing={1.5} alignItems="center" sx={{ textAlign: "center" }}>
            <Box
              sx={{
                width: 64,
                height: 64,
                borderRadius: 2,
                display: "grid",
                placeItems: "center",
                color: "primary.contrastText",
                background:
                  "linear-gradient(135deg, #667eea 0%, #764ba2 100%)",
              }}
            >
              <CCP4Icon sx={{ fontSize: 34 }} />
            </Box>
            <Typography variant="h4" fontWeight={700}>
              CCP4i2
            </Typography>
            <Typography variant="body1" color="text.secondary">
              Crystallographic computing on your desktop
            </Typography>
          </Stack>

          <ConfigContent />
        </Stack>
      </Container>
    </NavigationShortcutsProvider>
  );
}
