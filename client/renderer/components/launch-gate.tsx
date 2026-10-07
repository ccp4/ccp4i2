"use client";
import React, { useEffect, useState } from "react";
import {
  Box,
  Button,
  CircularProgress,
  Fade,
  Stack,
  Typography,
} from "@mui/material";
import { ErrorOutline } from "@mui/icons-material";
import { CCP4Icon } from "./General/CCP4i2Icons";
import { useServerReady } from "@/hooks/use-server-ready";

/**
 * Whether the backend has answered once already in this window's session.
 * The gate wraps the projects page, so without this every return to that
 * page (a client-side navigation, which remounts the gate) showed
 * "Starting CCP4i2..." again, as if the app were relaunching (#732). Module
 * state lives as long as the renderer does and so resets with the app; it is
 * deliberately not sessionStorage, which the server render cannot read and
 * would make a reload's first render disagree with the server's.
 */
let serverReadyThisSession = false;

/**
 * Waits for the backend to be reachable before rendering its children.
 *
 * This closes the gap where the app used to redirect to the projects list the
 * instant the server process was spawned — before it could answer — causing a
 * first-launch flicker or a transient "couldn't load projects". Instead we show
 * a calm "Starting CCP4i2…" state and only reveal the app once /health responds.
 */
export const LaunchGate: React.FC<{ children: React.ReactNode }> = ({
  children,
}) => {
  const [alreadyReady] = useState(() => serverReadyThisSession);
  const { status, attempts, retry } = useServerReady({
    intervalMs: 1000,
    timeoutMs: 60_000,
    enabled: !alreadyReady,
  });

  useEffect(() => {
    if (status === "ready") serverReadyThisSession = true;
  }, [status]);

  if (alreadyReady || status === "ready") {
    // A flex column, not a plain block. Every page inside this gate starts
    // its own height chain with flex: 1, and a flex item needs a flex
    // container to be one -- against a block parent it is inert, each child
    // sizes to its content, and the chain that was meant to end in a
    // scrollable pane ends in an element as tall as its contents. The project
    // list showed it: MuiTableContainer-root computed 1252 x 8409px with
    // flex: 1 1 0% already set, so the virtualiser measured a viewport
    // thousands of pixels tall, windowed nothing, and the ancestors'
    // overflow: hidden clipped the result -- rows to the bottom of the window
    // and no scrollbar anywhere.
    return (
      <Fade in>
        {
          <Box
            sx={{
              height: "100%",
              minHeight: 0,
              display: "flex",
              flexDirection: "column",
            }}
          >
            {children}
          </Box>
        }
      </Fade>
    );
  }

  return (
    <Stack
      alignItems="center"
      justifyContent="center"
      spacing={3}
      sx={{ height: "100vh", px: 3, textAlign: "center" }}
    >
      {status === "checking" ? (
        <>
          <Box sx={{ position: "relative", display: "grid", placeItems: "center" }}>
            <CircularProgress size={72} thickness={2.5} />
            <CCP4Icon
              sx={{
                position: "absolute",
                fontSize: 32,
                color: "primary.main",
              }}
            />
          </Box>
          <Stack spacing={0.5}>
            <Typography variant="h5" fontWeight={600}>
              Starting CCP4i2…
            </Typography>
            <Typography variant="body2" color="text.secondary">
              Bringing the crystallographic backend online
            </Typography>
          </Stack>
          {attempts > 5 && (
            <Typography variant="caption" color="text.secondary">
              Still starting — this can take a moment on first launch.
            </Typography>
          )}
        </>
      ) : (
        <>
          <ErrorOutline sx={{ fontSize: 56, color: "warning.main" }} />
          <Stack spacing={0.5}>
            <Typography variant="h5" fontWeight={600}>
              CCP4i2 didn&apos;t start
            </Typography>
            <Typography variant="body2" color="text.secondary" sx={{ maxWidth: 420 }}>
              The backend server hasn&apos;t responded. It may still be installing,
              or the setup may need attention.
            </Typography>
          </Stack>
          <Button variant="contained" onClick={retry}>
            Try again
          </Button>
        </>
      )}
    </Stack>
  );
};
