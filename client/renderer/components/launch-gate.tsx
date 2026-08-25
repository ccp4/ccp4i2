"use client";
import React from "react";
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
  const { status, attempts, retry } = useServerReady({
    intervalMs: 1000,
    timeoutMs: 60_000,
  });

  if (status === "ready") {
    return <Fade in>{<Box sx={{ height: "100%" }}>{children}</Box>}</Fade>;
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
