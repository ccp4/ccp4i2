"use client";
import { useEffect } from "react";
import { useRouter } from "next/navigation";
import { Box, CircularProgress, Typography } from "@mui/material";

/**
 * Root landing page - redirects to ccp4i2 app
 * 
 * In the future, this could be an app selector showing:
 * - CCP4i2 (crystallography workbench)
 * - Compounds (ligand database)
 * - etc.
 */
export default function RootPage() {
  const router = useRouter();

  useEffect(() => {
    // For now, redirect directly to ccp4i2
    // TODO: When compounds is added, show an app selector instead
    router.replace("/ccp4i2");
  }, [router]);

  return (
    <Box
      sx={{
        display: "flex",
        flexDirection: "column",
        alignItems: "center",
        justifyContent: "center",
        minHeight: "100vh",
        gap: 2,
      }}
    >
      <CircularProgress />
      <Typography variant="body2" color="text.secondary">
        Loading CCP4...
      </Typography>
    </Box>
  );
}
