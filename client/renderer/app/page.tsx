"use client";
import { useEffect, useState } from "react";
import { useRouter } from "next/navigation";
import Link from "next/link";
import {
  Box,
  CircularProgress,
  Typography,
  Container,
  Paper,
  Stack,
  Grid,
  Chip,
  Divider,
} from "@mui/material";
import {
  Science,
  Search,
  Assignment,
  TableChart,
  Folder,
  Logout,
  Person,
} from "@mui/icons-material";
import { useMsal } from "@azure/msal-react";
import { isElectron } from "../utils/platform";

const REQUIRE_AUTH = process.env.NEXT_PUBLIC_REQUIRE_AUTH === "true";

/**
 * Root landing page - App Selector
 *
 * In web mode: Shows app selector with CCP4i2 and Compounds
 * In Electron mode: Redirects directly to ccp4i2
 */
export default function RootPage() {
  const router = useRouter();
  const [isElectronMode, setIsElectronMode] = useState<boolean | null>(null);

  // Check if we're in Electron mode
  useEffect(() => {
    setIsElectronMode(isElectron());
  }, []);

  // In Electron mode, redirect directly to ccp4i2
  useEffect(() => {
    if (isElectronMode === true) {
      router.replace("/ccp4i2");
    }
  }, [isElectronMode, router]);

  // Show loading while determining mode
  if (isElectronMode === null || isElectronMode === true) {
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

  // Web mode - show app selector
  return <AppSelector />;
}

/**
 * App Selector component for web deployment
 */
function AppSelector() {
  const msalContext = REQUIRE_AUTH ? useMsalSafe() : null;
  const accounts = msalContext?.accounts || [];
  const instance = msalContext?.instance;
  const currentUser = accounts[0];

  const handleLogout = () => {
    if (instance) {
      instance.logoutRedirect();
    }
  };

  return (
    <Container maxWidth="lg" sx={{ py: 4 }}>
      {/* Header with user info */}
      <Box
        sx={{
          display: "flex",
          justifyContent: "space-between",
          alignItems: "center",
          mb: 4,
        }}
      >
        <Box>
          <Typography variant="h3" component="h1" gutterBottom>
            CCP4 Cloud
          </Typography>
          <Typography variant="subtitle1" color="text.secondary">
            Crystallographic computing and compound management
          </Typography>
        </Box>
        {REQUIRE_AUTH && currentUser && (
          <Paper sx={{ p: 2, display: "flex", alignItems: "center", gap: 2 }}>
            <Person color="action" />
            <Box>
              <Typography variant="body2" fontWeight={500}>
                {currentUser.name}
              </Typography>
              <Typography variant="caption" color="text.secondary">
                {currentUser.username}
              </Typography>
            </Box>
            <Chip
              icon={<Logout fontSize="small" />}
              label="Sign Out"
              size="small"
              variant="outlined"
              onClick={handleLogout}
              sx={{ cursor: "pointer" }}
            />
          </Paper>
        )}
      </Box>

      <Grid container spacing={4}>
        {/* CCP4i2 Section */}
        <Grid item xs={12} md={6}>
          <Paper sx={{ p: 3, height: "100%" }}>
            <Box sx={{ display: "flex", alignItems: "center", gap: 2, mb: 3 }}>
              <Folder sx={{ fontSize: 48, color: "primary.main" }} />
              <Box>
                <Typography variant="h4">CCP4i2</Typography>
                <Typography color="text.secondary">
                  Crystallography workbench
                </Typography>
              </Box>
            </Box>

            <Stack spacing={2}>
              <Paper
                elevation={1}
                sx={{
                  p: 2,
                  cursor: "pointer",
                  "&:hover": { bgcolor: "action.hover" },
                }}
                component={Link}
                href="/ccp4i2"
              >
                <Typography variant="h6">Projects</Typography>
                <Typography variant="body2" color="text.secondary">
                  View and manage crystallography projects
                </Typography>
              </Paper>
            </Stack>
          </Paper>
        </Grid>

        {/* Compounds Section */}
        <Grid item xs={12} md={6}>
          <Paper sx={{ p: 3, height: "100%" }}>
            <Box sx={{ display: "flex", alignItems: "center", gap: 2, mb: 3 }}>
              <Science sx={{ fontSize: 48, color: "secondary.main" }} />
              <Box>
                <Typography variant="h4">Compounds</Typography>
                <Typography color="text.secondary">
                  Ligand database and assay management
                </Typography>
              </Box>
            </Box>

            <Stack spacing={2}>
              <Paper
                elevation={1}
                sx={{
                  p: 2,
                  cursor: "pointer",
                  "&:hover": { bgcolor: "action.hover" },
                }}
                component={Link}
                href="/registry/targets"
              >
                <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
                  <Science fontSize="small" color="primary" />
                  <Typography variant="h6">Registry</Typography>
                </Box>
                <Typography variant="body2" color="text.secondary">
                  Browse targets, compounds, and batches
                </Typography>
              </Paper>

              <Paper
                elevation={1}
                sx={{
                  p: 2,
                  cursor: "pointer",
                  "&:hover": { bgcolor: "action.hover" },
                }}
                component={Link}
                href="/registry/search"
              >
                <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
                  <Search fontSize="small" color="warning" />
                  <Typography variant="h6">Compound Search</Typography>
                </Box>
                <Typography variant="body2" color="text.secondary">
                  Search by ID, supplier reference, or structure
                </Typography>
              </Paper>

              <Divider />

              <Paper
                elevation={1}
                sx={{
                  p: 2,
                  cursor: "pointer",
                  "&:hover": { bgcolor: "action.hover" },
                }}
                component={Link}
                href="/assays/protocols"
              >
                <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
                  <Assignment fontSize="small" color="info" />
                  <Typography variant="h6">Protocols</Typography>
                </Box>
                <Typography variant="body2" color="text.secondary">
                  Assay protocol definitions and methods
                </Typography>
              </Paper>

              <Paper
                elevation={1}
                sx={{
                  p: 2,
                  cursor: "pointer",
                  "&:hover": { bgcolor: "action.hover" },
                }}
                component={Link}
                href="/assays/aggregate"
              >
                <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
                  <TableChart fontSize="small" color="success" />
                  <Typography variant="h6">Data Aggregation</Typography>
                </Box>
                <Typography variant="body2" color="text.secondary">
                  Query and aggregate KPI values across compounds
                </Typography>
              </Paper>
            </Stack>
          </Paper>
        </Grid>
      </Grid>
    </Container>
  );
}

/**
 * Safe wrapper for useMsal that handles the case when MsalProvider is not present.
 */
function useMsalSafe() {
  try {
    return useMsal();
  } catch {
    return null;
  }
}
