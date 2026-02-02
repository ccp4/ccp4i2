'use client';

import { Container, Typography, Box, Stack, Paper, Collapse } from '@mui/material';
import { Science, Biotech, TableChart, Search, AccountTree, AdminPanelSettings, GridView, Visibility } from '@mui/icons-material';
import Link from 'next/link';
import { useRouter, useSearchParams } from 'next/navigation';
import { useState, useEffect, useRef, Suspense } from 'react';

interface VersionInfo {
  web?: {
    buildTimestamp?: string;
    gitCommit?: string;
  };
  server?: {
    buildTimestamp?: string;
    gitCommit?: string;
  };
}

const TEAMS_ROUTE_KEY = "ccp4i2-teams-last-route";

/**
 * Check if running in an iframe (Teams context).
 */
function isRunningInIframe(): boolean {
  if (typeof window === "undefined") return false;
  try {
    return window.self !== window.top;
  } catch {
    return true;
  }
}

/**
 * Clear the saved route from sessionStorage.
 */
function clearSavedTeamsRoute(): void {
  if (typeof window === "undefined") return;
  try {
    sessionStorage.removeItem(TEAMS_ROUTE_KEY);
  } catch {
    // Ignore
  }
}

/**
 * Get and clear the saved route from sessionStorage.
 * Returns null if not in Teams context or no route saved.
 */
function getAndClearSavedTeamsRoute(): string | null {
  if (typeof window === "undefined") return null;
  if (!isRunningInIframe()) return null;

  try {
    const route = sessionStorage.getItem(TEAMS_ROUTE_KEY);
    if (route) {
      sessionStorage.removeItem(TEAMS_ROUTE_KEY);
    }
    return route;
  } catch {
    return null;
  }
}

/**
 * Component that handles Teams route restoration.
 * Uses useSearchParams which requires Suspense boundary.
 */
function TeamsRouteHandler() {
  const router = useRouter();
  const searchParams = useSearchParams();
  const hasCheckedRoute = useRef(false);

  useEffect(() => {
    // Check for ?home parameter - this bypasses route restoration
    // Useful if saved route is broken and user needs to get back to app-selector
    if (searchParams.has('home')) {
      console.log("[APP-SELECTOR] ?home parameter detected, clearing saved route");
      clearSavedTeamsRoute();
      // Don't redirect, stay on app-selector
    } else if (!hasCheckedRoute.current) {
      // Check for saved route (Teams context only) - do this first
      hasCheckedRoute.current = true;
      const savedRoute = getAndClearSavedTeamsRoute();
      if (savedRoute && savedRoute !== "/") {
        console.log("[APP-SELECTOR] Restoring saved Teams route:", savedRoute);
        router.replace(savedRoute);
      }
    }
  }, [router, searchParams]);

  return null; // This component only handles side effects
}

/**
 * App selector page for the web deployment.
 * This page is overlaid onto the ccp4i2 client at build time
 * to provide navigation between available applications.
 *
 * In Teams context, if there's a saved route from a previous session,
 * this page will automatically redirect to that route.
 *
 * Use /?home to bypass route restoration and stay on this page
 * (useful if the saved route leads to a broken page).
 */
export default function AppSelectorPage() {
  const [versionInfo, setVersionInfo] = useState<VersionInfo>({});
  const [showVersion, setShowVersion] = useState(false);

  useEffect(() => {
    // Fetch web version
    fetch('/api/version')
      .then(res => res.json())
      .then(data => setVersionInfo(prev => ({ ...prev, web: data.web })))
      .catch(() => {});

    // Fetch server version
    fetch('/api/proxy/ccp4i2/version/')
      .then(res => res.json())
      .then(data => setVersionInfo(prev => ({ ...prev, server: data })))
      .catch(() => {});
  }, []);
  return (
    <Container maxWidth="md" sx={{ py: 6 }}>
      {/* Handle Teams route restoration - wrapped in Suspense for useSearchParams */}
      <Suspense fallback={null}>
        <TeamsRouteHandler />
      </Suspense>

      <Box sx={{ textAlign: 'center', mb: 5 }}>
        <Typography variant="h3" component="h1" gutterBottom>
          SBDD Database
        </Typography>
        <Typography variant="subtitle1" color="text.secondary">
          Structure Based Drug Design - Select an application to get started
        </Typography>
      </Box>

      <Stack spacing={3}>
        <Paper
          elevation={2}
          sx={{
            p: 3,
            display: 'flex',
            alignItems: 'center',
            gap: 3,
            cursor: 'pointer',
            textDecoration: 'none',
            color: 'inherit',
            '&:hover': { bgcolor: 'action.hover' },
          }}
          component={Link}
          href="/registry/targets"
        >
          <Science sx={{ fontSize: 56, color: 'secondary.main' }} />
          <Box>
            <Typography variant="h5">Drug Discovery Targets</Typography>
            <Typography color="text.secondary">
              Target dashboards with compounds, assays, and related projects
            </Typography>
          </Box>
        </Paper>

        <Paper
          elevation={2}
          sx={{
            p: 3,
            display: 'flex',
            alignItems: 'center',
            gap: 3,
            cursor: 'pointer',
            textDecoration: 'none',
            color: 'inherit',
            '&:hover': { bgcolor: 'action.hover' },
          }}
          component={Link}
          href="/ccp4i2"
        >
          <Science sx={{ fontSize: 56, color: 'primary.main' }} />
          <Box>
            <Typography variant="h5">Crystallography Workbench</Typography>
            <Typography color="text.secondary">
              Structure determination, refinement, and analysis pipelines
            </Typography>
          </Box>
        </Paper>

        <Paper
          elevation={2}
          sx={{
            p: 3,
            display: 'flex',
            alignItems: 'center',
            gap: 3,
            cursor: 'pointer',
            textDecoration: 'none',
            color: 'inherit',
            '&:hover': { bgcolor: 'action.hover' },
          }}
          component={Link}
          href="/ccp4i2/campaigns"
        >
          <GridView sx={{ fontSize: 56, color: 'success.main' }} />
          <Box>
            <Typography variant="h5">Fragment Screening</Typography>
            <Typography color="text.secondary">
              Manage fragment screening campaigns with batch processing
            </Typography>
          </Box>
        </Paper>

        <Paper
          elevation={2}
          sx={{
            p: 3,
            display: 'flex',
            alignItems: 'center',
            gap: 3,
            cursor: 'pointer',
            textDecoration: 'none',
            color: 'inherit',
            '&:hover': { bgcolor: 'action.hover' },
          }}
          component="a"
          href="/ccp4i2/moorhen-page"
          target="_blank"
          rel="noopener"
        >
          <Visibility sx={{ fontSize: 56, color: 'primary.light' }} />
          <Box>
            <Typography variant="h5">Moorhen Viewer</Typography>
            <Typography color="text.secondary">
              Opens in new tab (requires cross-origin isolation)
            </Typography>
          </Box>
        </Paper>

        <Paper
          elevation={2}
          sx={{
            p: 3,
            display: 'flex',
            alignItems: 'center',
            gap: 3,
            cursor: 'pointer',
            textDecoration: 'none',
            color: 'inherit',
            '&:hover': { bgcolor: 'action.hover' },
          }}
          component={Link}
          href="/registry/compounds"
        >
          <Biotech sx={{ fontSize: 56, color: 'info.main' }} />
          <Box>
            <Typography variant="h5">Compounds Registry</Typography>
            <Typography color="text.secondary">
              Browse all compounds, filter by target, and register new compounds
            </Typography>
          </Box>
        </Paper>

        <Paper
          elevation={2}
          sx={{
            p: 3,
            display: 'flex',
            alignItems: 'center',
            gap: 3,
            cursor: 'pointer',
            textDecoration: 'none',
            color: 'inherit',
            '&:hover': { bgcolor: 'action.hover' },
          }}
          component={Link}
          href="/registry/search"
        >
          <Search sx={{ fontSize: 56, color: 'warning.main' }} />
          <Box>
            <Typography variant="h5">Compound Search</Typography>
            <Typography color="text.secondary">
              Search compounds by ID, supplier reference, or structure
            </Typography>
          </Box>
        </Paper>

        <Paper
          elevation={2}
          sx={{
            p: 3,
            display: 'flex',
            alignItems: 'center',
            gap: 3,
            cursor: 'pointer',
            textDecoration: 'none',
            color: 'inherit',
            '&:hover': { bgcolor: 'action.hover' },
          }}
          component={Link}
          href="/assays/protocols"
        >
          <Biotech sx={{ fontSize: 56, color: 'info.main' }} />
          <Box>
            <Typography variant="h5">Assays</Typography>
            <Typography color="text.secondary">
              Protocols, experiments, and dose-response analysis
            </Typography>
          </Box>
        </Paper>

        <Paper
          elevation={2}
          sx={{
            p: 3,
            display: 'flex',
            alignItems: 'center',
            gap: 3,
            cursor: 'pointer',
            textDecoration: 'none',
            color: 'inherit',
            '&:hover': { bgcolor: 'action.hover' },
          }}
          component={Link}
          href="/assays/aggregate"
        >
          <TableChart sx={{ fontSize: 56, color: 'success.main' }} />
          <Box>
            <Typography variant="h5">Data Aggregation</Typography>
            <Typography color="text.secondary">
              Query and aggregate KPI values across compounds and protocols
            </Typography>
          </Box>
        </Paper>

        <Paper
          elevation={2}
          sx={{
            p: 3,
            display: 'flex',
            alignItems: 'center',
            gap: 3,
            cursor: 'pointer',
            textDecoration: 'none',
            color: 'inherit',
            '&:hover': { bgcolor: 'action.hover' },
          }}
          component={Link}
          href="/constructs"
        >
          <AccountTree sx={{ fontSize: 56, color: 'error.main' }} />
          <Box>
            <Typography variant="h5">Construct Database</Typography>
            <Typography color="text.secondary">
              Plasmid registry, cassette tracking, and sequencing results
            </Typography>
          </Box>
        </Paper>

        <Paper
          elevation={2}
          sx={{
            p: 3,
            display: 'flex',
            alignItems: 'center',
            gap: 3,
            cursor: 'pointer',
            textDecoration: 'none',
            color: 'inherit',
            '&:hover': { bgcolor: 'action.hover' },
          }}
          component={Link}
          href="/admin"
        >
          <AdminPanelSettings sx={{ fontSize: 56, color: 'grey.600' }} />
          <Box>
            <Typography variant="h5">Platform Admin</Typography>
            <Typography color="text.secondary">
              User management, platform settings, and data import tools
            </Typography>
          </Box>
        </Paper>
      </Stack>

      <Box
        sx={{
          mt: 6,
          pt: 3,
          borderTop: 1,
          borderColor: 'divider',
          textAlign: 'center',
        }}
      >
        <Typography variant="body2" color="text.secondary">
          Newcastle University
        </Typography>
        <Box sx={{ mt: 1, display: 'flex', justifyContent: 'center', gap: 3 }}>
          <Typography
            component={Link}
            href="/privacy"
            variant="body2"
            color="text.secondary"
            sx={{ textDecoration: 'none', '&:hover': { textDecoration: 'underline' } }}
          >
            Privacy Policy
          </Typography>
          <Typography
            component={Link}
            href="/terms"
            variant="body2"
            color="text.secondary"
            sx={{ textDecoration: 'none', '&:hover': { textDecoration: 'underline' } }}
          >
            Terms of Use
          </Typography>
          <Typography
            component="span"
            variant="body2"
            color="text.secondary"
            sx={{ cursor: 'pointer', '&:hover': { textDecoration: 'underline' } }}
            onClick={() => setShowVersion(!showVersion)}
          >
            About
          </Typography>
        </Box>

        <Collapse in={showVersion}>
          <Box sx={{ mt: 2, p: 2, bgcolor: 'action.hover', borderRadius: 1 }}>
            <Typography variant="caption" component="div" color="text.secondary">
              <strong>Web:</strong> {versionInfo.web?.buildTimestamp || 'dev'} ({versionInfo.web?.gitCommit || 'unknown'})
            </Typography>
            <Typography variant="caption" component="div" color="text.secondary">
              <strong>Server:</strong> {versionInfo.server?.buildTimestamp || 'dev'} ({versionInfo.server?.gitCommit || 'unknown'})
            </Typography>
          </Box>
        </Collapse>
      </Box>
    </Container>
  );
}
