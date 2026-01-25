'use client';

import { Container, Typography, Box, Stack, Paper } from '@mui/material';
import { Science, Biotech, TableChart, Search, AccountTree, AdminPanelSettings, GridView, Visibility } from '@mui/icons-material';
import Link from 'next/link';

/**
 * App selector page for the web deployment.
 * This page is overlaid onto the ccp4i2 client at build time
 * to provide navigation between available applications.
 */
export default function AppSelectorPage() {
  return (
    <Container maxWidth="md" sx={{ py: 6 }}>
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
        </Box>
      </Box>
    </Container>
  );
}
