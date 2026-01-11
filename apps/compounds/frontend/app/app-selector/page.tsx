'use client';

import { Container, Typography, Box, Stack, Paper } from '@mui/material';
import { Science, Biotech } from '@mui/icons-material';
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
          ccp4i2
        </Typography>
        <Typography variant="subtitle1" color="text.secondary">
          Select an application to get started
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
          href="/registry/targets"
        >
          <Biotech sx={{ fontSize: 56, color: 'secondary.main' }} />
          <Box>
            <Typography variant="h5">Compounds Registry</Typography>
            <Typography color="text.secondary">
              Compound registration, batch tracking, and assay management
            </Typography>
          </Box>
        </Paper>
      </Stack>
    </Container>
  );
}
