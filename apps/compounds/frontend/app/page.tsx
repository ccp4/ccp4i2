'use client';

import { Container, Typography, Box, Stack, Paper } from '@mui/material';
import { Science, Biotech, TableChart, Search, AccountTree } from '@mui/icons-material';
import Link from 'next/link';

/**
 * Root page for the compounds frontend.
 *
 * In standalone development mode, this shows the app selector.
 * In Docker deployment, this file is replaced by app-selector/page.tsx.
 */
export default function HomePage() {
  return (
    <Container maxWidth="md" sx={{ py: 6 }}>
      <Box sx={{ textAlign: 'center', mb: 5 }}>
        <Typography variant="h3" component="h1" gutterBottom>
          Compounds Apps
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
          href="/registry/targets"
        >
          <Biotech sx={{ fontSize: 56, color: 'secondary.main' }} />
          <Box>
            <Typography variant="h5">Compounds Registry</Typography>
            <Typography color="text.secondary">
              Compound registration, batch tracking, and target management
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
          href="/assays"
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
      </Stack>
    </Container>
  );
}
