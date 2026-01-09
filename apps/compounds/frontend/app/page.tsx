'use client';

import { Container, Typography, Box, Button, Stack, Paper } from '@mui/material';
import { Science, Assignment, Biotech, TableChart, Search } from '@mui/icons-material';
import Link from 'next/link';
import { routes } from '@/lib/routes';

export default function HomePage() {
  return (
    <Container maxWidth="md" sx={{ py: 4 }}>
      <Box sx={{ textAlign: 'center', mb: 4 }}>
        <Typography variant="h3" component="h1" gutterBottom>
          Compounds Management
        </Typography>
        <Typography variant="subtitle1" color="text.secondary">
          Compound registration, batch tracking, and assay management
        </Typography>
      </Box>

      <Stack spacing={3}>
        <Paper
          elevation={2}
          sx={{
            p: 3,
            display: 'flex',
            alignItems: 'center',
            gap: 2,
            cursor: 'pointer',
            '&:hover': { bgcolor: 'action.hover' },
          }}
          component={Link}
          href={routes.registry.targets()}
        >
          <Science sx={{ fontSize: 48, color: 'primary.main' }} />
          <Box>
            <Typography variant="h5">Registry</Typography>
            <Typography color="text.secondary">
              Browse targets, compounds, and batches
            </Typography>
          </Box>
        </Paper>

        <Paper
          elevation={2}
          sx={{
            p: 3,
            display: 'flex',
            alignItems: 'center',
            gap: 2,
            cursor: 'pointer',
            '&:hover': { bgcolor: 'action.hover' },
          }}
          component={Link}
          href={routes.registry.search()}
        >
          <Search sx={{ fontSize: 48, color: 'warning.main' }} />
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
            gap: 2,
            cursor: 'pointer',
            '&:hover': { bgcolor: 'action.hover' },
          }}
          component={Link}
          href={routes.assays.list()}
        >
          <Biotech sx={{ fontSize: 48, color: 'secondary.main' }} />
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
            gap: 2,
            cursor: 'pointer',
            '&:hover': { bgcolor: 'action.hover' },
          }}
          component={Link}
          href={routes.assays.protocols()}
        >
          <Assignment sx={{ fontSize: 48, color: 'info.main' }} />
          <Box>
            <Typography variant="h5">Protocols</Typography>
            <Typography color="text.secondary">
              Assay protocol definitions and methods
            </Typography>
          </Box>
        </Paper>

        <Paper
          elevation={2}
          sx={{
            p: 3,
            display: 'flex',
            alignItems: 'center',
            gap: 2,
            cursor: 'pointer',
            '&:hover': { bgcolor: 'action.hover' },
          }}
          component={Link}
          href={routes.assays.aggregate()}
        >
          <TableChart sx={{ fontSize: 48, color: 'success.main' }} />
          <Box>
            <Typography variant="h5">Data Aggregation</Typography>
            <Typography color="text.secondary">
              Query and aggregate KPI values across compounds and protocols
            </Typography>
          </Box>
        </Paper>

      </Stack>
    </Container>
  );
}
