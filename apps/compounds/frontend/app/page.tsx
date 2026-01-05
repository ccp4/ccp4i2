'use client';

import { Container, Typography, Box, Button, Stack, Paper } from '@mui/material';
import { Science, Assignment, Biotech } from '@mui/icons-material';
import Link from 'next/link';

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
          href="/registry/targets"
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
            opacity: 0.6,
          }}
        >
          <Biotech sx={{ fontSize: 48, color: 'secondary.main' }} />
          <Box>
            <Typography variant="h5">Assays</Typography>
            <Typography color="text.secondary">
              Coming soon - protocols, experiments, and analysis
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
            opacity: 0.6,
          }}
        >
          <Assignment sx={{ fontSize: 48, color: 'warning.main' }} />
          <Box>
            <Typography variant="h5">Hypotheses</Typography>
            <Typography color="text.secondary">
              Coming soon - compound design tracking
            </Typography>
          </Box>
        </Paper>
      </Stack>
    </Container>
  );
}
