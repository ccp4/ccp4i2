'use client';

import { Container, Typography, Box, Stack, Paper, Tooltip } from '@mui/material';
import { Science, Biotech, TableChart, Search, AccountTree, QuestionAnswer, DoNotDisturb } from '@mui/icons-material';
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
            position: 'relative',
            '&:hover': { bgcolor: 'action.hover' },
          }}
          component={Link}
          href="/nlp"
        >
          <QuestionAnswer sx={{ fontSize: 56, color: 'primary.main' }} />
          <Box>
            <Typography variant="h5">Ask (Natural-Language Query)</Typography>
            <Typography color="text.secondary">
              Describe the compounds you want and land on the aggregation page
            </Typography>
          </Box>
          <Tooltip
            title="Beta feature — may misinterpret prompts or return unexpected selections. Check the scope sentence before following the redirect."
            arrow
          >
            <DoNotDisturb
              sx={{
                position: 'absolute',
                top: 8,
                right: 8,
                fontSize: 40,
                color: 'error.main',
                opacity: 0.6,
              }}
            />
          </Tooltip>
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
