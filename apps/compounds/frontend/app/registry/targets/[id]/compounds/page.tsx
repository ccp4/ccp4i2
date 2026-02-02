'use client';

import { use } from 'react';
import { useRouter } from 'next/navigation';
import {
  Container,
  Typography,
  Box,
  Paper,
  Chip,
  Skeleton,
  Button,
} from '@mui/material';
import { Medication, Science, TableChart, Add, Dashboard } from '@mui/icons-material';
import Link from 'next/link';
import { PageHeader } from '@/components/compounds/PageHeader';
import { DataTable, Column } from '@/components/data-table';
import { MoleculeChip, CopyableSmiles } from '@/components/compounds/MoleculeView';
import { useCompoundsApi } from '@/lib/compounds/api';
import { routes } from '@/lib/compounds/routes';
import { useAdaptiveTableHeight } from '@/lib/compounds/useAdaptiveTableHeight';
import { Target, Compound } from '@/types/compounds/models';

interface PageProps {
  params: Promise<{ id: string }>;
}

export default function TargetCompoundsPage({ params }: PageProps) {
  const { id } = use(params);
  const router = useRouter();
  const api = useCompoundsApi();

  const { data: target, isLoading: targetLoading } = api.get<Target>(
    `targets/${id}/`
  );
  const { data: compounds, isLoading: compoundsLoading } = api.get<Compound[]>(
    `compounds/?target=${id}`
  );

  // Adaptive table height for compounds table
  const { headerRef, tableHeight } = useAdaptiveTableHeight({
    maxHeight: 600,
    minHeight: 400,
    bottomPadding: 32,
  });

  const columns: Column<Compound>[] = [
    {
      key: 'formatted_id',
      label: 'ID',
      sortable: true,
      searchable: true,
      width: 180,
      render: (value, row) => (
        <Box sx={{ display: 'flex', alignItems: 'flex-start', gap: 1 }}>
          <Medication fontSize="small" color="secondary" sx={{ mt: 0.3 }} />
          <Box>
            <Typography fontWeight={500} fontFamily="monospace">
              {value}
            </Typography>
            {row.supplier_ref && (
              <Typography
                variant="caption"
                color="text.secondary"
                sx={{ display: 'block', fontFamily: 'monospace' }}
              >
                {row.supplier_ref}
              </Typography>
            )}
          </Box>
        </Box>
      ),
    },
    {
      key: 'smiles',
      label: 'Structure',
      searchable: true,
      width: 175,
      render: (value) => <MoleculeChip smiles={value} size={140} />,
    },
    {
      key: 'smiles',
      label: 'SMILES',
      searchable: true,
      render: (value) => <CopyableSmiles smiles={value} />,
    },
    {
      key: 'molecular_weight',
      label: 'MW',
      sortable: true,
      width: 80,
      render: (value) => (value ? value.toFixed(1) : '-'),
    },
    {
      key: 'stereo_comment',
      label: 'Stereo',
      sortable: true,
      width: 100,
      render: (value) =>
        value && value !== 'unset' ? (
          <Chip label={value} size="small" variant="outlined" />
        ) : (
          '-'
        ),
    },
    {
      key: 'batch_count',
      label: 'Batches',
      sortable: true,
      width: 80,
      render: (value) =>
        value ? (
          <Chip label={value} size="small" color="info" variant="outlined" />
        ) : (
          '-'
        ),
    },
    {
      key: 'registered_at',
      label: 'Registered',
      sortable: true,
      width: 100,
      render: (value) =>
        value ? new Date(value).toLocaleDateString() : '-',
    },
  ];

  return (
    <Box sx={{ display: 'flex', flexDirection: 'column', height: '100vh', overflow: 'hidden' }}>
      <Container maxWidth="lg" sx={{ flex: 1, display: 'flex', flexDirection: 'column', overflow: 'hidden', py: 2 }}>
        {/* Scrollable header - inset shadow at bottom hints at scrollable content */}
        <Box
          ref={headerRef}
          sx={{
            flexShrink: 0,
            overflow: 'auto',
            maxHeight: '50vh',
            // Subtle inset shadow at bottom edge to hint at scroll
            boxShadow: 'inset 0 -12px 12px -12px rgba(0,0,0,0.08)',
          }}
        >
          <PageHeader
            breadcrumbs={[
              { label: 'Home', href: routes.home(), icon: 'home' },
              { label: 'Targets', href: routes.registry.targets(), icon: 'target' },
              { label: target?.name || 'Loading...', href: routes.registry.target(id), icon: 'target' },
              { label: 'Compounds', icon: 'compound' },
            ]}
          />

          {/* Target header */}
          <Paper sx={{ p: 3 }}>
            {targetLoading ? (
              <>
                <Skeleton variant="text" width={300} height={40} />
                <Skeleton variant="text" width={200} />
              </>
            ) : target ? (
              <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
                <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
                  <Science sx={{ fontSize: 48, color: 'primary.main' }} />
                  <Box>
                    <Typography variant="h4">{target.name}</Typography>
                    <Typography color="text.secondary">
                      {compounds?.length ?? 0} compounds registered
                    </Typography>
                  </Box>
                </Box>
                <Box sx={{ display: 'flex', gap: 1 }}>
                  <Button
                    component={Link}
                    href={routes.registry.target(id)}
                    variant="outlined"
                    startIcon={<Dashboard />}
                  >
                    Dashboard
                  </Button>
                  <Button
                    component={Link}
                    href={routes.assays.aggregate({ target: id })}
                    variant="outlined"
                    startIcon={<TableChart />}
                  >
                    View Assay Data
                  </Button>
                  <Button
                    component={Link}
                    href={`${routes.registry.new()}?target=${id}`}
                    variant="contained"
                    startIcon={<Add />}
                  >
                    New Compound
                  </Button>
                </Box>
              </Box>
            ) : (
              <Typography color="error">Target not found</Typography>
            )}
          </Paper>
        </Box>

        {/* Compounds table - adaptive height */}
        <Box sx={{ flex: 1, overflow: 'hidden', minHeight: 0, height: tableHeight, mt: 2 }}>
          <DataTable
            data={compounds}
            columns={columns}
            loading={compoundsLoading}
            onRowClick={(compound) =>
              router.push(routes.registry.compound(compound.id))
            }
            getRowKey={(row) => row.id}
            title={compounds ? `${compounds.length} compounds` : undefined}
            emptyMessage="No compounds registered for this target"
            additionalSearchFields={['supplier_ref', 'supplier_name', 'barcode', 'comments']}
            fillHeight
          />
        </Box>
      </Container>
    </Box>
  );
}
