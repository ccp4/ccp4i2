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
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <PageHeader
        breadcrumbs={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Targets', href: routes.registry.targets(), icon: 'target' },
          { label: target?.name || 'Loading...', href: routes.registry.target(id), icon: 'target' },
          { label: 'Compounds', icon: 'compound' },
        ]}
      />

      {/* Target header */}
      <Paper sx={{ p: 3, mb: 3 }}>
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

      {/* Compounds table */}
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
      />
    </Container>
  );
}
