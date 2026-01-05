'use client';

import { useRouter } from 'next/navigation';
import { Container, Typography, Box, Chip } from '@mui/material';
import { Science } from '@mui/icons-material';
import { Breadcrumbs } from '@/components/Breadcrumbs';
import { DataTable, Column } from '@/components/DataTable';
import { useCompoundsApi } from '@/lib/api';
import { Target } from '@/types/models';

export default function TargetsPage() {
  const router = useRouter();
  const api = useCompoundsApi();
  const { data: targets, isLoading } = api.get<Target[]>('targets/');

  const columns: Column<Target>[] = [
    {
      key: 'name',
      label: 'Target Name',
      sortable: true,
      searchable: true,
      render: (value) => (
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Science fontSize="small" color="primary" />
          <Typography fontWeight={500}>{value}</Typography>
        </Box>
      ),
    },
    {
      key: 'compound_count',
      label: 'Compounds',
      sortable: true,
      width: 120,
      render: (value) =>
        value ? (
          <Chip label={value} size="small" color="primary" variant="outlined" />
        ) : (
          '-'
        ),
    },
    {
      key: 'created_at',
      label: 'Created',
      sortable: true,
      width: 120,
      render: (value) =>
        value ? new Date(value).toLocaleDateString() : '-',
    },
  ];

  return (
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <Breadcrumbs
        items={[
          { label: 'Home', href: '/', icon: 'home' },
          { label: 'Registry', href: '/registry/targets' },
          { label: 'Targets', icon: 'target' },
        ]}
      />

      <Box sx={{ mb: 3 }}>
        <Typography variant="h4" gutterBottom>
          Targets
        </Typography>
        <Typography color="text.secondary">
          Drug discovery targets and campaigns
        </Typography>
      </Box>

      <DataTable
        data={targets}
        columns={columns}
        loading={isLoading}
        onRowClick={(target) => router.push(`/registry/targets/${target.id}`)}
        getRowKey={(row) => row.id}
        title={targets ? `${targets.length} targets` : undefined}
        emptyMessage="No targets found"
      />
    </Container>
  );
}
