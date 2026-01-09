'use client';

import { useRouter } from 'next/navigation';
import { Container, Typography, Box, Chip } from '@mui/material';
import { Science, Description } from '@mui/icons-material';
import { Breadcrumbs } from '@/components/Breadcrumbs';
import { DataTable, Column } from '@/components/DataTable';
import { useCompoundsApi } from '@/lib/api';
import { routes } from '@/lib/routes';
import { Protocol } from '@/types/models';

const ANALYSIS_METHOD_LABELS: Record<string, string> = {
  hill_langmuir: 'Hill-Langmuir',
  hill_langmuir_fix_hill: 'Hill-Langmuir (fixed Hill)',
  hill_langmuir_fix_hill_minmax: 'Hill-Langmuir (fixed Hill/min/max)',
  hill_langmuir_fix_minmax: 'Hill-Langmuir (fixed min/max)',
  ms_intact: 'MS-Intact',
  table_of_values: 'Table of values',
};

export default function ProtocolsPage() {
  const router = useRouter();
  const api = useCompoundsApi();
  const { data: protocols, isLoading } = api.get<Protocol[]>('protocols/');

  const columns: Column<Protocol>[] = [
    {
      key: 'name',
      label: 'Protocol Name',
      sortable: true,
      searchable: true,
      render: (value) => (
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Description fontSize="small" color="primary" />
          <Typography fontWeight={500}>{value}</Typography>
        </Box>
      ),
    },
    {
      key: 'analysis_method',
      label: 'Analysis Method',
      sortable: true,
      render: (value) => (
        <Chip
          label={ANALYSIS_METHOD_LABELS[value] || value}
          size="small"
          variant="outlined"
        />
      ),
    },
    {
      key: 'preferred_dilutions_display',
      label: 'Dilutions',
      render: (value) =>
        value ? (
          <Typography variant="body2" sx={{ fontFamily: 'monospace', fontSize: '0.8rem' }}>
            {value}
          </Typography>
        ) : (
          '-'
        ),
    },
    {
      key: 'assays_count',
      label: 'Assays',
      sortable: true,
      width: 100,
      render: (value) =>
        value !== undefined ? (
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
          { label: 'Assays', href: '/assays/protocols' },
          { label: 'Protocols', icon: 'protocol' },
        ]}
      />

      <Box sx={{ mb: 3 }}>
        <Typography variant="h4" gutterBottom>
          Protocols
        </Typography>
        <Typography color="text.secondary">
          Assay protocols define experimental methods and analysis approaches
        </Typography>
      </Box>

      <DataTable
        data={protocols}
        columns={columns}
        loading={isLoading}
        onRowClick={(protocol) => router.push(routes.assays.protocol(protocol.id))}
        getRowKey={(row) => row.id}
        title={protocols ? `${protocols.length} protocols` : undefined}
        emptyMessage="No protocols found"
      />
    </Container>
  );
}
