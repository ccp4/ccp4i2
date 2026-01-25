'use client';

import { useState } from 'react';
import { useRouter } from 'next/navigation';
import { Container, Typography, Box, Chip, Button, Tooltip } from '@mui/material';
import { Science, Description, Add, Settings, GridOn } from '@mui/icons-material';
import { useSWRConfig } from 'swr';
import { PageHeader } from '@/components/compounds/PageHeader';
import { DataTable, Column } from '@/components/compounds/DataTable';
import { ProtocolCreateDialog } from '@/components/compounds/ProtocolCreateDialog';
import { useCompoundsApi } from '@/lib/compounds/api';
import { useAuth } from '@/lib/compounds/auth-context';
import { routes } from '@/lib/compounds/routes';
import { Protocol } from '@/types/compounds/models';

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
  const { mutate } = useSWRConfig();
  const { canContribute } = useAuth();
  const [createDialogOpen, setCreateDialogOpen] = useState(false);
  const api = useCompoundsApi();
  const { data: protocols, isLoading } = api.get<Protocol[]>('protocols/');

  const handleProtocolCreated = (protocol: Protocol) => {
    // Refresh the protocols list
    mutate('/api/proxy/compounds/protocols/');
    // Navigate to the new protocol
    router.push(routes.assays.protocol(protocol.id));
  };

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
      <PageHeader
        breadcrumbs={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Assays', href: routes.assays.protocols() },
          { label: 'Protocols', icon: 'protocol' },
        ]}
      />

      <Box sx={{ mb: 3, display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
        <Box>
          <Typography variant="h4" gutterBottom>
            Protocols
          </Typography>
          <Typography color="text.secondary">
            Assay protocols define experimental methods and analysis approaches
          </Typography>
        </Box>
        <Box sx={{ display: 'flex', gap: 1 }}>
          <Button
            variant="outlined"
            startIcon={<Settings />}
            onClick={() => router.push(routes.assays.dilutionSeries())}
          >
            Manage Dilutions
          </Button>
          <Button
            variant="outlined"
            startIcon={<GridOn />}
            onClick={() => router.push(routes.assays.plateLayouts())}
          >
            Manage Plate Layouts
          </Button>
          <Tooltip title={canContribute ? '' : 'Requires Contributor or Admin operating level'} arrow>
            <span>
              <Button
                variant="contained"
                startIcon={<Add />}
                onClick={() => setCreateDialogOpen(true)}
                disabled={!canContribute}
              >
                Add Protocol
              </Button>
            </span>
          </Tooltip>
        </Box>
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

      <ProtocolCreateDialog
        open={createDialogOpen}
        onClose={() => setCreateDialogOpen(false)}
        onCreated={handleProtocolCreated}
      />
    </Container>
  );
}
