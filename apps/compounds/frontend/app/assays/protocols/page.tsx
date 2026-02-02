'use client';

import { useState } from 'react';
import { useRouter } from 'next/navigation';
import { Container, Typography, Box, Chip, Button, Tooltip } from '@mui/material';
import { Science, Description, Add, Settings, GridOn, FiberNew } from '@mui/icons-material';
import { useSWRConfig } from 'swr';
import { PageHeader } from '@/components/compounds/PageHeader';
import { DataTable, Column } from '@/components/data-table';
import { ProtocolCreateDialog } from '@/components/compounds/ProtocolCreateDialog';
import { useCompoundsApi } from '@/lib/compounds/api';
import { useAuth } from '@/lib/compounds/auth-context';
import { routes } from '@/lib/compounds/routes';
import { Protocol, ImportType } from '@/types/compounds/models';

const IMPORT_TYPE_LABELS: Record<ImportType, string> = {
  raw_data: 'Raw Data (Dose-Response)',
  ms_intact: 'MS-Intact',
  table_of_values: 'Table of Values',
  pharmaron_adme: 'Pharmaron ADME',
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
      render: (value, row) => (
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Description fontSize="small" color="primary" />
          <Typography fontWeight={500}>{value}</Typography>
          {row.has_recent_assays && (
            <Tooltip title="New assays in the last 7 days">
              <FiberNew fontSize="small" color="secondary" />
            </Tooltip>
          )}
        </Box>
      ),
    },
    {
      key: 'import_type',
      label: 'Data Type',
      sortable: true,
      hiddenOnMobile: true,
      render: (value: ImportType) => (
        <Chip
          label={IMPORT_TYPE_LABELS[value] || value}
          size="small"
          variant="outlined"
        />
      ),
    },
    {
      key: 'preferred_dilutions_display',
      label: 'Dilutions',
      hiddenOnMobile: true,
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
      hiddenOnMobile: true,
      render: (value) =>
        value ? new Date(value).toLocaleDateString() : '-',
    },
  ];

  return (
    <Box sx={{ display: 'flex', flexDirection: 'column', height: '100vh', overflow: 'hidden' }}>
      <Container maxWidth="lg" sx={{ flex: 1, display: 'flex', flexDirection: 'column', overflow: 'hidden', py: 2 }}>
        <Box sx={{ flexShrink: 0 }}>
          <PageHeader
            breadcrumbs={[
              { label: 'Home', href: routes.home(), icon: 'home' },
              { label: 'Assays', href: routes.assays.protocols() },
              { label: 'Protocols', icon: 'protocol' },
            ]}
          />

          <Box sx={{ mb: 2, display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
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
        </Box>

        <Box sx={{ flex: 1, overflow: 'hidden', minHeight: 0 }}>
          <DataTable
            data={protocols}
            columns={columns}
            loading={isLoading}
            onRowClick={(protocol) => router.push(routes.assays.protocol(protocol.id))}
            getRowKey={(row) => row.id}
            title={protocols ? `${protocols.length} protocols` : undefined}
            emptyMessage="No protocols found"
            fillHeight
          />
        </Box>
      </Container>

      <ProtocolCreateDialog
        open={createDialogOpen}
        onClose={() => setCreateDialogOpen(false)}
        onCreated={handleProtocolCreated}
      />
    </Box>
  );
}
