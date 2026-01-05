'use client';

import { use } from 'react';
import { useRouter } from 'next/navigation';
import {
  Container,
  Typography,
  Box,
  Paper,
  Grid,
  Chip,
  Skeleton,
  Divider,
} from '@mui/material';
import { Description, Science, Assessment } from '@mui/icons-material';
import { Breadcrumbs } from '@/components/Breadcrumbs';
import { DataTable, Column } from '@/components/DataTable';
import { useCompoundsApi } from '@/lib/api';
import { Protocol, Assay } from '@/types/models';

interface PageProps {
  params: Promise<{ id: string }>;
}

const ANALYSIS_METHOD_LABELS: Record<string, string> = {
  hill_langmuir: 'Hill-Langmuir',
  hill_langmuir_fix_hill: 'Hill-Langmuir (fixed Hill)',
  hill_langmuir_fix_hill_minmax: 'Hill-Langmuir (fixed Hill/min/max)',
  hill_langmuir_fix_minmax: 'Hill-Langmuir (fixed min/max)',
  ms_intact: 'MS-Intact',
  table_of_values: 'Table of values',
};

function InfoRow({ label, value }: { label: string; value: React.ReactNode }) {
  return (
    <Box sx={{ display: 'flex', py: 0.5 }}>
      <Typography
        color="text.secondary"
        sx={{ minWidth: 140, fontWeight: 500 }}
      >
        {label}:
      </Typography>
      <Typography component="div">{value ?? '-'}</Typography>
    </Box>
  );
}

export default function ProtocolDetailPage({ params }: PageProps) {
  const { id } = use(params);
  const router = useRouter();
  const api = useCompoundsApi();

  const { data: protocol, isLoading: protocolLoading } = api.get<Protocol>(
    `protocols/${id}/`
  );
  const { data: assays, isLoading: assaysLoading } = api.get<Assay[]>(
    `assays/?protocol=${id}`
  );

  const columns: Column<Assay>[] = [
    {
      key: 'data_filename',
      label: 'Data File',
      sortable: true,
      searchable: true,
      render: (value) => (
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Assessment fontSize="small" color="info" />
          <Typography fontWeight={500}>{value || 'No file'}</Typography>
        </Box>
      ),
    },
    {
      key: 'target_name',
      label: 'Target',
      sortable: true,
      searchable: true,
      render: (value, row) =>
        value ? (
          <Chip
            icon={<Science fontSize="small" />}
            label={value}
            size="small"
            variant="outlined"
            onClick={(e) => {
              e.stopPropagation();
              if (row.target) router.push(`/registry/targets/${row.target}`);
            }}
          />
        ) : (
          '-'
        ),
    },
    {
      key: 'data_series_count',
      label: 'Series',
      sortable: true,
      width: 80,
      render: (value) =>
        value !== undefined ? (
          <Chip label={value} size="small" color="secondary" variant="outlined" />
        ) : (
          '-'
        ),
    },
    {
      key: 'created_by_email',
      label: 'Created By',
      searchable: true,
      width: 180,
      render: (value) => value || '-',
    },
    {
      key: 'created_at',
      label: 'Date',
      sortable: true,
      width: 100,
      render: (value) =>
        value ? new Date(value).toLocaleDateString() : '-',
    },
  ];

  return (
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <Breadcrumbs
        items={[
          { label: 'Home', href: '/', icon: 'home' },
          { label: 'Protocols', href: '/assays/protocols', icon: 'protocol' },
          { label: protocol?.name || 'Loading...', icon: 'protocol' },
        ]}
      />

      {/* Protocol header */}
      <Paper sx={{ p: 3, mb: 3 }}>
        {protocolLoading ? (
          <>
            <Skeleton variant="text" width={300} height={40} />
            <Skeleton variant="text" width={200} />
          </>
        ) : protocol ? (
          <>
            <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 2 }}>
              <Description sx={{ fontSize: 48, color: 'primary.main' }} />
              <Box>
                <Typography variant="h4">{protocol.name}</Typography>
                <Box sx={{ display: 'flex', gap: 1, mt: 0.5 }}>
                  <Chip
                    label={ANALYSIS_METHOD_LABELS[protocol.analysis_method] || protocol.analysis_method}
                    size="small"
                    color="primary"
                    variant="outlined"
                  />
                </Box>
              </Box>
            </Box>

            <Divider sx={{ my: 2 }} />

            <Grid container spacing={3}>
              <Grid item xs={12} md={6}>
                <Typography variant="h6" gutterBottom>
                  Configuration
                </Typography>
                <InfoRow
                  label="Analysis"
                  value={ANALYSIS_METHOD_LABELS[protocol.analysis_method]}
                />
                <InfoRow
                  label="Dilutions"
                  value={
                    protocol.preferred_dilutions_display ? (
                      <Typography
                        component="span"
                        sx={{ fontFamily: 'monospace', fontSize: '0.85rem' }}
                      >
                        {protocol.preferred_dilutions_display}
                      </Typography>
                    ) : null
                  }
                />
                <InfoRow
                  label="PHERAstar Table"
                  value={protocol.pherastar_table}
                />
              </Grid>
              <Grid item xs={12} md={6}>
                <Typography variant="h6" gutterBottom>
                  Metadata
                </Typography>
                <InfoRow label="Created By" value={protocol.created_by_email} />
                <InfoRow
                  label="Created"
                  value={
                    protocol.created_at
                      ? new Date(protocol.created_at).toLocaleString()
                      : null
                  }
                />
                <InfoRow
                  label="Assays"
                  value={assays ? `${assays.length} experiments` : 'Loading...'}
                />
              </Grid>
            </Grid>

            {protocol.comments && (
              <>
                <Divider sx={{ my: 2 }} />
                <Typography variant="h6" gutterBottom>
                  Comments
                </Typography>
                <Typography>{protocol.comments}</Typography>
              </>
            )}
          </>
        ) : (
          <Typography color="error">Protocol not found</Typography>
        )}
      </Paper>

      {/* Assays table */}
      <DataTable
        data={assays}
        columns={columns}
        loading={assaysLoading}
        onRowClick={(assay) => router.push(`/assays/${assay.id}`)}
        getRowKey={(row) => row.id}
        title={assays ? `${assays.length} assays` : undefined}
        emptyMessage="No assays using this protocol"
      />
    </Container>
  );
}
