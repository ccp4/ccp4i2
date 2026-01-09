'use client';

import { useRouter } from 'next/navigation';
import Link from 'next/link';
import { Container, Typography, Box, Chip, Button } from '@mui/material';
import { Assessment, Science, Description, Upload, Functions } from '@mui/icons-material';
import { Breadcrumbs } from '@/components/Breadcrumbs';
import { DataTable, Column } from '@/components/DataTable';
import { useCompoundsApi } from '@/lib/api';
import { routes } from '@/lib/routes';
import { Assay } from '@/types/models';

export default function AssaysPage() {
  const router = useRouter();
  const api = useCompoundsApi();
  const { data: assays, isLoading } = api.get<Assay[]>('assays/');

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
      key: 'protocol_name',
      label: 'Protocol',
      sortable: true,
      searchable: true,
      render: (value, row) => (
        <Chip
          icon={<Description fontSize="small" />}
          label={value}
          size="small"
          variant="outlined"
          onClick={(e) => {
            e.stopPropagation();
            router.push(routes.assays.protocol(row.protocol));
          }}
        />
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
              if (row.target) router.push(routes.registry.target(row.target));
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
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Assays', icon: 'assay' },
        ]}
      />

      <Box sx={{ mb: 3, display: 'flex', alignItems: 'flex-start', gap: 2 }}>
        <Box sx={{ flex: 1 }}>
          <Typography variant="h4" gutterBottom>
            Assays
          </Typography>
          <Typography color="text.secondary">
            Dose-response experiments and compound activity data
          </Typography>
        </Box>
        <Box sx={{ display: 'flex', gap: 1 }}>
          <Button
            component={Link}
            href={routes.assays.aggregate()}
            variant="outlined"
            startIcon={<Functions />}
          >
            Aggregate
          </Button>
          <Button
            component={Link}
            href={routes.assays.import()}
            variant="contained"
            startIcon={<Upload />}
          >
            Import Data
          </Button>
        </Box>
      </Box>

      <DataTable
        data={assays}
        columns={columns}
        loading={isLoading}
        onRowClick={(assay) => router.push(routes.assays.detail(assay.id))}
        getRowKey={(row) => row.id}
        title={assays ? `${assays.length} assays` : undefined}
        emptyMessage="No assays found"
      />
    </Container>
  );
}
