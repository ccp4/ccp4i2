'use client';

import { useState } from 'react';
import { useRouter } from 'next/navigation';
import Link from 'next/link';
import { Container, Typography, Box, Chip, Button } from '@mui/material';
import { Science, Add, LocalShipping, Apps } from '@mui/icons-material';
import { useSWRConfig } from 'swr';
import { PageHeader } from '@/components/compounds/PageHeader';
import { DataTable, Column } from '@/components/compounds/DataTable';
import { TargetCreateDialog } from '@/components/compounds/TargetCreateDialog';
import { useCompoundsApi } from '@/lib/compounds/api';
import { routes } from '@/lib/compounds/routes';
import { Target } from '@/types/compounds/models';

export default function TargetsPage() {
  const router = useRouter();
  const { mutate } = useSWRConfig();
  const api = useCompoundsApi();
  const { data: targets, isLoading } = api.get<Target[]>('targets/');
  const [createDialogOpen, setCreateDialogOpen] = useState(false);

  const handleTargetCreated = () => {
    // Invalidate the targets cache to refresh the list
    mutate('/api/proxy/compounds/targets/');
  };

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
      <PageHeader
        breadcrumbs={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Registry', href: routes.registry.targets() },
          { label: 'Targets', icon: 'target' },
        ]}
      />

      <Box sx={{ mb: 3, display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
        <Box>
          <Typography variant="h4" gutterBottom>
            Targets
          </Typography>
          <Typography color="text.secondary">
            Drug discovery targets and campaigns
          </Typography>
        </Box>
        <Box sx={{ display: 'flex', gap: 1 }}>
          <Button
            component={Link}
            href="/"
            variant="outlined"
            startIcon={<Apps />}
          >
            All Apps
          </Button>
          <Button
            component={Link}
            href={routes.registry.suppliers()}
            variant="outlined"
            startIcon={<LocalShipping />}
          >
            Suppliers
          </Button>
          <Button
            variant="contained"
            startIcon={<Add />}
            onClick={() => setCreateDialogOpen(true)}
          >
            New Target
          </Button>
        </Box>
      </Box>

      <DataTable
        data={targets}
        columns={columns}
        loading={isLoading}
        onRowClick={(target) => router.push(routes.registry.target(target.id))}
        getRowKey={(row) => row.id}
        title={targets ? `${targets.length} targets` : undefined}
        emptyMessage="No targets found"
      />

      <TargetCreateDialog
        open={createDialogOpen}
        onClose={() => setCreateDialogOpen(false)}
        onCreated={handleTargetCreated}
      />
    </Container>
  );
}
