'use client';

import { useState } from 'react';
import { useRouter } from 'next/navigation';
import { Container, Typography, Box, Chip, Button } from '@mui/material';
import { LocalShipping, Add, Person } from '@mui/icons-material';
import { useSWRConfig } from 'swr';
import { PageHeader } from '@/components/compounds/PageHeader';
import { DataTable, Column } from '@/components/compounds/DataTable';
import { SupplierCreateDialog } from '@/components/compounds/SupplierCreateDialog';
import { useCompoundsApi } from '@/lib/compounds/api';
import { routes } from '@/lib/compounds/routes';
import { Supplier } from '@/types/compounds/models';

export default function SuppliersPage() {
  const router = useRouter();
  const { mutate } = useSWRConfig();
  const api = useCompoundsApi();
  const { data: suppliers, isLoading } = api.get<Supplier[]>('suppliers/');
  const [createDialogOpen, setCreateDialogOpen] = useState(false);

  const handleSupplierCreated = () => {
    // Invalidate the suppliers cache to refresh the list
    mutate('/api/proxy/compounds/suppliers/');
  };

  const columns: Column<Supplier>[] = [
    {
      key: 'name',
      label: 'Supplier Name',
      sortable: true,
      searchable: true,
      render: (value, row) => (
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <LocalShipping fontSize="small" color="primary" />
          <Typography fontWeight={500}>{value}</Typography>
          {row.is_current_user && (
            <Chip
              icon={<Person sx={{ fontSize: 14 }} />}
              label="You"
              size="small"
              color="primary"
              variant="outlined"
              sx={{ ml: 1 }}
            />
          )}
        </Box>
      ),
    },
    {
      key: 'initials',
      label: 'Initials',
      sortable: true,
      width: 100,
      render: (value) => value || '-',
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
          <Chip label="0" size="small" variant="outlined" />
        ),
    },
    {
      key: 'batch_count',
      label: 'Batches',
      sortable: true,
      width: 100,
      render: (value) =>
        value ? (
          <Chip label={value} size="small" color="secondary" variant="outlined" />
        ) : (
          <Chip label="0" size="small" variant="outlined" />
        ),
    },
  ];

  return (
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <PageHeader
        breadcrumbs={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Registry', href: routes.registry.targets() },
          { label: 'Suppliers', icon: 'supplier' },
        ]}
      />

      <Box sx={{ mb: 3, display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
        <Box>
          <Typography variant="h4" gutterBottom>
            Suppliers
          </Typography>
          <Typography color="text.secondary">
            Chemical suppliers and synthesis sources
          </Typography>
        </Box>
        <Button
          variant="contained"
          startIcon={<Add />}
          onClick={() => setCreateDialogOpen(true)}
        >
          New Supplier
        </Button>
      </Box>

      <DataTable
        data={suppliers}
        columns={columns}
        loading={isLoading}
        onRowClick={(supplier) => router.push(routes.registry.supplier(supplier.id))}
        getRowKey={(row) => row.id}
        title={suppliers ? `${suppliers.length} suppliers` : undefined}
        emptyMessage="No suppliers found"
      />

      <SupplierCreateDialog
        open={createDialogOpen}
        onClose={() => setCreateDialogOpen(false)}
        onCreated={handleSupplierCreated}
      />
    </Container>
  );
}
