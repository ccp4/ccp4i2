'use client';

import { useRouter } from 'next/navigation';
import {
  Container,
  Typography,
  Box,
  Chip,
  Button,
  IconButton,
  Tooltip,
} from '@mui/material';
import { Science, Add, Delete } from '@mui/icons-material';
import { useState } from 'react';
import { useSWRConfig } from 'swr';
import { PageHeader } from '@/components/compounds/PageHeader';
import { DataTable, Column } from '@/components/data-table';
import { PlasmidCreateDialog } from '@/components/compounds/PlasmidCreateDialog';
import { ConfirmDialog } from '@/components/compounds/ConfirmDialog';
import { useCompoundsApi } from '@/lib/compounds/api';
import { useAuth } from '@/lib/compounds/auth-context';
import { routes } from '@/lib/compounds/routes';
import { Plasmid } from '@/types/compounds/constructs';

export default function ConstructsPage() {
  const router = useRouter();
  const api = useCompoundsApi();
  const { mutate } = useSWRConfig();
  const { canContribute } = useAuth();
  const { data: plasmids, isLoading } = api.get<Plasmid[]>('plasmids/');
  const [createDialogOpen, setCreateDialogOpen] = useState(false);
  const [deleteTarget, setDeleteTarget] = useState<Plasmid | null>(null);
  const [deleting, setDeleting] = useState(false);

  const handleDelete = async () => {
    if (!deleteTarget) return;
    setDeleting(true);
    try {
      await api.delete(`plasmids/${deleteTarget.id}/`);
      mutate('/api/proxy/compounds/plasmids/');
      setDeleteTarget(null);
    } catch (error) {
      console.error('Failed to delete plasmid:', error);
    } finally {
      setDeleting(false);
    }
  };

  const columns: Column<Plasmid>[] = [
    {
      key: 'formatted_id',
      label: 'Construct ID',
      sortable: true,
      searchable: true,
      render: (value) => (
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Science fontSize="small" color="primary" />
          <Typography fontFamily="monospace" fontWeight={500}>
            {value}
          </Typography>
        </Box>
      ),
    },
    {
      key: 'name',
      label: 'Name',
      sortable: true,
      searchable: true,
    },
    {
      key: 'project_name',
      label: 'Project',
      sortable: true,
      searchable: true,
      render: (value) => value || '-',
    },
    {
      key: 'cassette_count',
      label: 'Cassettes',
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
    {
      key: 'id',
      label: '',
      width: 50,
      render: (_, row) => (
        <Tooltip title={canContribute ? 'Delete' : 'Requires Contributor or Admin operating level'}>
          <span>
            <IconButton
              size="small"
              color="error"
              disabled={!canContribute}
              onClick={(e) => {
                e.stopPropagation();
                setDeleteTarget(row);
              }}
            >
              <Delete fontSize="small" />
            </IconButton>
          </span>
        </Tooltip>
      ),
    },
  ];

  return (
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <PageHeader
        breadcrumbs={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Constructs', icon: 'construct' },
        ]}
      />

      <Box sx={{ mb: 3, display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
        <Box>
          <Typography variant="h4" gutterBottom>
            Construct Database
          </Typography>
          <Typography color="text.secondary">
            Plasmid constructs, cassettes, and sequencing results
          </Typography>
        </Box>
        <Tooltip title={canContribute ? '' : 'Requires Contributor or Admin operating level'} arrow>
          <span>
            <Button
              variant="contained"
              startIcon={<Add />}
              onClick={() => setCreateDialogOpen(true)}
              disabled={!canContribute}
            >
              Add Construct
            </Button>
          </span>
        </Tooltip>
      </Box>

      <DataTable
        data={plasmids}
        columns={columns}
        loading={isLoading}
        onRowClick={(plasmid) => router.push(routes.constructs.plasmid(plasmid.id))}
        getRowKey={(row) => row.id}
        title={plasmids ? `${plasmids.length} plasmids` : undefined}
        emptyMessage="No plasmids found. Click 'Add Construct' to create one."
      />

      <PlasmidCreateDialog
        open={createDialogOpen}
        onClose={() => setCreateDialogOpen(false)}
        onCreated={(newPlasmid) => {
          // Refresh the plasmids list
          mutate('/api/proxy/compounds/plasmids/');
          // Navigate to the new plasmid
          router.push(routes.constructs.plasmid(newPlasmid.id));
        }}
      />

      <ConfirmDialog
        open={!!deleteTarget}
        title="Delete Plasmid"
        message={
          <>
            Are you sure you want to delete{' '}
            <strong>{deleteTarget?.formatted_id}</strong> ({deleteTarget?.name})?
            <br />
            <br />
            This action cannot be undone.
          </>
        }
        confirmLabel="Delete"
        loading={deleting}
        onConfirm={handleDelete}
        onCancel={() => setDeleteTarget(null)}
      />
    </Container>
  );
}
