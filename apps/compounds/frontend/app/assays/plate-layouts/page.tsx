'use client';

import { useState } from 'react';
import {
  Container,
  Typography,
  Box,
  Chip,
  Button,
  IconButton,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Alert,
  Tooltip,
} from '@mui/material';
import { GridOn, Add, Delete, ArrowBack } from '@mui/icons-material';
import { useSWRConfig } from 'swr';
import { useRouter } from 'next/navigation';
import { PageHeader } from '@/components/compounds/PageHeader';
import { DataTable, Column } from '@/components/compounds/DataTable';
import { PlateLayoutCreateDialog } from '@/components/compounds/PlateLayoutCreateDialog';
import { useCompoundsApi } from '@/lib/compounds/api';
import { useAuth } from '@/lib/compounds/auth-context';
import { routes } from '@/lib/compounds/routes';
import type { PlateLayoutRecord } from '@/types/compounds/models';

export default function PlateLayoutsPage() {
  const router = useRouter();
  const { mutate } = useSWRConfig();
  const { canContribute } = useAuth();
  const [createDialogOpen, setCreateDialogOpen] = useState(false);
  const [deleteDialogOpen, setDeleteDialogOpen] = useState(false);
  const [selectedLayout, setSelectedLayout] = useState<PlateLayoutRecord | null>(null);
  const [deleting, setDeleting] = useState(false);
  const [deleteError, setDeleteError] = useState<string | null>(null);

  const api = useCompoundsApi();
  const { data: plateLayouts, isLoading } = api.get<PlateLayoutRecord[]>('plate-layouts/');

  const handleCreated = () => {
    mutate('/api/proxy/compounds/plate-layouts/');
  };

  const handleRowClick = (layout: PlateLayoutRecord) => {
    router.push(routes.assays.plateLayout(layout.id));
  };

  const handleDeleteClick = (e: React.MouseEvent, layout: PlateLayoutRecord) => {
    e.stopPropagation();
    setSelectedLayout(layout);
    setDeleteError(null);
    setDeleteDialogOpen(true);
  };

  const handleDelete = async () => {
    if (!selectedLayout) return;

    // Check if layout is in use
    if (selectedLayout.protocols_count && selectedLayout.protocols_count > 0) {
      setDeleteError(`Cannot delete: ${selectedLayout.protocols_count} protocol(s) are using this layout.`);
      return;
    }

    setDeleting(true);
    setDeleteError(null);
    try {
      await api.delete(`plate-layouts/${selectedLayout.id}/`);
      mutate('/api/proxy/compounds/plate-layouts/');
      setDeleteDialogOpen(false);
      setSelectedLayout(null);
    } catch (err) {
      console.error('Delete failed:', err);
      setDeleteError('Failed to delete plate layout.');
    } finally {
      setDeleting(false);
    }
  };

  const columns: Column<PlateLayoutRecord>[] = [
    {
      key: 'name',
      label: 'Name',
      sortable: true,
      render: (value, row) => (
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <GridOn fontSize="small" color="primary" />
          <Box>
            <Typography variant="body2" fontWeight="medium">
              {value}
            </Typography>
            {row.description && (
              <Typography variant="caption" color="text.secondary">
                {row.description}
              </Typography>
            )}
          </Box>
        </Box>
      ),
    },
    {
      key: 'plate_format',
      label: 'Format',
      sortable: true,
      width: 100,
      render: (value) => value ? (
        <Chip label={`${value}-well`} size="small" variant="outlined" />
      ) : (
        <Typography color="text.secondary" variant="body2">-</Typography>
      ),
    },
    {
      key: 'protocols_count',
      label: 'Protocols',
      sortable: true,
      width: 100,
      render: (value) => (
        <Chip
          label={value || 0}
          size="small"
          color={value ? 'primary' : 'default'}
          variant="outlined"
        />
      ),
    },
    {
      key: 'is_builtin',
      label: 'Type',
      sortable: true,
      width: 100,
      render: (value) => (
        <Chip
          label={value ? 'Built-in' : 'Custom'}
          size="small"
          color={value ? 'secondary' : 'default'}
          variant="outlined"
        />
      ),
    },
    {
      key: 'created_by_email',
      label: 'Created By',
      sortable: true,
      width: 150,
      render: (value) => (
        <Typography variant="body2" color="text.secondary">
          {value || '-'}
        </Typography>
      ),
    },
    {
      key: 'id',
      label: '',
      width: 50,
      render: (_value, row) => {
        const isDisabled = !canContribute || row.is_builtin || (row.protocols_count ?? 0) > 0;
        const disabledReason = !canContribute
          ? 'Requires Contributor or Admin operating level'
          : row.is_builtin
          ? 'Built-in layouts cannot be deleted'
          : (row.protocols_count ?? 0) > 0
          ? 'Layout is in use by protocols'
          : 'Delete';
        return (
          <Tooltip title={disabledReason}>
            <span>
              <IconButton
                size="small"
                onClick={(e) => handleDeleteClick(e, row)}
                disabled={isDisabled}
                sx={{ opacity: 0.5, '&:hover': { opacity: 1 } }}
              >
                <Delete fontSize="small" />
              </IconButton>
            </span>
          </Tooltip>
        );
      },
    },
  ];

  return (
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <PageHeader
        breadcrumbs={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Assays', href: routes.assays.protocols() },
          { label: 'Protocols', href: routes.assays.protocols() },
          { label: 'Plate Layouts', icon: 'protocol' },
        ]}
      />

      <Box sx={{ mb: 3, display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
        <Box>
          <Typography variant="h4" gutterBottom>
            Plate Layouts
          </Typography>
          <Typography color="text.secondary">
            Reusable plate configurations for assay protocols
          </Typography>
        </Box>
        <Box sx={{ display: 'flex', gap: 1 }}>
          <Button
            variant="outlined"
            startIcon={<ArrowBack />}
            onClick={() => router.push(routes.assays.protocols())}
          >
            Back to Protocols
          </Button>
          <Tooltip title={canContribute ? '' : 'Requires Contributor or Admin operating level'} arrow>
            <span>
              <Button
                variant="contained"
                startIcon={<Add />}
                onClick={() => setCreateDialogOpen(true)}
                disabled={!canContribute}
              >
                Add Plate Layout
              </Button>
            </span>
          </Tooltip>
        </Box>
      </Box>

      <DataTable
        data={plateLayouts}
        columns={columns}
        loading={isLoading}
        onRowClick={handleRowClick}
        getRowKey={(row) => row.id}
        title={plateLayouts ? `${plateLayouts.length} plate layouts` : undefined}
        emptyMessage="No plate layouts found. Click 'Add Plate Layout' to create one."
      />

      <PlateLayoutCreateDialog
        open={createDialogOpen}
        onClose={() => setCreateDialogOpen(false)}
        onCreated={handleCreated}
      />

      <Dialog open={deleteDialogOpen} onClose={() => setDeleteDialogOpen(false)}>
        <DialogTitle>Delete Plate Layout?</DialogTitle>
        <DialogContent>
          {deleteError && (
            <Alert severity="error" sx={{ mb: 2 }}>
              {deleteError}
            </Alert>
          )}
          <Typography>
            Are you sure you want to delete this plate layout?
          </Typography>
          {selectedLayout && (
            <Typography sx={{ fontWeight: 'medium', mt: 1 }}>
              {selectedLayout.name}
            </Typography>
          )}
          {selectedLayout?.protocols_count && selectedLayout.protocols_count > 0 && (
            <Alert severity="warning" sx={{ mt: 2 }}>
              This layout is used by {selectedLayout.protocols_count} protocol(s) and cannot be deleted.
            </Alert>
          )}
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setDeleteDialogOpen(false)} disabled={deleting}>
            Cancel
          </Button>
          <Button
            onClick={handleDelete}
            color="error"
            variant="contained"
            disabled={deleting || (selectedLayout?.protocols_count ?? 0) > 0}
          >
            {deleting ? 'Deleting...' : 'Delete'}
          </Button>
        </DialogActions>
      </Dialog>
    </Container>
  );
}
