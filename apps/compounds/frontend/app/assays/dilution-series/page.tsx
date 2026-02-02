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
  Tooltip,
} from '@mui/material';
import { Science, Add, Delete, ArrowBack } from '@mui/icons-material';
import { useSWRConfig } from 'swr';
import { useRouter } from 'next/navigation';
import { PageHeader } from '@/components/compounds/PageHeader';
import { DataTable, Column } from '@/components/data-table';
import { DilutionSeriesCreateDialog } from '@/components/compounds/DilutionSeriesCreateDialog';
import { DilutionSeriesEditDialog } from '@/components/compounds/DilutionSeriesEditDialog';
import { useCompoundsApi } from '@/lib/compounds/api';
import { useAuth } from '@/lib/compounds/auth-context';
import { routes } from '@/lib/compounds/routes';
import type { DilutionSeries } from '@/types/compounds/models';

export default function DilutionSeriesPage() {
  const router = useRouter();
  const { mutate } = useSWRConfig();
  const { canContribute } = useAuth();
  const [createDialogOpen, setCreateDialogOpen] = useState(false);
  const [editDialogOpen, setEditDialogOpen] = useState(false);
  const [deleteDialogOpen, setDeleteDialogOpen] = useState(false);
  const [selectedSeries, setSelectedSeries] = useState<DilutionSeries | null>(null);
  const [deleting, setDeleting] = useState(false);

  const api = useCompoundsApi();
  const { data: dilutionSeries, isLoading } = api.get<DilutionSeries[]>('dilution-series/');

  const handleCreated = () => {
    mutate('/api/proxy/compounds/dilution-series/');
  };

  const handleSaved = () => {
    mutate('/api/proxy/compounds/dilution-series/');
  };

  const handleRowClick = (series: DilutionSeries) => {
    setSelectedSeries(series);
    setEditDialogOpen(true);
  };

  const handleDeleteClick = (e: React.MouseEvent, series: DilutionSeries) => {
    e.stopPropagation();
    setSelectedSeries(series);
    setDeleteDialogOpen(true);
  };

  const handleDelete = async () => {
    if (!selectedSeries) return;

    setDeleting(true);
    try {
      await api.delete(`dilution-series/${selectedSeries.id}/`);
      mutate('/api/proxy/compounds/dilution-series/');
      setDeleteDialogOpen(false);
      setSelectedSeries(null);
    } catch (err) {
      console.error('Delete failed:', err);
    } finally {
      setDeleting(false);
    }
  };

  const columns: Column<DilutionSeries>[] = [
    {
      key: 'concentrations',
      label: 'Concentrations',
      sortable: false,
      render: (value: number[]) => (
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Science fontSize="small" color="primary" />
          <Typography sx={{ fontFamily: 'monospace', fontSize: '0.85rem' }}>
            {value.join(', ')}
          </Typography>
        </Box>
      ),
    },
    {
      key: 'unit',
      label: 'Unit',
      sortable: true,
      width: 100,
      render: (value) => (
        <Chip label={value} size="small" variant="outlined" />
      ),
    },
    {
      key: 'concentrations',
      label: 'Points',
      sortable: false,
      width: 80,
      render: (value: number[]) => (
        <Chip label={value.length} size="small" color="primary" variant="outlined" />
      ),
    },
    {
      key: 'id',
      label: '',
      width: 50,
      render: (_value, row) => (
        <Tooltip title={canContribute ? 'Delete' : 'Requires Contributor or Admin operating level'}>
          <span>
            <IconButton
              size="small"
              onClick={(e) => handleDeleteClick(e, row)}
              disabled={!canContribute}
              sx={{ opacity: 0.5, '&:hover': { opacity: 1 } }}
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
          { label: 'Assays', href: routes.assays.protocols() },
          { label: 'Protocols', href: routes.assays.protocols() },
          { label: 'Dilution Series', icon: 'protocol' },
        ]}
      />

      <Box sx={{ mb: 3, display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
        <Box>
          <Typography variant="h4" gutterBottom>
            Dilution Series
          </Typography>
          <Typography color="text.secondary">
            Standard concentration series for dose-response experiments
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
                Add Dilution Series
              </Button>
            </span>
          </Tooltip>
        </Box>
      </Box>

      <DataTable
        data={dilutionSeries}
        columns={columns}
        loading={isLoading}
        onRowClick={handleRowClick}
        getRowKey={(row) => row.id}
        title={dilutionSeries ? `${dilutionSeries.length} dilution series` : undefined}
        emptyMessage="No dilution series found"
      />

      <DilutionSeriesCreateDialog
        open={createDialogOpen}
        onClose={() => setCreateDialogOpen(false)}
        onCreated={handleCreated}
      />

      {selectedSeries && (
        <DilutionSeriesEditDialog
          open={editDialogOpen}
          onClose={() => {
            setEditDialogOpen(false);
            setSelectedSeries(null);
          }}
          dilutionSeries={selectedSeries}
          onSave={handleSaved}
        />
      )}

      <Dialog open={deleteDialogOpen} onClose={() => setDeleteDialogOpen(false)}>
        <DialogTitle>Delete Dilution Series?</DialogTitle>
        <DialogContent>
          <Typography>
            Are you sure you want to delete this dilution series?
          </Typography>
          {selectedSeries && (
            <Typography sx={{ fontFamily: 'monospace', mt: 1 }}>
              {selectedSeries.concentrations.join(', ')} {selectedSeries.unit}
            </Typography>
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
            disabled={deleting}
          >
            {deleting ? 'Deleting...' : 'Delete'}
          </Button>
        </DialogActions>
      </Dialog>
    </Container>
  );
}
