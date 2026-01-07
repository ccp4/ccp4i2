'use client';

import { use, useState } from 'react';
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
  Button,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  IconButton,
  Alert,
} from '@mui/material';
import { Description, Science, Assessment, Edit, GridOn, Close, Add } from '@mui/icons-material';
import { Breadcrumbs } from '@/components/Breadcrumbs';
import { DataTable, Column } from '@/components/DataTable';
import { PlatePreview } from '@/components/PlatePreview';
import { PlateLayoutEditor } from '@/components/PlateLayoutEditor';
import { AssayUploadDrawer } from '@/components/AssayUploadDrawer';
import { useCompoundsApi } from '@/lib/api';
import { Protocol, Assay, PlateLayout } from '@/types/models';

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

  const [layoutEditorOpen, setLayoutEditorOpen] = useState(false);
  const [editedLayout, setEditedLayout] = useState<PlateLayout | null>(null);
  const [saveError, setSaveError] = useState<string | null>(null);
  const [saving, setSaving] = useState(false);
  const [uploadDrawerOpen, setUploadDrawerOpen] = useState(false);

  const { data: protocol, isLoading: protocolLoading, mutate } = api.get<Protocol>(
    `protocols/${id}/`
  );
  const { data: assays, isLoading: assaysLoading, mutate: mutateAssays } = api.get<Assay[]>(
    `assays/?protocol=${id}`
  );

  const handleOpenLayoutEditor = () => {
    setEditedLayout(protocol?.plate_layout as PlateLayout || null);
    setLayoutEditorOpen(true);
    setSaveError(null);
  };

  const handleSaveLayout = async () => {
    if (!editedLayout) return;

    setSaving(true);
    setSaveError(null);

    try {
      const response = await fetch(`/api/proxy/compounds/protocols/${id}/`, {
        method: 'PATCH',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ plate_layout: editedLayout }),
      });

      if (!response.ok) {
        throw new Error('Failed to save plate layout');
      }

      await mutate();
      setLayoutEditorOpen(false);
    } catch (err) {
      setSaveError(err instanceof Error ? err.message : 'Failed to save');
    } finally {
      setSaving(false);
    }
  };

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

            {/* Plate Layout Section */}
            <Divider sx={{ my: 2 }} />
            <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', mb: 2 }}>
              <Typography variant="h6">
                Plate Layout
              </Typography>
              <Button
                size="small"
                startIcon={protocol.plate_layout ? <Edit /> : <GridOn />}
                onClick={handleOpenLayoutEditor}
              >
                {protocol.plate_layout ? 'Edit Layout' : 'Configure Layout'}
              </Button>
            </Box>

            {protocol.plate_layout && Object.keys(protocol.plate_layout).length > 0 ? (
              <Box sx={{ display: 'flex', justifyContent: 'center' }}>
                <PlatePreview
                  layout={protocol.plate_layout}
                  width={450}
                  height={300}
                />
              </Box>
            ) : (
              <Paper sx={{ p: 3, bgcolor: 'grey.50', textAlign: 'center' }}>
                <GridOn sx={{ fontSize: 48, color: 'grey.400', mb: 1 }} />
                <Typography color="text.secondary">
                  No plate layout configured
                </Typography>
                <Typography variant="body2" color="text.secondary" sx={{ mt: 0.5 }}>
                  Configure a plate layout to define control positions, sample regions, and dilution patterns
                </Typography>
              </Paper>
            )}
          </>
        ) : (
          <Typography color="error">Protocol not found</Typography>
        )}
      </Paper>

      {/* Plate Layout Editor Dialog */}
      <Dialog
        open={layoutEditorOpen}
        onClose={() => setLayoutEditorOpen(false)}
        maxWidth="lg"
        fullWidth
      >
        <DialogTitle sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
          <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
            <GridOn color="primary" />
            Configure Plate Layout
          </Box>
          <IconButton onClick={() => setLayoutEditorOpen(false)} size="small">
            <Close />
          </IconButton>
        </DialogTitle>
        <DialogContent dividers>
          {saveError && (
            <Alert severity="error" sx={{ mb: 2 }}>
              {saveError}
            </Alert>
          )}
          <PlateLayoutEditor
            value={editedLayout || {}}
            onChange={(layout) => setEditedLayout(layout)}
          />
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setLayoutEditorOpen(false)}>
            Cancel
          </Button>
          <Button
            variant="contained"
            onClick={handleSaveLayout}
            disabled={saving || !editedLayout}
          >
            {saving ? 'Saving...' : 'Save Layout'}
          </Button>
        </DialogActions>
      </Dialog>

      {/* Assays section header with Add button */}
      <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', mb: 2 }}>
        <Typography variant="h6">
          Assays
        </Typography>
        <Button
          variant="contained"
          startIcon={<Add />}
          onClick={() => setUploadDrawerOpen(true)}
          disabled={!protocol}
        >
          Add Assay
        </Button>
      </Box>

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

      {/* Assay Upload Drawer */}
      {protocol && (
        <AssayUploadDrawer
          open={uploadDrawerOpen}
          onClose={() => setUploadDrawerOpen(false)}
          protocolId={id}
          protocolName={protocol.name}
          onAssayCreated={() => mutateAssays()}
        />
      )}
    </Container>
  );
}
