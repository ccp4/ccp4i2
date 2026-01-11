'use client';

import { use, useState } from 'react';
import { useRouter } from 'next/navigation';
import {
  Container,
  Typography,
  Box,
  Paper,
  Grid2 as Grid,
  Chip,
  Skeleton,
  Divider,
  Link as MuiLink,
  Button,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogContentText,
  DialogActions,
  CircularProgress,
} from '@mui/material';
import {
  Assessment,
  Science,
  Description,
  Medication,
  CheckCircle,
  Cancel,
  HelpOutline,
  Delete,
} from '@mui/icons-material';
import { Breadcrumbs } from '@/components/compounds/Breadcrumbs';
import { DataTable, Column } from '@/components/compounds/DataTable';
import { DoseResponseThumb } from '@/components/compounds/DoseResponseChart';
import { CompoundStructureCell } from '@/components/compounds/CompoundStructureCell';
import { useCompoundsApi } from '@/lib/compounds/api';
import { routes } from '@/lib/compounds/routes';
import { Assay, DataSeries, Protocol, Target } from '@/types/compounds/models';

interface PageProps {
  params: Promise<{ id: string }>;
}

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

function StatusChip({ status }: { status: string | undefined }) {
  if (!status) {
    return <Chip icon={<HelpOutline />} label="Pending" size="small" />;
  }

  switch (status) {
    case 'valid':
      return (
        <Chip
          icon={<CheckCircle />}
          label="Valid"
          size="small"
          color="success"
        />
      );
    case 'invalid':
      return (
        <Chip icon={<Cancel />} label="Invalid" size="small" color="error" />
      );
    default:
      return <Chip icon={<HelpOutline />} label="Unassigned" size="small" />;
  }
}

export default function AssayDetailPage({ params }: PageProps) {
  const { id } = use(params);
  const router = useRouter();
  const api = useCompoundsApi();
  const [deleteDialogOpen, setDeleteDialogOpen] = useState(false);
  const [deleting, setDeleting] = useState(false);

  const { data: assay, isLoading: assayLoading } = api.get<Assay>(
    `assays/${id}/`
  );
  const { data: dataSeries, isLoading: dataSeriesLoading } = api.get<DataSeries[]>(
    `data-series/?assay=${id}`
  );
  const { data: protocol } = api.get<Protocol>(
    assay?.protocol ? `protocols/${assay.protocol}/` : null
  );
  const { data: target } = api.get<Target>(
    assay?.target ? `targets/${assay.target}/` : null
  );

  const handleDelete = async () => {
    setDeleting(true);
    try {
      const response = await fetch(`/api/proxy/compounds/assays/${id}/`, {
        method: 'DELETE',
      });
      if (response.ok || response.status === 204) {
        // Navigate back to protocol page or assays list
        if (protocol) {
          router.push(routes.assays.protocol(protocol.id));
        } else {
          router.push(routes.assays.list());
        }
      } else {
        console.error('Delete failed:', await response.text());
      }
    } catch (err) {
      console.error('Delete error:', err);
    } finally {
      setDeleting(false);
      setDeleteDialogOpen(false);
    }
  };

  // Helper to get chart data from a data series - requires dilution_series
  const getChartData = (row: DataSeries) => {
    if (!row.dilution_series?.concentrations || !row.extracted_data) return null;

    const concentrations = row.dilution_series.concentrations;
    let responses = Array.isArray(row.extracted_data) ? row.extracted_data : [];
    const unit = row.dilution_series.unit || 'nM';

    // Detect format: if extracted_data has 2 more elements than concentrations,
    // it has embedded controls at first and last positions
    if (responses.length === concentrations.length + 2) {
      responses = responses.slice(1, -1);
    }

    return { concentrations, responses, unit };
  };

  const columns: Column<DataSeries>[] = [
    {
      key: 'chart',
      label: 'Curve',
      width: 140,
      render: (_, row) => {
        const chartData = getChartData(row);
        if (!chartData || chartData.concentrations.length === 0) {
          return <Box sx={{ width: 120, height: 120, bgcolor: 'grey.100', borderRadius: 1 }} />;
        }
        const fitParams = row.analysis?.results ? {
          ec50: row.analysis.results.EC50,
          hill: row.analysis.results.Hill,
          minVal: row.analysis.results.minVal,
          maxVal: row.analysis.results.maxVal,
          status: row.analysis.status,
        } : undefined;
        return <DoseResponseThumb data={chartData} fit={fitParams} size={120} />;
      },
    },
    {
      key: 'structure',
      label: 'Structure',
      width: 110,
      render: (_, row) => (
        <CompoundStructureCell
          compoundId={row.compound}
          compoundName={row.compound_name}
          size={100}
        />
      ),
    },
    {
      key: 'compound_name',
      label: 'Compound Name',
      sortable: true,
      searchable: true,
      render: (value) => (
        <Typography fontWeight={500}>{value || 'Unknown'}</Typography>
      ),
    },
    {
      key: 'compound_formatted_id',
      label: 'Compound ID',
      sortable: true,
      searchable: true,
      render: (value, row) =>
        value ? (
          <Chip
            icon={<Medication fontSize="small" />}
            label={value}
            size="small"
            variant="outlined"
            onClick={(e) => {
              e.stopPropagation();
              if (row.compound) router.push(routes.registry.compound(row.compound));
            }}
          />
        ) : (
          <Typography color="text.secondary" variant="body2">
            Not matched
          </Typography>
        ),
    },
    {
      key: 'analysis_status',
      label: 'Status',
      sortable: true,
      width: 120,
      render: (value) => <StatusChip status={value} />,
    },
    {
      key: 'analysis_kpi',
      label: 'KPI',
      sortable: true,
      width: 100,
      render: (value) =>
        value !== null && value !== undefined ? (
          <Typography fontFamily="monospace" fontWeight={500}>
            {typeof value === 'number' ? value.toFixed(2) : value}
          </Typography>
        ) : (
          '-'
        ),
    },
    {
      key: 'row',
      label: 'Position',
      width: 120,
      render: (value, row) => (
        <Typography variant="body2" color="text.secondary" fontFamily="monospace">
          Row {value}, Col {row.start_column}-{row.end_column}
        </Typography>
      ),
    },
  ];

  return (
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <Breadcrumbs
        items={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Protocols', href: routes.assays.protocols(), icon: 'protocol' },
          ...(protocol ? [{ label: protocol.name, href: routes.assays.protocol(protocol.id), icon: 'protocol' as const }] : []),
          { label: assay?.data_filename || 'Loading...', icon: 'assay' },
        ]}
      />

      {/* Assay header */}
      <Paper sx={{ p: 3, mb: 3 }}>
        {assayLoading ? (
          <>
            <Skeleton variant="text" width={300} height={40} />
            <Skeleton variant="text" width={200} />
          </>
        ) : assay ? (
          <>
            <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 2 }}>
              <Assessment sx={{ fontSize: 48, color: 'info.main' }} />
              <Box sx={{ flex: 1 }}>
                <Typography variant="h4">{assay.data_filename || 'Assay'}</Typography>
                <Box sx={{ display: 'flex', gap: 1, mt: 0.5, flexWrap: 'wrap' }}>
                  {protocol && (
                    <Chip
                      icon={<Description fontSize="small" />}
                      label={protocol.name}
                      size="small"
                      onClick={() => router.push(routes.assays.protocol(protocol.id))}
                    />
                  )}
                  {target && (
                    <Chip
                      icon={<Science fontSize="small" />}
                      label={target.name}
                      size="small"
                      variant="outlined"
                      onClick={() => router.push(routes.registry.target(target.id))}
                    />
                  )}
                </Box>
              </Box>
              <Button
                variant="outlined"
                color="error"
                startIcon={<Delete />}
                onClick={() => setDeleteDialogOpen(true)}
                size="small"
              >
                Delete
              </Button>
            </Box>

            <Divider sx={{ my: 2 }} />

            <Grid container spacing={3}>
              <Grid size={{ xs: 12, md: 6 }}>
                <Typography variant="h6" gutterBottom>
                  Experiment Details
                </Typography>
                <InfoRow label="Protocol" value={protocol?.name} />
                <InfoRow label="Target" value={target?.name} />
                <InfoRow label="Lab Book" value={assay.labbook_number} />
                <InfoRow label="Page" value={assay.page_number} />
              </Grid>
              <Grid size={{ xs: 12, md: 6 }}>
                <Typography variant="h6" gutterBottom>
                  Metadata
                </Typography>
                <InfoRow label="Created By" value={assay.created_by_email} />
                <InfoRow
                  label="Created"
                  value={
                    assay.created_at
                      ? new Date(assay.created_at).toLocaleString()
                      : null
                  }
                />
                <InfoRow
                  label="Data Series"
                  value={dataSeries ? `${dataSeries.length} compounds` : 'Loading...'}
                />
                {assay.data_file && (
                  <InfoRow
                    label="Data File"
                    value={
                      <MuiLink href={assay.data_file} target="_blank">
                        Download
                      </MuiLink>
                    }
                  />
                )}
              </Grid>
            </Grid>

            {assay.comments && (
              <>
                <Divider sx={{ my: 2 }} />
                <Typography variant="h6" gutterBottom>
                  Comments
                </Typography>
                <Typography>{assay.comments}</Typography>
              </>
            )}
          </>
        ) : (
          <Typography color="error">Assay not found</Typography>
        )}
      </Paper>

      {/* Data series table */}
      <DataTable
        data={dataSeries}
        columns={columns}
        loading={dataSeriesLoading}
        getRowKey={(row) => row.id}
        title={dataSeries ? `${dataSeries.length} data series` : undefined}
        emptyMessage="No data series in this assay"
        onRowClick={(row) => router.push(routes.assays.dataSeries(row.id))}
      />

      {/* Delete confirmation dialog */}
      <Dialog open={deleteDialogOpen} onClose={() => setDeleteDialogOpen(false)}>
        <DialogTitle>Delete Assay?</DialogTitle>
        <DialogContent>
          <DialogContentText>
            This will permanently delete the assay &quot;{assay?.data_filename}&quot;
            and all {dataSeries?.length || 0} data series with their analysis results.
            Dilution series will not be affected.
          </DialogContentText>
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
            startIcon={deleting ? <CircularProgress size={16} /> : <Delete />}
          >
            {deleting ? 'Deleting...' : 'Delete'}
          </Button>
        </DialogActions>
      </Dialog>
    </Container>
  );
}
