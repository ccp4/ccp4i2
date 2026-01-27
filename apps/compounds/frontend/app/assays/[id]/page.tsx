'use client';

import { use, useState, useCallback } from 'react';
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
  Tooltip,
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
  Edit,
  Image as ImageIcon,
  Palette,
  Refresh,
} from '@mui/icons-material';
import { PageHeader } from '@/components/compounds/PageHeader';
import { DataTable, Column } from '@/components/compounds/DataTable';
import { DoseResponseThumb } from '@/components/compounds/DoseResponseChart';
import { CompoundStructureCell } from '@/components/compounds/CompoundStructureCell';
import { ImageBatchUpload } from '@/components/compounds/ImageBatchUpload';
import { PlateHeatMapDialog } from '@/components/compounds/PlateHeatMap';
import { AssayEditDialog } from '@/components/compounds/AssayEditDialog';
import { AuthenticatedImage } from '@/components/compounds/AuthenticatedImage';
import { useCompoundsApi, getAuthenticatedDownloadUrl, authFetch } from '@/lib/compounds/api';
import { formatKpiUnit } from '@/lib/compounds/aggregation-api';
import { useAuth } from '@/lib/compounds/auth-context';
import { routes } from '@/lib/compounds/routes';
import { Assay, DataSeries, Protocol, Target, PlateLayout } from '@/types/compounds/models';

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
  const { canContribute } = useAuth();
  const [deleteDialogOpen, setDeleteDialogOpen] = useState(false);
  const [deleting, setDeleting] = useState(false);
  const [editDialogOpen, setEditDialogOpen] = useState(false);
  const [imageUploadOpen, setImageUploadOpen] = useState(false);
  const [heatMapOpen, setHeatMapOpen] = useState(false);
  const [heatMapCells, setHeatMapCells] = useState<(string | number | null)[][] | null>(null);
  const [heatMapLoading, setHeatMapLoading] = useState(false);
  const [reanalysing, setReanalysing] = useState(false);

  const { data: assay, isLoading: assayLoading, mutate: mutateAssay } = api.get<Assay>(
    `assays/${id}/`
  );
  const { data: dataSeries, isLoading: dataSeriesLoading, mutate: mutateDataSeries } = api.get<DataSeries[]>(
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
      await api.delete(`assays/${id}/`);
      // Navigate back to protocol page or assays list
      if (protocol) {
        router.push(routes.assays.protocol(protocol.id));
      } else {
        router.push(routes.assays.list());
      }
    } catch (err) {
      console.error('Delete error:', err);
    } finally {
      setDeleting(false);
      setDeleteDialogOpen(false);
    }
  };

  const handleReanalyseAll = async () => {
    setReanalysing(true);
    try {
      await api.post(`assays/${id}/analyse_all/`, {});
      // Refresh data series to show updated analysis results
      mutateDataSeries();
    } catch (err) {
      console.error('Re-analyse error:', err);
    } finally {
      setReanalysing(false);
    }
  };

  const handleOpenHeatMap = useCallback(async () => {
    if (!assay?.data_file) return;

    // If we already have the cells loaded, just open the dialog
    if (heatMapCells) {
      setHeatMapOpen(true);
      return;
    }

    setHeatMapLoading(true);
    try {
      // Fetch the Excel file (through protected endpoint with auth)
      const response = await authFetch(assay.data_file);
      const buffer = await response.arrayBuffer();

      // Parse with XLSX
      const XLSX = await import('xlsx');
      const workbook = XLSX.read(buffer, { type: 'array' });
      const sheetName = workbook.SheetNames[0];
      const worksheet = workbook.Sheets[sheetName];

      if (!worksheet) {
        throw new Error('No worksheet found');
      }

      // Get sheet as 2D array
      const range = XLSX.utils.decode_range(worksheet['!ref'] || 'A1');
      const cells: (string | number | null)[][] = [];

      for (let row = range.s.r; row <= range.e.r; row++) {
        const rowData: (string | number | null)[] = [];
        for (let col = range.s.c; col <= range.e.c; col++) {
          const cellAddress = XLSX.utils.encode_cell({ r: row, c: col });
          const cell = worksheet[cellAddress];
          if (cell) {
            rowData.push(cell.v as string | number | null);
          } else {
            rowData.push(null);
          }
        }
        cells.push(rowData);
      }

      setHeatMapCells(cells);
      setHeatMapOpen(true);
    } catch (err) {
      console.error('Failed to load Excel file for heat map:', err);
    } finally {
      setHeatMapLoading(false);
    }
  }, [assay?.data_file, heatMapCells]);

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

  // Check if this is a table_of_values assay
  const isTableOfValues = protocol?.analysis_method === 'table_of_values';

  const columns: Column<DataSeries>[] = [
    {
      key: 'chart',
      label: isTableOfValues ? 'Plot' : 'Curve',
      width: 140,
      render: (_, row) => {
        // For table_of_values, show the plot image if available
        if (isTableOfValues) {
          const hasImageFile = row.analysis?.results?.['Image File'];
          const plotImageUrl = row.plot_image
            ? `/api/proxy/compounds/media/data-series/${row.id}/plot/`
            : null;

          if (hasImageFile && plotImageUrl) {
            return (
              <AuthenticatedImage
                src={plotImageUrl}
                alt={`Plot for ${row.compound_name || 'compound'}`}
                width={120}
                height={120}
                objectFit="contain"
                sx={{
                  borderRadius: 1,
                  bgcolor: 'background.paper',
                  border: '1px solid',
                  borderColor: 'divider',
                }}
              />
            );
          }
          // No image available
          return (
            <Box
              sx={{
                width: 120,
                height: 120,
                bgcolor: 'grey.100',
                borderRadius: 1,
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
              }}
            >
              <ImageIcon sx={{ color: 'grey.400', fontSize: 32 }} />
            </Box>
          );
        }

        // For dose-response assays, show the chart
        const chartData = getChartData(row);
        if (!chartData || chartData.concentrations.length === 0) {
          return <Box sx={{ width: 120, height: 120, bgcolor: 'grey.100', borderRadius: 1 }} />;
        }
        const fitParams = row.analysis?.results ? (() => {
          const results = row.analysis.results;
          // Use KPI field to find the correct ec50 value
          const kpiKey = results.KPI || results.kpi;
          const ec50Value = kpiKey
            ? (results[kpiKey] ?? results[kpiKey.toLowerCase?.()] ?? results[kpiKey.toUpperCase?.()])
            : null;
          return {
            ec50: ec50Value ?? results.EC50 ?? results.ec50 ?? results.IC50 ?? results.ic50 ?? null,
            hill: results.Hill ?? results.hill ?? results.hill_slope ?? null,
            minVal: results.minVal ?? results.bottom ?? null,
            maxVal: results.maxVal ?? results.top ?? null,
            status: row.analysis.status,
          };
        })() : undefined;
        return <DoseResponseThumb data={chartData} fit={fitParams} size={120} />;
      },
    },
    {
      key: 'structure',
      label: 'Structure',
      width: 125,
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
      width: 120,
      render: (value, row) => {
        if (value === null || value === undefined) return '-';
        // Get unit from analysis results or dilution series
        const unit = row.analysis?.results?.kpi_unit || row.dilution_series?.unit;
        const formattedUnit = formatKpiUnit(unit);
        return (
          <Typography fontFamily="monospace" fontWeight={500}>
            {typeof value === 'number' ? value.toFixed(2) : value}
            {formattedUnit && (
              <Typography
                component="span"
                color="text.secondary"
                sx={{ ml: 0.5, fontWeight: 400 }}
              >
                {formattedUnit}
              </Typography>
            )}
          </Typography>
        );
      },
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
      <PageHeader
        breadcrumbs={[
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
              <Box sx={{ display: 'flex', gap: 1 }}>
                {protocol?.plate_layout_config && assay?.data_file && (
                  <Button
                    variant="outlined"
                    startIcon={heatMapLoading ? <CircularProgress size={16} /> : <Palette />}
                    onClick={handleOpenHeatMap}
                    disabled={heatMapLoading}
                    size="small"
                  >
                    Heat Map
                  </Button>
                )}
                {protocol?.analysis_method === 'table_of_values' && (
                  <Tooltip title={canContribute ? '' : 'Requires Contributor or Admin operating level'} arrow>
                    <span>
                      <Button
                        variant="outlined"
                        startIcon={<ImageIcon />}
                        onClick={() => setImageUploadOpen(true)}
                        size="small"
                        disabled={!canContribute}
                      >
                        Upload Images
                      </Button>
                    </span>
                  </Tooltip>
                )}
                <Tooltip title={canContribute ? 'Re-analyse all data series' : 'Requires Contributor or Admin operating level'} arrow>
                  <span>
                    <Button
                      variant="outlined"
                      startIcon={reanalysing ? <CircularProgress size={16} /> : <Refresh />}
                      onClick={handleReanalyseAll}
                      size="small"
                      disabled={!canContribute || reanalysing}
                    >
                      {reanalysing ? 'Analysing...' : 'Re-analyse All'}
                    </Button>
                  </span>
                </Tooltip>
                <Tooltip title={canContribute ? '' : 'Requires Contributor or Admin operating level'} arrow>
                  <span>
                    <Button
                      variant="outlined"
                      startIcon={<Edit />}
                      onClick={() => setEditDialogOpen(true)}
                      size="small"
                      disabled={!canContribute}
                    >
                      Edit
                    </Button>
                  </span>
                </Tooltip>
                <Tooltip title={canContribute ? '' : 'Requires Contributor or Admin operating level'} arrow>
                  <span>
                    <Button
                      variant="outlined"
                      color="error"
                      startIcon={<Delete />}
                      onClick={() => setDeleteDialogOpen(true)}
                      size="small"
                      disabled={!canContribute}
                    >
                      Delete
                    </Button>
                  </span>
                </Tooltip>
              </Box>
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
                      <MuiLink
                        component="button"
                        onClick={async () => {
                          const url = await getAuthenticatedDownloadUrl(assay.data_file!);
                          window.open(url, '_blank');
                        }}
                        sx={{ cursor: 'pointer' }}
                      >
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

      {/* Image batch upload dialog for Table of Values */}
      <ImageBatchUpload
        open={imageUploadOpen}
        onClose={() => setImageUploadOpen(false)}
        assayId={id}
        onUploaded={() => mutateDataSeries()}
      />

      {/* Assay edit dialog */}
      {assay && (
        <AssayEditDialog
          open={editDialogOpen}
          onClose={() => setEditDialogOpen(false)}
          assay={assay}
          onSave={() => mutateAssay()}
        />
      )}

      {/* Plate Heat Map dialog */}
      {heatMapCells && protocol?.plate_layout_config && (
        <PlateHeatMapDialog
          open={heatMapOpen}
          onClose={() => setHeatMapOpen(false)}
          cells={heatMapCells}
          plateLayout={protocol.plate_layout_config}
        />
      )}
    </Container>
  );
}
