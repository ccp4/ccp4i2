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
  Link as MuiLink,
} from '@mui/material';
import {
  Assessment,
  Science,
  Description,
  Medication,
  CheckCircle,
  Cancel,
  HelpOutline,
} from '@mui/icons-material';
import { Breadcrumbs } from '@/components/Breadcrumbs';
import { DataTable, Column } from '@/components/DataTable';
import { DoseResponseThumb } from '@/components/DoseResponseChart';
import { CompoundStructureCell } from '@/components/CompoundStructureCell';
import { useCompoundsApi } from '@/lib/api';
import { Assay, DataSeries, Protocol, Target } from '@/types/models';

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

  // Helper to get chart data from a data series
  const getChartData = (row: DataSeries) => {
    if (!row.dilution_series?.concentrations || !row.extracted_data) return null;
    return {
      concentrations: row.dilution_series.concentrations,
      responses: Array.isArray(row.extracted_data) ? row.extracted_data : [],
      unit: row.dilution_series.unit,
    };
  };

  const columns: Column<DataSeries>[] = [
    {
      key: 'chart',
      label: 'Curve',
      width: 80,
      render: (_, row) => {
        const chartData = getChartData(row);
        if (!chartData || chartData.concentrations.length === 0) {
          return <Box sx={{ width: 60, height: 60, bgcolor: 'grey.100', borderRadius: 1 }} />;
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
      width: 90,
      render: (_, row) => (
        <CompoundStructureCell
          compoundId={row.compound}
          compoundName={row.compound_name}
          size={80}
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
              if (row.compound) router.push(`/registry/compounds/${row.compound}`);
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
          { label: 'Home', href: '/', icon: 'home' },
          { label: 'Protocols', href: '/assays/protocols', icon: 'protocol' },
          ...(protocol ? [{ label: protocol.name, href: `/assays/protocols/${protocol.id}`, icon: 'protocol' as const }] : []),
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
              <Box>
                <Typography variant="h4">{assay.data_filename || 'Assay'}</Typography>
                <Box sx={{ display: 'flex', gap: 1, mt: 0.5, flexWrap: 'wrap' }}>
                  {protocol && (
                    <Chip
                      icon={<Description fontSize="small" />}
                      label={protocol.name}
                      size="small"
                      onClick={() => router.push(`/assays/protocols/${protocol.id}`)}
                    />
                  )}
                  {target && (
                    <Chip
                      icon={<Science fontSize="small" />}
                      label={target.name}
                      size="small"
                      variant="outlined"
                      onClick={() => router.push(`/registry/targets/${target.id}`)}
                    />
                  )}
                </Box>
              </Box>
            </Box>

            <Divider sx={{ my: 2 }} />

            <Grid container spacing={3}>
              <Grid item xs={12} md={6}>
                <Typography variant="h6" gutterBottom>
                  Experiment Details
                </Typography>
                <InfoRow label="Protocol" value={protocol?.name} />
                <InfoRow label="Target" value={target?.name} />
                <InfoRow label="Lab Book" value={assay.labbook_number} />
                <InfoRow label="Page" value={assay.page_number} />
              </Grid>
              <Grid item xs={12} md={6}>
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
        onRowClick={(row) => router.push(`/assays/data-series/${row.id}`)}
      />
    </Container>
  );
}
