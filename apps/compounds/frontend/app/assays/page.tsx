'use client';

import { useRouter } from 'next/navigation';
import Link from 'next/link';
import { Container, Typography, Box, Chip, Button, Tooltip } from '@mui/material';
import { Assessment, Science, Description, Upload, Functions, TableChart, Biotech, FiberNew } from '@mui/icons-material';
import { PageHeader } from '@/components/compounds/PageHeader';
import { DataTable, Column } from '@/components/data-table';
import { useCompoundsApi } from '@/lib/compounds/api';
import { routes } from '@/lib/compounds/routes';
import { Assay } from '@/types/compounds/models';

/**
 * Strip Django's uniqueifying suffix from a filename.
 * Django adds a suffix like `_fm29oOG` (underscore + 7 alphanumeric chars) before
 * the extension when saving files with duplicate names.
 *
 * Examples:
 * - "2026_01_19_B_fm29oOG.xlsx" → "2026_01_19_B.xlsx"
 * - "data_abc1234.csv" → "data.csv"
 * - "report.pdf" → "report.pdf" (no suffix)
 */
function stripDjangoSuffix(filename: string): string {
  // Match: underscore + exactly 7 alphanumeric chars + extension
  // The pattern is conservative to avoid stripping legitimate filename parts
  const match = filename.match(/^(.+)_[a-zA-Z0-9]{7}(\.[^.]+)$/);
  if (match) {
    return match[1] + match[2];
  }
  return filename;
}

/**
 * Get a display-friendly version of a filename.
 * Strips Django's uniqueifying suffix and optionally truncates.
 */
function getDisplayFilename(
  filename: string | null | undefined,
  maxLength: number = 30
): { display: string; full: string; truncated: boolean } {
  if (!filename) {
    return { display: 'No file', full: '', truncated: false };
  }

  const cleaned = stripDjangoSuffix(filename);
  const truncated = cleaned.length > maxLength;
  const display = truncated
    ? cleaned.slice(0, maxLength - 3) + '...'
    : cleaned;

  return { display, full: filename, truncated: truncated || cleaned !== filename };
}

export default function AssaysPage() {
  const router = useRouter();
  const api = useCompoundsApi();
  const { data: assays, isLoading } = api.get<Assay[]>('assays/');

  // Helper to check if assay is new (created in last 7 days)
  const isRecentAssay = (createdAt: string | undefined) => {
    if (!createdAt) return false;
    const sevenDaysAgo = new Date();
    sevenDaysAgo.setDate(sevenDaysAgo.getDate() - 7);
    return new Date(createdAt) >= sevenDaysAgo;
  };

  const columns: Column<Assay>[] = [
    {
      key: 'data_filename',
      label: 'Data File',
      sortable: true,
      searchable: true,
      width: 220,
      render: (value, row) => {
        const { display, full, truncated } = getDisplayFilename(value as string, 25);
        const content = (
          <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, minWidth: 0 }}>
            <Assessment fontSize="small" color="info" sx={{ flexShrink: 0 }} />
            <Typography
              fontWeight={500}
              sx={{
                overflow: 'hidden',
                textOverflow: 'ellipsis',
                whiteSpace: 'nowrap',
              }}
            >
              {display}
            </Typography>
            {isRecentAssay(row.created_at) && (
              <Tooltip title="New in the last 7 days">
                <FiberNew fontSize="small" color="secondary" sx={{ flexShrink: 0 }} />
              </Tooltip>
            )}
          </Box>
        );

        // Wrap in tooltip if filename was cleaned or truncated
        return truncated ? (
          <Tooltip title={full} placement="top-start">
            {content}
          </Tooltip>
        ) : (
          content
        );
      },
    },
    {
      key: 'protocol_name',
      label: 'Protocol',
      sortable: true,
      searchable: true,
      hiddenOnMobile: true,
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
      hiddenOnMobile: true,
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
      hiddenOnMobile: true,
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
      hiddenOnMobile: true,
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
    <Box sx={{ display: 'flex', flexDirection: 'column', height: '100vh', overflow: 'hidden' }}>
      <Container maxWidth="lg" sx={{ flex: 1, display: 'flex', flexDirection: 'column', overflow: 'hidden', py: 2 }}>
        <Box sx={{ flexShrink: 0 }}>
          <PageHeader
            breadcrumbs={[
              { label: 'Home', href: routes.home(), icon: 'home' },
              { label: 'Assays', icon: 'assay' },
            ]}
          />

          <Box sx={{ mb: 2, display: 'flex', alignItems: 'flex-start', gap: 2 }}>
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
                href={routes.assays.importTableOfValues()}
                variant="outlined"
                startIcon={<TableChart />}
              >
                Import Table of Values
              </Button>
              <Button
                component={Link}
                href={routes.assays.importAdme()}
                variant="outlined"
                startIcon={<Biotech />}
              >
                Import ADME
              </Button>
              <Button
                component={Link}
                href={routes.assays.import()}
                variant="contained"
                startIcon={<Upload />}
              >
                Import Plate Data
              </Button>
            </Box>
          </Box>
        </Box>

        <Box sx={{ flex: 1, overflow: 'hidden', minHeight: 0 }}>
          <DataTable
            data={assays}
            columns={columns}
            loading={isLoading}
            onRowClick={(assay) => router.push(routes.assays.detail(assay.id))}
            getRowKey={(row) => row.id}
            title={assays ? `${assays.length} assays` : undefined}
            emptyMessage="No assays found"
            fillHeight
          />
        </Box>
      </Container>
    </Box>
  );
}
