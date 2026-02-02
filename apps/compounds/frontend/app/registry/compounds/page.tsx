'use client';

import { Suspense, useState, useMemo } from 'react';
import { useRouter, useSearchParams } from 'next/navigation';
import Link from 'next/link';
import {
  Container,
  Typography,
  Box,
  Chip,
  Button,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  SelectChangeEvent,
  Skeleton,
  Tooltip,
} from '@mui/material';
import {
  Science,
  Add,
  FilterList,
  Apps,
  Error as ErrorIcon,
} from '@mui/icons-material';
import { PageHeader } from '@/components/compounds/PageHeader';
import { DataTable, Column } from '@/components/data-table';
import { MoleculeChip } from '@/components/compounds/MoleculeView';
import { useCompoundsApi } from '@/lib/compounds/api';
import { useAuth } from '@/lib/compounds/auth-context';
import { routes } from '@/lib/compounds/routes';
import { Compound, Target, PaginatedResponse } from '@/types/compounds/models';

function CompoundsPageContent() {
  const router = useRouter();
  const searchParams = useSearchParams();
  const api = useCompoundsApi();
  const { canContribute } = useAuth();

  // Get initial target filter from URL
  const initialTarget = searchParams.get('target') || '';
  const [targetFilter, setTargetFilter] = useState(initialTarget);

  // Fetch targets for filter dropdown
  const { data: targets, error: targetsError } = api.get<Target[]>('targets/');

  // Build API URL with target filter
  const compoundsUrl = useMemo(() => {
    const params = new URLSearchParams();
    if (targetFilter) {
      params.set('target', targetFilter);
    }
    const queryString = params.toString();
    return `compounds/${queryString ? `?${queryString}` : ''}`;
  }, [targetFilter]);

  // Fetch compounds with optional target filter
  // Backend may return either paginated response (with results array) or plain array
  const { data: compoundsResponse, isLoading, error: compoundsError } = api.get<PaginatedResponse<Compound> | Compound[]>(compoundsUrl);
  const compounds = Array.isArray(compoundsResponse)
    ? compoundsResponse
    : (compoundsResponse?.results || []);

  // Combine errors for display - identify which endpoint failed
  const error = compoundsError || targetsError;
  const errorSource = compoundsError ? 'compounds' : targetsError ? 'targets' : null;

  // Log error for debugging - this helps identify if the large response is failing
  if (error) {
    console.error(`[Compounds Page] API error (${errorSource}):`, error);
  }

  const handleTargetChange = (event: SelectChangeEvent) => {
    const newTarget = event.target.value;
    setTargetFilter(newTarget);

    // Update URL to reflect filter state
    const newUrl = newTarget
      ? routes.registry.compounds({ target: newTarget })
      : routes.registry.compounds();
    router.replace(newUrl, { scroll: false });
  };

  const columns: Column<Compound>[] = [
    {
      key: 'formatted_id',
      label: 'Compound ID',
      sortable: true,
      searchable: true,
      width: 140,
      render: (value) => (
        <Typography fontWeight={600} fontFamily="monospace" color="primary.main">
          {value}
        </Typography>
      ),
    },
    {
      key: 'smiles',
      label: 'Structure',
      width: 125,
      render: (value) => value ? <MoleculeChip smiles={value} size={100} /> : '-',
    },
    {
      key: 'target_name',
      label: 'Target',
      sortable: true,
      searchable: true,
      width: 120,
      hiddenOnMobile: true,
      render: (value, row) => (
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 0.5 }}>
          <Science fontSize="small" color="action" />
          <Typography variant="body2">{value || row.target}</Typography>
        </Box>
      ),
    },
    {
      key: 'supplier_name',
      label: 'Supplier',
      sortable: true,
      searchable: true,
      hiddenOnMobile: true,
      render: (value, row) => value || row.supplier_ref || '-',
    },
    {
      key: 'molecular_weight',
      label: 'MW',
      sortable: true,
      width: 90,
      hiddenOnMobile: true,
      render: (value) => value ? value.toFixed(1) : '-',
    },
    {
      key: 'batch_count',
      label: 'Batches',
      sortable: true,
      width: 90,
      hiddenOnMobile: true,
      render: (value) =>
        value ? (
          <Chip label={value} size="small" variant="outlined" />
        ) : (
          '-'
        ),
    },
    {
      key: 'registered_at',
      label: 'Registered',
      sortable: true,
      width: 110,
      render: (value) =>
        value ? new Date(value).toLocaleDateString() : '-',
    },
  ];

  // Get selected target name for breadcrumb
  const selectedTargetName = targetFilter
    ? targets?.find(t => t.id === targetFilter)?.name
    : null;

  return (
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <PageHeader
        breadcrumbs={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Registry', href: routes.registry.targets() },
          { label: 'Compounds', icon: 'compound' },
          ...(selectedTargetName ? [{ label: selectedTargetName }] : []),
        ]}
      />

      <Box sx={{ mb: 3, display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', flexWrap: 'wrap', gap: 2 }}>
        <Box>
          <Typography variant="h4" gutterBottom>
            Compounds
          </Typography>
          <Typography color="text.secondary">
            {compounds.length > 0
              ? `${compounds.length} compounds${selectedTargetName ? ` for ${selectedTargetName}` : ' registered'}`
              : 'Browse registered compounds'}
          </Typography>
        </Box>
        <Box sx={{ display: 'flex', gap: 1, alignItems: 'center', flexWrap: 'wrap' }}>
          {/* Target filter */}
          <FormControl size="small" sx={{ minWidth: 180 }}>
            <InputLabel id="target-filter-label">
              <Box sx={{ display: 'flex', alignItems: 'center', gap: 0.5 }}>
                <FilterList fontSize="small" />
                Target
              </Box>
            </InputLabel>
            <Select
              labelId="target-filter-label"
              value={targetFilter}
              label="Target"
              onChange={handleTargetChange}
            >
              <MenuItem value="">
                <em>All Targets</em>
              </MenuItem>
              {targets?.map((target) => (
                <MenuItem key={target.id} value={target.id}>
                  {target.name}
                </MenuItem>
              ))}
            </Select>
          </FormControl>

          <Button
            component={Link}
            href="/"
            variant="outlined"
            startIcon={<Apps />}
          >
            All Apps
          </Button>

          <Tooltip title={canContribute ? '' : 'Requires Contributor or Admin operating level'} arrow>
            <span>
              <Button
                component={Link}
                href={targetFilter ? `${routes.registry.new()}?target=${targetFilter}` : routes.registry.new()}
                variant="contained"
                startIcon={<Add />}
                disabled={!canContribute}
              >
                Register Compound
              </Button>
            </span>
          </Tooltip>
        </Box>
      </Box>

      {error && (
        <Box sx={{ mb: 2, p: 2, bgcolor: 'error.light', borderRadius: 1, display: 'flex', alignItems: 'center', gap: 1 }}>
          <ErrorIcon color="error" />
          <Typography color="error.dark">
            Failed to load {errorSource}: {error.message || 'Unknown error'}
          </Typography>
        </Box>
      )}

      <DataTable
        data={compounds}
        columns={columns}
        loading={isLoading}
        onRowClick={(compound) => router.push(routes.registry.compound(compound.id))}
        getRowKey={(row) => row.id}
        emptyMessage={targetFilter ? "No compounds found for this target" : "No compounds registered yet"}
        comfortable
      />
    </Container>
  );
}

function CompoundsPageFallback() {
  return (
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <Skeleton variant="rectangular" height={40} sx={{ mb: 2 }} />
      <Skeleton variant="rectangular" height={60} sx={{ mb: 3 }} />
      <Skeleton variant="rectangular" height={400} />
    </Container>
  );
}

export default function CompoundsPage() {
  return (
    <Suspense fallback={<CompoundsPageFallback />}>
      <CompoundsPageContent />
    </Suspense>
  );
}
