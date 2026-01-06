'use client';

import { useState, useEffect, useCallback } from 'react';
import { useRouter, useSearchParams } from 'next/navigation';
import {
  Container,
  Typography,
  Box,
  Paper,
  TextField,
  Button,
  Autocomplete,
  Chip,
  CircularProgress,
  InputAdornment,
  Alert,
} from '@mui/material';
import { Search, Medication, Science, TableChart, Clear } from '@mui/icons-material';
import Link from 'next/link';
import { Breadcrumbs } from '@/components/Breadcrumbs';
import { DataTable, Column } from '@/components/DataTable';
import { MoleculeChip } from '@/components/MoleculeView';
import { useCompoundsApi } from '@/lib/api';
import { Target, Compound } from '@/types/models';

export default function CompoundSearchPage() {
  const router = useRouter();
  const searchParams = useSearchParams();
  const api = useCompoundsApi();

  // URL params
  const initialQuery = searchParams.get('q') || '';
  const initialTarget = searchParams.get('target') || '';

  // Search state
  const [searchQuery, setSearchQuery] = useState(initialQuery);
  const [debouncedQuery, setDebouncedQuery] = useState(initialQuery);
  const [selectedTarget, setSelectedTarget] = useState<Target | null>(null);
  const [loading, setLoading] = useState(false);
  const [results, setResults] = useState<Compound[]>([]);
  const [searched, setSearched] = useState(false);

  // Target options
  const { data: targets } = api.get<Target[]>('targets/');

  // Load initial target from URL
  useEffect(() => {
    if (initialTarget && targets) {
      const target = targets.find((t) => t.id === initialTarget);
      if (target) setSelectedTarget(target);
    }
  }, [initialTarget, targets]);

  // Debounce search query
  useEffect(() => {
    const timer = setTimeout(() => {
      setDebouncedQuery(searchQuery);
    }, 300);
    return () => clearTimeout(timer);
  }, [searchQuery]);

  // Perform search
  const doSearch = useCallback(async () => {
    if (!debouncedQuery.trim() && !selectedTarget) {
      setResults([]);
      setSearched(false);
      return;
    }

    setLoading(true);
    setSearched(true);

    try {
      // Build query params
      const params = new URLSearchParams();
      if (debouncedQuery.trim()) {
        params.set('search', debouncedQuery.trim());
      }
      if (selectedTarget) {
        params.set('target', selectedTarget.id);
      }

      const response = await fetch(`/api/proxy/compounds/compounds/?${params}`);
      if (response.ok) {
        const data = await response.json();
        setResults(data);
      } else {
        setResults([]);
      }
    } catch (error) {
      console.error('Search failed:', error);
      setResults([]);
    } finally {
      setLoading(false);
    }
  }, [debouncedQuery, selectedTarget]);

  // Auto-search when query changes
  useEffect(() => {
    doSearch();
  }, [doSearch]);

  // Update URL with search params
  useEffect(() => {
    const params = new URLSearchParams();
    if (debouncedQuery.trim()) params.set('q', debouncedQuery.trim());
    if (selectedTarget) params.set('target', selectedTarget.id);

    const newUrl = params.toString()
      ? `/registry/search?${params}`
      : '/registry/search';

    window.history.replaceState({}, '', newUrl);
  }, [debouncedQuery, selectedTarget]);

  const clearSearch = () => {
    setSearchQuery('');
    setSelectedTarget(null);
    setResults([]);
    setSearched(false);
  };

  const columns: Column<Compound>[] = [
    {
      key: 'formatted_id',
      label: 'ID',
      sortable: true,
      searchable: true,
      width: 180,
      render: (value, row) => (
        <Box sx={{ display: 'flex', alignItems: 'flex-start', gap: 1 }}>
          <Medication fontSize="small" color="secondary" sx={{ mt: 0.3 }} />
          <Box>
            <Typography fontWeight={500} fontFamily="monospace">
              {value}
            </Typography>
            {row.supplier_ref && (
              <Typography
                variant="caption"
                color="text.secondary"
                sx={{ display: 'block', fontFamily: 'monospace' }}
              >
                {row.supplier_ref}
              </Typography>
            )}
          </Box>
        </Box>
      ),
    },
    {
      key: 'smiles',
      label: 'Structure',
      width: 100,
      render: (value) => <MoleculeChip smiles={value} size={80} />,
    },
    {
      key: 'target_name',
      label: 'Target',
      sortable: true,
      width: 120,
      render: (value) => value || '-',
    },
    {
      key: 'molecular_weight',
      label: 'MW',
      sortable: true,
      width: 80,
      render: (value) => (value ? value.toFixed(1) : '-'),
    },
    {
      key: 'stereo_comment',
      label: 'Stereo',
      sortable: true,
      width: 100,
      render: (value) =>
        value && value !== 'unset' ? (
          <Chip label={value} size="small" variant="outlined" />
        ) : (
          '-'
        ),
    },
    {
      key: 'id',
      label: 'Actions',
      width: 140,
      render: (_, row) => (
        <Button
          component={Link}
          href={`/assays/aggregate?compound=${row.formatted_id}`}
          size="small"
          startIcon={<TableChart />}
          onClick={(e) => e.stopPropagation()}
        >
          Assay Data
        </Button>
      ),
    },
  ];

  return (
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <Breadcrumbs
        items={[
          { label: 'Home', href: '/', icon: 'home' },
          { label: 'Registry', href: '/registry/targets' },
          { label: 'Search', icon: 'search' },
        ]}
      />

      <Box sx={{ mb: 3 }}>
        <Typography variant="h4" gutterBottom>
          Compound Search
        </Typography>
        <Typography color="text.secondary">
          Search compounds by ID, supplier reference, or structure
        </Typography>
      </Box>

      {/* Search form */}
      <Paper sx={{ p: 3, mb: 3 }}>
        <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
          <Box sx={{ display: 'flex', gap: 2, flexWrap: 'wrap' }}>
            <TextField
              label="Search"
              placeholder="NCL-00026..., supplier ref, or SMILES"
              value={searchQuery}
              onChange={(e) => setSearchQuery(e.target.value)}
              sx={{ flex: 1, minWidth: 300 }}
              InputProps={{
                startAdornment: (
                  <InputAdornment position="start">
                    <Search />
                  </InputAdornment>
                ),
                endAdornment: loading ? (
                  <InputAdornment position="end">
                    <CircularProgress size={20} />
                  </InputAdornment>
                ) : null,
              }}
            />
            <Autocomplete
              options={targets || []}
              value={selectedTarget}
              onChange={(_, newValue) => setSelectedTarget(newValue)}
              getOptionLabel={(option) => option.name}
              isOptionEqualToValue={(option, value) => option.id === value.id}
              sx={{ minWidth: 200 }}
              renderInput={(params) => (
                <TextField
                  {...params}
                  label="Filter by Target"
                  InputProps={{
                    ...params.InputProps,
                    startAdornment: (
                      <>
                        <Science sx={{ mr: 1, color: 'action.active' }} />
                        {params.InputProps.startAdornment}
                      </>
                    ),
                  }}
                />
              )}
            />
            <Button
              variant="outlined"
              onClick={clearSearch}
              startIcon={<Clear />}
              disabled={!searchQuery && !selectedTarget}
            >
              Clear
            </Button>
          </Box>

          {selectedTarget && (
            <Box sx={{ display: 'flex', gap: 1 }}>
              <Chip
                icon={<Science fontSize="small" />}
                label={selectedTarget.name}
                onDelete={() => setSelectedTarget(null)}
                color="primary"
                variant="outlined"
              />
            </Box>
          )}
        </Box>
      </Paper>

      {/* Results */}
      {searched && !loading && results.length === 0 && (
        <Alert severity="info" sx={{ mb: 3 }}>
          No compounds found matching your search criteria.
        </Alert>
      )}

      {results.length > 0 && (
        <>
          <Box sx={{ mb: 2, display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
            <Typography variant="h6">
              {results.length} compound{results.length !== 1 ? 's' : ''} found
            </Typography>
            {results.length > 0 && (
              <Button
                component={Link}
                href={`/assays/aggregate?compound=${results.map(r => r.formatted_id).join(',')}`}
                variant="contained"
                startIcon={<TableChart />}
              >
                View All Assay Data
              </Button>
            )}
          </Box>

          <DataTable
            data={results}
            columns={columns}
            loading={loading}
            onRowClick={(compound) =>
              router.push(`/registry/compounds/${compound.id}`)
            }
            getRowKey={(row) => row.id}
            emptyMessage="No compounds found"
            additionalSearchFields={['supplier_ref', 'supplier_name', 'barcode', 'comments']}
          />
        </>
      )}

      {/* Help text when no search */}
      {!searched && !loading && (
        <Paper sx={{ p: 4, textAlign: 'center' }}>
          <Search sx={{ fontSize: 64, color: 'action.disabled', mb: 2 }} />
          <Typography variant="h6" color="text.secondary" gutterBottom>
            Search for compounds
          </Typography>
          <Typography color="text.secondary">
            Enter a compound ID (e.g., NCL-00026123), supplier reference, or SMILES string.
            <br />
            You can also filter by target to narrow your search.
          </Typography>
        </Paper>
      )}
    </Container>
  );
}
