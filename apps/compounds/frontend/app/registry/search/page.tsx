'use client';

import { useState, useEffect, useCallback, Suspense } from 'react';
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
  Accordion,
  AccordionSummary,
  AccordionDetails,
  ToggleButton,
  ToggleButtonGroup,
  Tooltip,
} from '@mui/material';
import {
  Search,
  Medication,
  Science,
  TableChart,
  Clear,
  ExpandMore,
  Draw,
  AccountTree,
  Hub,
  Download,
  Add,
  Upload,
} from '@mui/icons-material';
import Link from 'next/link';
import { PageHeader } from '@/components/compounds/PageHeader';
import { DataTable, Column } from '@/components/compounds/DataTable';
import { MoleculeChip } from '@/components/compounds/MoleculeView';
import { JSMEEditor } from '@/components/compounds/JSMEEditor';
import { useCompoundsApi } from '@/lib/compounds/api';
import { routes } from '@/lib/compounds/routes';
import { Target, Compound } from '@/types/compounds/models';

// =============================================================================
// Authentication Integration
// =============================================================================

// For Docker integration: Try to import auth helpers from ccp4i2 client's auth-token
// Falls back to no-op for standalone development
let getAccessToken: () => Promise<string | null>;
let getUserEmail: () => string | null;

try {
  const authModule = require('../../../utils/auth-token');
  getAccessToken = authModule.getAccessToken;
  getUserEmail = authModule.getUserEmail || (() => null);
} catch {
  getAccessToken = async () => null;
  getUserEmail = () => null;
}

/**
 * Authenticated fetch wrapper
 */
async function authFetch(url: string, options: RequestInit = {}): Promise<Response> {
  const headers: Record<string, string> = {};

  const token = await getAccessToken();
  if (token) {
    headers['Authorization'] = `Bearer ${token}`;
  }

  const email = getUserEmail();
  if (email) {
    headers['X-User-Email'] = email;
  }

  if (options.headers) {
    if (options.headers instanceof Headers) {
      options.headers.forEach((value, key) => {
        headers[key] = value;
      });
    } else if (Array.isArray(options.headers)) {
      options.headers.forEach(([key, value]) => {
        headers[key] = value;
      });
    } else {
      Object.assign(headers, options.headers);
    }
  }

  return fetch(url, { ...options, headers });
}

type StructureSearchMode = 'substructure' | 'superstructure';

export default function CompoundSearchPage() {
  return (
    <Suspense fallback={<Container maxWidth="lg" sx={{ py: 3 }}><CircularProgress /></Container>}>
      <CompoundSearchContent />
    </Suspense>
  );
}

function CompoundSearchContent() {
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

  // Structure search state
  const [structureExpanded, setStructureExpanded] = useState(false);
  const [drawnSmiles, setDrawnSmiles] = useState('');
  const [initialSmiles, setInitialSmiles] = useState('');
  const [structureSearchMode, setStructureSearchMode] = useState<StructureSearchMode>('substructure');
  const [structureSearchLoading, setStructureSearchLoading] = useState(false);
  const [structureSearchError, setStructureSearchError] = useState<string | null>(null);

  // Compound lookup state for loading structures
  const [compoundLookupQuery, setCompoundLookupQuery] = useState('');
  const [compoundLookupResults, setCompoundLookupResults] = useState<Compound[]>([]);
  const [compoundLookupLoading, setCompoundLookupLoading] = useState(false);
  const [selectedLookupCompound, setSelectedLookupCompound] = useState<Compound | null>(null);

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

  // Debounced compound lookup for loading structures
  useEffect(() => {
    if (!compoundLookupQuery.trim()) {
      setCompoundLookupResults([]);
      return;
    }

    const timer = setTimeout(async () => {
      setCompoundLookupLoading(true);
      try {
        const params = new URLSearchParams();
        params.set('search', compoundLookupQuery.trim());
        const response = await authFetch(`/api/proxy/compounds/compounds/?${params}`);
        if (response.ok) {
          const data = await response.json();
          // Limit to first 20 results for the dropdown
          setCompoundLookupResults(data.slice(0, 20));
        }
      } catch (error) {
        console.error('Compound lookup failed:', error);
      } finally {
        setCompoundLookupLoading(false);
      }
    }, 300);

    return () => clearTimeout(timer);
  }, [compoundLookupQuery]);

  // Perform text search
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

      const response = await authFetch(`/api/proxy/compounds/compounds/?${params}`);
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
      ? `${routes.registry.search()}?${params}`
      : routes.registry.search();

    window.history.replaceState({}, '', newUrl);
  }, [debouncedQuery, selectedTarget]);

  // Perform structure search
  const doStructureSearch = useCallback(async () => {
    if (!drawnSmiles.trim()) {
      setStructureSearchError('Please draw a structure first');
      return;
    }

    setStructureSearchLoading(true);
    setStructureSearchError(null);
    setSearched(true);

    try {
      const params = new URLSearchParams();
      params.set('smiles', drawnSmiles.trim());
      params.set('mode', structureSearchMode);
      if (selectedTarget) {
        params.set('target', selectedTarget.id);
      }

      const response = await authFetch(`/api/proxy/compounds/compounds/structure_search/?${params}`);
      if (response.ok) {
        const data = await response.json();
        setResults(data.matches || []);
      } else {
        const errorData = await response.json().catch(() => ({}));
        setStructureSearchError(errorData.error || 'Structure search failed');
        setResults([]);
      }
    } catch (error) {
      console.error('Structure search failed:', error);
      setStructureSearchError('Structure search failed. Please try again.');
      setResults([]);
    } finally {
      setStructureSearchLoading(false);
    }
  }, [drawnSmiles, structureSearchMode, selectedTarget]);

  const clearSearch = () => {
    setSearchQuery('');
    setSelectedTarget(null);
    setResults([]);
    setSearched(false);
    setStructureSearchError(null);
  };

  // Load a compound's structure into the JSME editor
  const loadStructure = useCallback((compound: Compound | null) => {
    if (compound?.smiles) {
      setInitialSmiles(compound.smiles);
      setDrawnSmiles(compound.smiles);
      setStructureSearchError(null);
      // Expand the structure search panel if not already expanded
      setStructureExpanded(true);
    }
  }, []);

  const handleStructureChange = useCallback((smiles: string) => {
    setDrawnSmiles(smiles);
    setStructureSearchError(null);
  }, []);

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
          href={routes.assays.aggregate({ compound: row.formatted_id })}
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
      <PageHeader
        breadcrumbs={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Registry', href: routes.registry.targets() },
          { label: 'Search', icon: 'search' },
        ]}
      />

      <Box sx={{ mb: 3, display: 'flex', alignItems: 'flex-start', gap: 2 }}>
        <Box sx={{ flex: 1 }}>
          <Typography variant="h4" gutterBottom>
            Compound Search
          </Typography>
          <Typography color="text.secondary">
            Search compounds by ID, supplier reference, or structure
          </Typography>
        </Box>
        <Box sx={{ display: 'flex', gap: 1 }}>
          <Button
            component={Link}
            href={routes.registry.import()}
            variant="outlined"
            startIcon={<Upload />}
          >
            Bulk Import
          </Button>
          <Button
            component={Link}
            href={routes.registry.new()}
            variant="contained"
            startIcon={<Add />}
          >
            New Compound
          </Button>
        </Box>
      </Box>

      {/* Text Search form */}
      <Paper sx={{ p: 3, mb: 3 }}>
        <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
          <Box sx={{ display: 'flex', gap: 2, flexWrap: 'wrap' }}>
            <TextField
              label="Text Search"
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
              disabled={!searchQuery && !selectedTarget && !drawnSmiles}
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

      {/* Structure Search */}
      <Accordion
        expanded={structureExpanded}
        onChange={(_, expanded) => setStructureExpanded(expanded)}
        sx={{ mb: 3 }}
      >
        <AccordionSummary expandIcon={<ExpandMore />}>
          <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
            <Draw />
            <Typography fontWeight={500}>Structure Search</Typography>
            {drawnSmiles && (
              <Chip
                size="small"
                label="Structure drawn"
                color="success"
                variant="outlined"
              />
            )}
          </Box>
        </AccordionSummary>
        <AccordionDetails>
          <Box sx={{ display: 'flex', flexDirection: { xs: 'column', md: 'row' }, gap: 3 }}>
            {/* JSME Editor */}
            <Box sx={{ flexShrink: 0 }}>
              {/* Compound lookup to load structure */}
              <Box sx={{ mb: 2 }}>
                <Autocomplete
                  options={compoundLookupResults}
                  value={selectedLookupCompound}
                  onChange={(_, newValue) => {
                    setSelectedLookupCompound(newValue);
                    loadStructure(newValue);
                  }}
                  onInputChange={(_, value) => setCompoundLookupQuery(value)}
                  getOptionLabel={(option) => `${option.formatted_id}${option.supplier_ref ? ` (${option.supplier_ref})` : ''}`}
                  isOptionEqualToValue={(option, value) => option.id === value.id}
                  loading={compoundLookupLoading}
                  filterOptions={(x) => x} // Disable client-side filtering since we do server-side
                  renderOption={(props, option) => (
                    <Box component="li" {...props} key={option.id} sx={{ display: 'flex', gap: 1, alignItems: 'center' }}>
                      <MoleculeChip smiles={option.smiles} size={40} />
                      <Box>
                        <Typography variant="body2" fontFamily="monospace">
                          {option.formatted_id}
                        </Typography>
                        {option.supplier_ref && (
                          <Typography variant="caption" color="text.secondary">
                            {option.supplier_ref}
                          </Typography>
                        )}
                      </Box>
                    </Box>
                  )}
                  renderInput={(params) => (
                    <TextField
                      {...params}
                      label="Load structure from compound"
                      placeholder="Search by ID or supplier ref..."
                      size="small"
                      InputProps={{
                        ...params.InputProps,
                        startAdornment: (
                          <>
                            <Download sx={{ mr: 1, color: 'action.active' }} />
                            {params.InputProps.startAdornment}
                          </>
                        ),
                        endAdornment: (
                          <>
                            {compoundLookupLoading ? <CircularProgress size={16} /> : null}
                            {params.InputProps.endAdornment}
                          </>
                        ),
                      }}
                    />
                  )}
                  sx={{ width: 400 }}
                />
              </Box>

              <JSMEEditor
                id="search-jsme"
                onChange={handleStructureChange}
                editable={true}
                query={true}
                initialSmiles={initialSmiles}
                width={400}
                height={350}
                showPreview={true}
              />
            </Box>

            {/* Search controls */}
            <Box sx={{ flex: 1, display: 'flex', flexDirection: 'column', gap: 2 }}>
              <Typography variant="subtitle2" color="text.secondary">
                Search Mode
              </Typography>
              <ToggleButtonGroup
                value={structureSearchMode}
                exclusive
                onChange={(_, value) => value && setStructureSearchMode(value)}
                size="small"
              >
                <ToggleButton value="substructure">
                  <Tooltip title="Find compounds containing the drawn structure">
                    <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                      <Hub fontSize="small" />
                      Substructure
                    </Box>
                  </Tooltip>
                </ToggleButton>
                <ToggleButton value="superstructure">
                  <Tooltip title="Find compounds that are substructures of the drawn structure">
                    <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                      <AccountTree fontSize="small" />
                      Superstructure
                    </Box>
                  </Tooltip>
                </ToggleButton>
              </ToggleButtonGroup>

              <Typography variant="body2" color="text.secondary" sx={{ mt: 1 }}>
                {structureSearchMode === 'substructure'
                  ? 'Find all compounds that contain the drawn structure as a fragment (the drawn structure is a substructure of the results).'
                  : 'Find all compounds that are contained within the drawn structure (the results are substructures of the drawn structure).'}
              </Typography>

              {structureSearchError && (
                <Alert severity="error" sx={{ mt: 1 }}>
                  {structureSearchError}
                </Alert>
              )}

              <Button
                variant="contained"
                onClick={doStructureSearch}
                disabled={!drawnSmiles || structureSearchLoading}
                startIcon={structureSearchLoading ? <CircularProgress size={20} /> : <Search />}
                sx={{ mt: 2, alignSelf: 'flex-start' }}
              >
                {structureSearchLoading ? 'Searching...' : `Find ${structureSearchMode === 'substructure' ? 'Superstructures' : 'Substructures'}`}
              </Button>
            </Box>
          </Box>
        </AccordionDetails>
      </Accordion>

      {/* Results */}
      {searched && !loading && !structureSearchLoading && results.length === 0 && (
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
                href={routes.assays.aggregate({ compound: results.map(r => r.formatted_id).join(',') })}
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
            loading={loading || structureSearchLoading}
            onRowClick={(compound) =>
              router.push(routes.registry.compound(compound.id))
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
            You can also filter by target or use the structure editor to draw a query structure.
          </Typography>
        </Paper>
      )}
    </Container>
  );
}
