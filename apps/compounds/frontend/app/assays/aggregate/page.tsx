'use client';

import { useState, useCallback, useRef, Suspense, useEffect } from 'react';
import { useSearchParams, useRouter, usePathname } from 'next/navigation';
import { Typography, Box, Alert, CircularProgress, IconButton, Tooltip, Snackbar } from '@mui/material';
import { TableChart, Link as LinkIcon, Save } from '@mui/icons-material';
import { DetailPageLayout } from '@/components/compounds/DetailPageLayout';
import { routes } from '@/lib/compounds/routes';
import { PredicateBuilder, PredicateBuilderState } from '@/components/compounds/PredicateBuilder';
import { AggregationTable } from '@/components/compounds/AggregationTable';
import {
  Predicates,
  AggregationType,
  OutputFormat,
  AggregationResponse,
  ConcentrationDisplayMode,
  MolecularPropertyName,
} from '@/types/compounds/aggregation';
import { fetchAggregation, saveAggregationView } from '@/lib/compounds/aggregation-api';
import { useAuth } from '@/lib/compounds/auth-context';

export default function AggregationPage() {
  return (
    <Suspense fallback={<Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100vh' }}><CircularProgress /></Box>}>
      <AggregationPageContent />
    </Suspense>
  );
}

function AggregationPageContent() {
  const searchParams = useSearchParams();
  const router = useRouter();
  const pathname = usePathname();
  const { canAdminister } = useAuth();

  // Read initial values from URL params for deep linking
  // Support both single 'compound' param and multiple 'compounds' params
  const initialCompoundSearch = searchParams.get('compound') || searchParams.get('compounds') || undefined;
  const initialTargetId = searchParams.get('target') || undefined;
  // Support multiple target names via comma-separated 'targets' param
  const initialTargetNames = searchParams.get('targets')?.split(',').map((s) => s.trim()).filter(Boolean) || undefined;
  // Support multiple protocol names via comma-separated 'protocols' param
  const initialProtocolNames = searchParams.get('protocols')?.split(',').map((s) => s.trim()).filter(Boolean) || undefined;
  // Support output format via 'format' param
  const formatParam = searchParams.get('format');
  const initialOutputFormat = (formatParam === 'compact' || formatParam === 'medium' || formatParam === 'long' || formatParam === 'pivot' || formatParam === 'cards') ? formatParam : undefined;
  // Support aggregations via comma-separated 'aggregations' param
  const aggregationsParam = searchParams.get('aggregations');
  const validAggregations = ['geomean', 'count', 'stdev', 'list'] as const;
  const initialAggregations = aggregationsParam
    ? aggregationsParam.split(',').map((s) => s.trim()).filter((a): a is AggregationType => validAggregations.includes(a as any))
    : undefined;
  // Support status filter via 'status' param (simplified to 'valid' or '' for any)
  const statusParam = searchParams.get('status');
  const initialStatus = (statusParam === 'valid' || statusParam === '') ? statusParam : undefined;
  // Support concentration display mode via 'concentrationDisplay' param
  const concentrationDisplayParam = searchParams.get('concentrationDisplay');
  const initialConcentrationDisplay = (concentrationDisplayParam === 'natural' || concentrationDisplayParam === 'nM' || concentrationDisplayParam === 'uM' || concentrationDisplayParam === 'mM' || concentrationDisplayParam === 'pConc')
    ? concentrationDisplayParam as ConcentrationDisplayMode
    : undefined;

  // Determine if this is a fresh aggregation (no query filters) - should start expanded
  const isFreshAggregation = !initialCompoundSearch && !initialTargetId && !initialTargetNames?.length && !initialProtocolNames?.length;
  // Support group by batch via 'groupByBatch' param
  const groupByBatchParam = searchParams.get('groupByBatch');
  const initialGroupByBatch = groupByBatchParam === 'true';
  // Support include tested no data via 'includeTestedNoData' param
  const includeTestedNoDataParam = searchParams.get('includeTestedNoData');
  const initialIncludeTestedNoData = includeTestedNoDataParam === 'true';
  // Support molecular properties via comma-separated 'properties' param
  const propertiesParam = searchParams.get('properties');
  const validProperties = ['molecular_weight', 'heavy_atom_count', 'hbd', 'hba', 'clogp', 'tpsa', 'rotatable_bonds', 'fraction_sp3'] as const;
  const initialIncludeProperties = propertiesParam
    ? propertiesParam.split(',').map((s) => s.trim()).filter((p): p is MolecularPropertyName => validProperties.includes(p as any))
    : undefined;

  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [data, setData] = useState<AggregationResponse | null>(null);
  const [currentAggregations, setCurrentAggregations] = useState<AggregationType[]>(initialAggregations || ['geomean', 'count']);
  const [currentOutputFormat, setCurrentOutputFormat] = useState<OutputFormat>(initialOutputFormat || 'compact');
  const [currentGroupByBatch, setCurrentGroupByBatch] = useState<boolean>(initialGroupByBatch);
  const [currentIncludeTestedNoData, setCurrentIncludeTestedNoData] = useState<boolean>(initialIncludeTestedNoData);
  const [currentIncludeProperties, setCurrentIncludeProperties] = useState<MolecularPropertyName[]>(initialIncludeProperties || []);
  const [currentState, setCurrentState] = useState<PredicateBuilderState | null>(null);
  const [concentrationDisplay, setConcentrationDisplay] = useState<ConcentrationDisplayMode>(initialConcentrationDisplay || 'natural');
  const [snackbarOpen, setSnackbarOpen] = useState(false);
  const [snackbarMessage, setSnackbarMessage] = useState('');
  const [saving, setSaving] = useState(false);

  // Track current request to avoid race conditions
  const requestIdRef = useRef(0);

  const handleStateChange = useCallback((state: PredicateBuilderState) => {
    setCurrentState(state);
  }, []);

  const handleCopyLink = useCallback(() => {
    if (!currentState) return;

    const params = new URLSearchParams();
    if (currentState.targetNames.length > 0) {
      params.set('targets', currentState.targetNames.join(','));
    }
    if (currentState.protocolNames.length > 0) {
      params.set('protocols', currentState.protocolNames.join(','));
    }
    if (currentState.compoundSearch) {
      params.set('compound', currentState.compoundSearch);
    }
    if (currentState.outputFormat !== 'compact') {
      params.set('format', currentState.outputFormat);
    }
    if (currentState.aggregations.length > 0 && !(currentState.aggregations.length === 2 && currentState.aggregations.includes('geomean') && currentState.aggregations.includes('count'))) {
      params.set('aggregations', currentState.aggregations.join(','));
    }
    if (currentState.status !== undefined && currentState.status !== 'valid') {
      params.set('status', currentState.status);
    }
    if (concentrationDisplay !== 'natural') {
      params.set('concentrationDisplay', concentrationDisplay);
    }
    if (currentState.groupByBatch) {
      params.set('groupByBatch', 'true');
    }
    if (currentState.includeTestedNoData) {
      params.set('includeTestedNoData', 'true');
    }
    if (currentState.includeProperties.length > 0) {
      params.set('properties', currentState.includeProperties.join(','));
    }

    const url = `${window.location.origin}${window.location.pathname}${params.toString() ? '?' + params.toString() : ''}`;
    navigator.clipboard.writeText(url).then(() => {
      setSnackbarMessage('Link copied to clipboard');
      setSnackbarOpen(true);
    });
  }, [currentState, concentrationDisplay]);

  const handleSaveToTarget = useCallback(async () => {
    if (!currentState || currentState.targets.length !== 1) return;

    const target = currentState.targets[0];
    setSaving(true);

    try {
      await saveAggregationView(target.id, {
        protocol_names: currentState.protocolNames,
        compound_search: currentState.compoundSearch,
        output_format: currentState.outputFormat,
        aggregations: currentState.aggregations,
        status: 'valid', // Default to valid status
        concentration_display: concentrationDisplay,
        include_properties: currentState.includeProperties,
      });
      setSnackbarMessage(`Saved view to ${target.name}`);
      setSnackbarOpen(true);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to save view');
    } finally {
      setSaving(false);
    }
  }, [currentState, concentrationDisplay]);

  // Update URL to persist query state for back navigation
  const updateUrlState = useCallback((state: PredicateBuilderState, outputFormat: OutputFormat, aggregations: AggregationType[], groupByBatch: boolean, includeTestedNoData: boolean, includeProperties: MolecularPropertyName[]) => {
    const params = new URLSearchParams();
    if (state.targetNames.length > 0) {
      params.set('targets', state.targetNames.join(','));
    }
    if (state.protocolNames.length > 0) {
      params.set('protocols', state.protocolNames.join(','));
    }
    if (state.compoundSearch) {
      params.set('compound', state.compoundSearch);
    }
    if (outputFormat !== 'compact') {
      params.set('format', outputFormat);
    }
    if (aggregations.length > 0 && !(aggregations.length === 2 && aggregations.includes('geomean') && aggregations.includes('count'))) {
      // Only include aggregations if different from default
      params.set('aggregations', aggregations.join(','));
    }
    // Always include status in URL for complete state preservation
    if (state.status !== undefined) {
      params.set('status', state.status);
    }
    if (concentrationDisplay !== 'natural') {
      params.set('concentrationDisplay', concentrationDisplay);
    }
    if (groupByBatch) {
      params.set('groupByBatch', 'true');
    }
    if (includeTestedNoData) {
      params.set('includeTestedNoData', 'true');
    }
    if (includeProperties.length > 0) {
      params.set('properties', includeProperties.join(','));
    }

    const queryString = params.toString();
    const newUrl = `${pathname}${queryString ? '?' + queryString : ''}`;
    router.replace(newUrl, { scroll: false });
  }, [pathname, router, concentrationDisplay]);

  // Update URL when concentration display changes (if we have data)
  useEffect(() => {
    if (data && currentState) {
      updateUrlState(currentState, currentState.outputFormat, currentAggregations, currentGroupByBatch, currentIncludeTestedNoData, currentIncludeProperties);
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [concentrationDisplay]);

  const handleChange = useCallback(async (
    predicates: Predicates,
    outputFormat: OutputFormat,
    aggregations: AggregationType[],
    groupByBatch: boolean,
    includeTestedNoData: boolean,
    includeProperties: MolecularPropertyName[]
  ) => {
    // Increment request ID to track latest request
    const requestId = ++requestIdRef.current;

    setLoading(true);
    setError(null);
    setCurrentAggregations(aggregations);
    setCurrentOutputFormat(outputFormat);
    setCurrentGroupByBatch(groupByBatch);
    setCurrentIncludeTestedNoData(includeTestedNoData);
    setCurrentIncludeProperties(includeProperties);

    // Update URL to persist state for back navigation
    if (currentState) {
      updateUrlState(currentState, outputFormat, aggregations, groupByBatch, includeTestedNoData, includeProperties);
    }

    try {
      // Map pivot and cards to compact for the API - they use the same data structure
      const apiOutputFormat = (outputFormat === 'pivot' || outputFormat === 'cards') ? 'compact' : outputFormat;
      const result = await fetchAggregation({
        predicates,
        output_format: apiOutputFormat,
        aggregations,
        group_by_batch: groupByBatch,
        include_tested_no_data: includeTestedNoData,
        include_properties: includeProperties.length > 0 ? includeProperties : undefined,
      });

      // Only update if this is still the latest request
      if (requestId === requestIdRef.current) {
        setData(result);
      }
    } catch (err) {
      // Only update error if this is still the latest request
      if (requestId === requestIdRef.current) {
        setError(err instanceof Error ? err.message : 'Failed to run aggregation');
        setData(null);
      }
    } finally {
      // Only clear loading if this is still the latest request
      if (requestId === requestIdRef.current) {
        setLoading(false);
      }
    }
  }, [currentState, updateUrlState]);

  // Build summary description from current state
  const getSummaryDescription = () => {
    const parts: string[] = [];
    if (currentState?.targetNames.length) {
      parts.push(`${currentState.targetNames.length} target${currentState.targetNames.length > 1 ? 's' : ''}`);
    }
    if (currentState?.protocolNames.length) {
      parts.push(`${currentState.protocolNames.length} protocol${currentState.protocolNames.length > 1 ? 's' : ''}`);
    }
    if (currentState?.compoundSearch) {
      parts.push(`"${currentState.compoundSearch}"`);
    }
    return parts.length > 0 ? parts.join(', ') : 'No filters';
  };

  // Summary configuration for collapsed header
  const summary = {
    title: 'Data Aggregation',
    titleIcon: <TableChart sx={{ fontSize: 'inherit' }} />,
    fields: [
      { label: 'Query', value: getSummaryDescription() },
      { label: 'Results', value: data ? `${data.data?.length || 0} rows` : loading ? 'Loading...' : '-' },
      { label: 'Format', value: currentOutputFormat },
    ],
    actions: (
      <>
        <Tooltip title="Copy shareable link">
          <IconButton
            onClick={handleCopyLink}
            disabled={!currentState || (!currentState.targetNames.length && !currentState.protocolNames.length && !currentState.compoundSearch)}
            size="small"
          >
            <LinkIcon />
          </IconButton>
        </Tooltip>
        {canAdminister && (
          <Tooltip title={currentState?.targets.length === 1 ? `Save view to ${currentState.targets[0].name}` : 'Select exactly one target to save'}>
            <span>
              <IconButton
                onClick={handleSaveToTarget}
                disabled={saving || !currentState || currentState.targets.length !== 1}
                size="small"
              >
                <Save />
              </IconButton>
            </span>
          </Tooltip>
        )}
      </>
    ),
  };

  // Detail content: the query builder
  const detailContent = (
    <Box>
      <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
        Aggregate KPI values across compounds and protocols
      </Typography>

      <PredicateBuilder
        onRunQuery={handleChange}
        onStateChange={handleStateChange}
        loading={loading}
        initialCompoundSearch={initialCompoundSearch}
        initialTargetId={initialTargetId}
        initialTargetNames={initialTargetNames}
        initialProtocolNames={initialProtocolNames}
        initialOutputFormat={initialOutputFormat}
        initialAggregations={initialAggregations}
        initialStatus={initialStatus}
        initialGroupByBatch={initialGroupByBatch}
        initialIncludeTestedNoData={initialIncludeTestedNoData}
        initialIncludeProperties={initialIncludeProperties}
      />

      {error && (
        <Alert severity="error" sx={{ mt: 2 }}>
          {error}
        </Alert>
      )}
    </Box>
  );

  return (
    <>
      <DetailPageLayout
        breadcrumbs={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Assays', href: routes.assays.list(), icon: 'assay' },
          { label: 'Data Aggregation', icon: 'aggregate' },
        ]}
        summary={summary}
        detailContent={detailContent}
        defaultCollapsed={!isFreshAggregation}
      >
        <AggregationTable
          data={data}
          loading={loading}
          aggregations={currentAggregations}
          outputFormat={currentOutputFormat}
          concentrationDisplay={concentrationDisplay}
          onConcentrationDisplayChange={setConcentrationDisplay}
          fillHeight
        />
      </DetailPageLayout>

      <Snackbar
        open={snackbarOpen}
        autoHideDuration={2000}
        onClose={() => setSnackbarOpen(false)}
        message={snackbarMessage}
      />
    </>
  );
}
