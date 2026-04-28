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
  CardContent,
  OutputFormat,
  AggregationResponse,
  ConcentrationDisplayMode,
  MolecularPropertyName,
} from '@/types/compounds/aggregation';
import { fetchAggregation, saveAggregationView } from '@/lib/compounds/aggregation-api';
import { scorecardDataNeeds } from '@/lib/compounds/scorecard';
import { useAuth } from '@/lib/compounds/auth-context';
import { useCompoundsApi } from '@/lib/compounds/api';
import { getSelection } from '@/lib/compounds/selections-api';
import type { CategorisationSelection } from '@/components/compounds/ProtocolScatterPlot';
import type { Target as TargetRecord } from '@/types/compounds/models';

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
  const api = useCompoundsApi();

  // Read initial values from URL params for deep linking
  // Support both single 'compound' param and multiple 'compounds' params,
  // plus the slice-20 token-addressed `selection` param (NLP redirects
  // through this when the compound list would exceed safe URL length).
  const inlineCompoundSearch = searchParams?.get('compound') || searchParams?.get('compounds') || undefined;
  const selectionToken = searchParams?.get('selection') || undefined;
  // Selection fetch — only fires when ?selection=<uuid> is in the URL.
  // The fetched compound_ids are joined into the same comma-separated
  // shape as the inline ?compound= form, so everything downstream
  // treats them identically.
  const { data: selectionData, error: selectionError } = api.get<{
    id: string;
    name: string;
    compound_ids: string[];
  }>(selectionToken ? `selections/${selectionToken}/` : null);
  const selectionPending = !!selectionToken && !selectionData && !selectionError;
  const initialCompoundSearch =
    inlineCompoundSearch ?? (selectionData?.compound_ids?.length
      ? selectionData.compound_ids.join(',')
      : undefined);
  const initialTargetId = searchParams?.get('target') || undefined;
  // Support multiple target names via comma-separated 'targets' param
  const initialTargetNames = searchParams?.get('targets')?.split(',').map((s) => s.trim()).filter(Boolean) || undefined;
  // Support multiple protocol names via comma-separated 'protocols' param
  const initialProtocolNames = searchParams?.get('protocols')?.split(',').map((s) => s.trim()).filter(Boolean) || undefined;
  // Support output format via 'format' param
  const formatParam = searchParams?.get('format');
  const initialOutputFormat = (formatParam === 'compact' || formatParam === 'medium' || formatParam === 'long' || formatParam === 'pivot' || formatParam === 'cards' || formatParam === 'bullets' || formatParam === 'scatter') ? formatParam : undefined;
  // Support aggregations via comma-separated 'aggregations' param
  const aggregationsParam = searchParams?.get('aggregations');
  const validAggregations = ['geomean', 'count', 'stdev', 'list'] as const;
  const initialAggregations = aggregationsParam
    ? aggregationsParam.split(',').map((s) => s.trim()).filter((a): a is AggregationType => validAggregations.includes(a as any))
    : undefined;
  // Support status filter via 'status' param (simplified to 'valid' or '' for any)
  const statusParam = searchParams?.get('status');
  const initialStatus = (statusParam === 'valid' || statusParam === '') ? statusParam : undefined;
  // Support concentration display mode via 'concentrationDisplay' param
  const concentrationDisplayParam = searchParams?.get('concentrationDisplay');
  const initialConcentrationDisplay = (concentrationDisplayParam === 'natural' || concentrationDisplayParam === 'nM' || concentrationDisplayParam === 'uM' || concentrationDisplayParam === 'mM' || concentrationDisplayParam === 'pConc')
    ? concentrationDisplayParam as ConcentrationDisplayMode
    : undefined;

  // Determine if this is a fresh aggregation (no query filters) - should start expanded
  const isFreshAggregation = !initialCompoundSearch && !initialTargetId && !initialTargetNames?.length && !initialProtocolNames?.length;
  // Support group by batch via 'groupByBatch' param
  const groupByBatchParam = searchParams?.get('groupByBatch');
  const initialGroupByBatch = groupByBatchParam === 'true';
  // Support include tested no data via 'includeTestedNoData' param
  const includeTestedNoDataParam = searchParams?.get('includeTestedNoData');
  const initialIncludeTestedNoData = includeTestedNoDataParam === 'true';
  // Support molecular properties via comma-separated 'properties' param
  const propertiesParam = searchParams?.get('properties');
  const validProperties = ['molecular_weight', 'heavy_atom_count', 'hbd', 'hba', 'clogp', 'tpsa', 'rotatable_bonds', 'fraction_sp3'] as const;
  const initialIncludeProperties = propertiesParam
    ? propertiesParam.split(',').map((s) => s.trim()).filter((p): p is MolecularPropertyName => validProperties.includes(p as any))
    : undefined;
  // Include identifiers defaults to true; explicit '0' disables it via URL.
  const identifiersParam = searchParams?.get('identifiers');
  const initialIncludeIdentifiers = identifiersParam === '0' ? false : true;
  // Cards-view body selector (Cards format only). Default chosen at render
  // time once we know whether the target has a scorecard, so we leave the
  // initial value undefined here and let CardsView pick a sensible default.
  const cardContentParam = searchParams?.get('cardContent');
  const initialCardContent: CardContent | undefined =
    (cardContentParam === 'protocols' || cardContentParam === 'spider' || cardContentParam === 'both' || cardContentParam === 'compact')
      ? cardContentParam
      : undefined;

  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [data, setData] = useState<AggregationResponse | null>(null);
  const [currentAggregations, setCurrentAggregations] = useState<AggregationType[]>(initialAggregations || ['geomean', 'count']);
  const [currentOutputFormat, setCurrentOutputFormat] = useState<OutputFormat>(initialOutputFormat || 'compact');
  const [currentGroupByBatch, setCurrentGroupByBatch] = useState<boolean>(initialGroupByBatch);
  const [currentIncludeTestedNoData, setCurrentIncludeTestedNoData] = useState<boolean>(initialIncludeTestedNoData);
  const [currentIncludeProperties, setCurrentIncludeProperties] = useState<MolecularPropertyName[]>(initialIncludeProperties || []);
  const [currentIncludeIdentifiers, setCurrentIncludeIdentifiers] = useState<boolean>(initialIncludeIdentifiers);
  const [cardContent, setCardContent] = useState<CardContent | undefined>(initialCardContent);
  const [currentState, setCurrentState] = useState<PredicateBuilderState | null>(null);

  // Categorisation overlay for scatter — derive directly from the URL
  // so the chip strip's router.replace flows through without extra
  // mirror state. Comma-joined string is the dep so React only refires
  // when the actual list changes (not on unrelated URL edits).
  const colourByIdsKey = searchParams?.get('colour_by') ?? '';
  const [categorisationSelections, setCategorisationSelections] = useState<CategorisationSelection[]>([]);
  useEffect(() => {
    const ids = colourByIdsKey.split(',').map((s) => s.trim()).filter(Boolean);
    if (ids.length === 0) {
      setCategorisationSelections([]);
      return;
    }
    let cancelled = false;
    Promise.all(ids.map(async (id) => {
      try {
        const sel = await getSelection(id);
        return { id: sel.id, name: sel.name, compoundIds: new Set(sel.compound_ids) } as CategorisationSelection;
      } catch {
        // 404 / 410 / network — drop silently; chip strip lets the
        // chemist re-pick. Surfacing every failed lookup as an error
        // would be noisy, especially after a Selection expired.
        return null;
      }
    })).then((rows) => {
      if (cancelled) return;
      setCategorisationSelections(rows.filter((r): r is CategorisationSelection => r !== null));
    });
    return () => { cancelled = true; };
  }, [colourByIdsKey]);

  // When exactly one target is selected, fetch its full record so we can pass
  // its scorecard_config down to the Cards view for per-card spider rendering.
  const singleTargetId =
    currentState?.targets.length === 1 ? currentState.targets[0].id : null;
  const { data: singleTargetRecord } = api.get<TargetRecord>(
    singleTargetId ? `targets/${singleTargetId}/` : null,
  );
  const singleTargetScorecard = singleTargetRecord?.scorecard_config ?? null;
  const [concentrationDisplay, setConcentrationDisplay] = useState<ConcentrationDisplayMode>(initialConcentrationDisplay || 'natural');
  const [snackbarOpen, setSnackbarOpen] = useState(false);
  const [snackbarMessage, setSnackbarMessage] = useState('');
  const [saving, setSaving] = useState(false);

  // Track current request to avoid race conditions
  const requestIdRef = useRef(0);

  // Track which scorecard was in scope at the time of the last fetch,
  // so we can detect when the SWR scorecard fetch resolves AFTER the
  // initial aggregation request has already gone out (without it, the
  // Lipinski axis arrives empty because molecular_weight / clogp / hbd
  // / hba weren't auto-injected into include_properties — see the
  // useEffect below that re-runs the query on transition).
  const requestScorecardRef = useRef<unknown>(undefined);

  // Cached args from the last handleChange call. Used when the SWR
  // scorecard fetch resolves AFTER the initial query has gone out and
  // we need to re-fire with the exact same predicates / format / etc.
  // — predicates are built inside PredicateBuilder and aren't directly
  // re-derivable here, so we stash a copy when the first call comes in.
  const lastFetchArgsRef = useRef<{
    predicates: Predicates;
    outputFormat: OutputFormat;
    aggregations: AggregationType[];
    groupByBatch: boolean;
    includeTestedNoData: boolean;
    includeProperties: MolecularPropertyName[];
    includeIdentifiers: boolean;
  } | null>(null);

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
    if (currentState.includeIdentifiers === false) {
      params.set('identifiers', '0');
    }
    // Cards body selector — only emit when in cards format and explicitly set,
    // so default-state share links stay short.
    if (currentState.outputFormat === 'cards' && cardContent) {
      params.set('cardContent', cardContent);
    }

    const url = `${window.location.origin}${window.location.pathname}${params.toString() ? '?' + params.toString() : ''}`;
    navigator.clipboard.writeText(url).then(() => {
      setSnackbarMessage('Link copied to clipboard');
      setSnackbarOpen(true);
    });
  }, [currentState, concentrationDisplay, cardContent]);

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
        include_identifiers: currentState.includeIdentifiers,
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
  const updateUrlState = useCallback((state: PredicateBuilderState, outputFormat: OutputFormat, aggregations: AggregationType[], groupByBatch: boolean, includeTestedNoData: boolean, includeProperties: MolecularPropertyName[], includeIdentifiers: boolean) => {
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
    if (includeIdentifiers === false) {
      params.set('identifiers', '0');
    }
    if (outputFormat === 'cards' && cardContent) {
      params.set('cardContent', cardContent);
    }

    const queryString = params.toString();
    const newUrl = `${pathname}${queryString ? '?' + queryString : ''}`;
    router.replace(newUrl, { scroll: false });
  }, [pathname, router, concentrationDisplay, cardContent]);

  // Update URL when concentration display changes (if we have data)
  useEffect(() => {
    if (data && currentState) {
      updateUrlState(currentState, currentState.outputFormat, currentAggregations, currentGroupByBatch, currentIncludeTestedNoData, currentIncludeProperties, currentIncludeIdentifiers);
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [concentrationDisplay]);

  // Update URL when the cards-view body selector changes
  useEffect(() => {
    if (data && currentState) {
      updateUrlState(currentState, currentState.outputFormat, currentAggregations, currentGroupByBatch, currentIncludeTestedNoData, currentIncludeProperties, currentIncludeIdentifiers);
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [cardContent]);

  const handleChange = useCallback(async (
    predicates: Predicates,
    outputFormat: OutputFormat,
    aggregations: AggregationType[],
    groupByBatch: boolean,
    includeTestedNoData: boolean,
    includeProperties: MolecularPropertyName[],
    includeIdentifiers: boolean
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
    setCurrentIncludeIdentifiers(includeIdentifiers);

    // Update URL to persist state for back navigation
    if (currentState) {
      updateUrlState(currentState, outputFormat, aggregations, groupByBatch, includeTestedNoData, includeProperties, includeIdentifiers);
    }

    try {
      // Map pivot/cards/bullets to compact for the API - they use the same data structure
      const apiOutputFormat = (outputFormat === 'pivot' || outputFormat === 'cards' || outputFormat === 'bullets' || outputFormat === 'scatter') ? 'compact' : outputFormat;

      // If a scorecard is in play for a single target, auto-augment the
      // request with the protocols and molecular properties the scorecard
      // axes need. This avoids the "my Lipinski axis is empty because I
      // didn't tick Include Properties" footgun.
      // Stamp the scorecard we used so the post-fetch effect can detect
      // a stale fetch when SWR resolves the scorecard later.
      requestScorecardRef.current = singleTargetScorecard;
      lastFetchArgsRef.current = {
        predicates,
        outputFormat,
        aggregations,
        groupByBatch,
        includeTestedNoData,
        includeProperties,
        includeIdentifiers,
      };
      const scorecardNeeds = scorecardDataNeeds(singleTargetScorecard);
      const effectivePredicates: Predicates = scorecardNeeds.protocolIds.length > 0
        ? {
            ...predicates,
            protocols: Array.from(
              new Set([...(predicates.protocols ?? []), ...scorecardNeeds.protocolIds]),
            ),
          }
        : predicates;
      const effectiveProperties = Array.from(
        new Set([...includeProperties, ...scorecardNeeds.properties]),
      ) as MolecularPropertyName[];

      const result = await fetchAggregation({
        predicates: effectivePredicates,
        output_format: apiOutputFormat,
        aggregations,
        group_by_batch: groupByBatch,
        include_tested_no_data: includeTestedNoData,
        include_properties: effectiveProperties.length > 0 ? effectiveProperties : undefined,
        include_identifiers: includeIdentifiers,
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
  }, [currentState, updateUrlState, singleTargetScorecard]);

  // SWR fetches the target's scorecard_config asynchronously. If
  // PredicateBuilder fired the initial query before the scorecard
  // resolved, the request went out without the scorecard's auto-
  // injected protocols/properties — most visibly, a Lipinski axis
  // shows "no data" because molecular_weight / clogp / hbd / hba
  // weren't requested. When the scorecard arrives, re-run the query
  // so the Cards view can populate those values.
  useEffect(() => {
    if (!data || !lastFetchArgsRef.current) return;
    if (requestScorecardRef.current === singleTargetScorecard) return;
    const args = lastFetchArgsRef.current;
    handleChange(
      args.predicates,
      args.outputFormat,
      args.aggregations,
      args.groupByBatch,
      args.includeTestedNoData,
      args.includeProperties,
      args.includeIdentifiers,
    );
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [singleTargetScorecard]);

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

  // Detail content: the query builder. When ?selection=<token> is in
  // the URL, gate the builder behind the selection fetch so it mounts
  // with the right initialCompoundSearch (slice 20 — the URL-length
  // fix for large NLP results).
  const detailContent = selectionPending ? (
    <Box sx={{ py: 4, display: 'flex', alignItems: 'center', justifyContent: 'center', gap: 2 }}>
      <CircularProgress size={20} />
      <Typography variant="body2" color="text.secondary">
        Loading selection…
      </Typography>
    </Box>
  ) : selectionError ? (
    <Alert severity="error" sx={{ mb: 2 }}>
      Couldn&apos;t load the selection — it may have expired or been
      created by a different user. Try re-running the original query.
    </Alert>
  ) : (
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
        initialIncludeIdentifiers={initialIncludeIdentifiers}
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
          scorecardConfig={singleTargetScorecard}
          cardContent={cardContent}
          onCardContentChange={setCardContent}
          fillHeight
          categorisationSelections={categorisationSelections}
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
