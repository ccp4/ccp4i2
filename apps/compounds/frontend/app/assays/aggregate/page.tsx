'use client';

import { useState, useCallback, useRef, Suspense } from 'react';
import { useSearchParams } from 'next/navigation';
import { Container, Typography, Box, Alert, CircularProgress, IconButton, Tooltip, Snackbar } from '@mui/material';
import { TableChart, Link as LinkIcon, Save } from '@mui/icons-material';
import { PageHeader } from '@/components/compounds/PageHeader';
import { routes } from '@/lib/compounds/routes';
import { PredicateBuilder, PredicateBuilderState } from '@/components/compounds/PredicateBuilder';
import { AggregationTable } from '@/components/compounds/AggregationTable';
import {
  Predicates,
  AggregationType,
  OutputFormat,
  AggregationResponse,
  ConcentrationDisplayMode,
} from '@/types/compounds/aggregation';
import { fetchAggregation, saveAggregationView } from '@/lib/compounds/aggregation-api';
import { useAuth } from '@/lib/compounds/auth-context';

export default function AggregationPage() {
  return (
    <Suspense fallback={<Container maxWidth="xl" sx={{ py: 3 }}><CircularProgress /></Container>}>
      <AggregationPageContent />
    </Suspense>
  );
}

function AggregationPageContent() {
  const searchParams = useSearchParams();
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
  const initialOutputFormat = (formatParam === 'compact' || formatParam === 'medium' || formatParam === 'long') ? formatParam : undefined;

  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [data, setData] = useState<AggregationResponse | null>(null);
  const [currentAggregations, setCurrentAggregations] = useState<AggregationType[]>(['geomean', 'count']);
  const [currentState, setCurrentState] = useState<PredicateBuilderState | null>(null);
  const [concentrationDisplay, setConcentrationDisplay] = useState<ConcentrationDisplayMode>('natural');
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

    const url = `${window.location.origin}${window.location.pathname}${params.toString() ? '?' + params.toString() : ''}`;
    navigator.clipboard.writeText(url).then(() => {
      setSnackbarMessage('Link copied to clipboard');
      setSnackbarOpen(true);
    });
  }, [currentState]);

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
      });
      setSnackbarMessage(`Saved view to ${target.name}`);
      setSnackbarOpen(true);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to save view');
    } finally {
      setSaving(false);
    }
  }, [currentState, concentrationDisplay]);

  const handleChange = useCallback(async (
    predicates: Predicates,
    outputFormat: OutputFormat,
    aggregations: AggregationType[]
  ) => {
    // Increment request ID to track latest request
    const requestId = ++requestIdRef.current;

    setLoading(true);
    setError(null);
    setCurrentAggregations(aggregations);

    try {
      const result = await fetchAggregation({
        predicates,
        output_format: outputFormat,
        aggregations,
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
  }, []);

  return (
    <Container maxWidth="xl" sx={{ py: 2 }}>
      <PageHeader
        breadcrumbs={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Assays', href: routes.assays.list(), icon: 'assay' },
          { label: 'Data Aggregation', icon: 'aggregate' },
        ]}
      />

      <Box sx={{ mb: 2, display: 'flex', alignItems: 'center', gap: 2 }}>
        <TableChart sx={{ fontSize: 32, color: 'primary.main' }} />
        <Box sx={{ flex: 1 }}>
          <Typography variant="h5" sx={{ mb: 0 }}>
            Data Aggregation
          </Typography>
          <Typography variant="body2" color="text.secondary">
            Aggregate KPI values across compounds and protocols
          </Typography>
        </Box>
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
      </Box>

      <PredicateBuilder
        onRunQuery={handleChange}
        onStateChange={handleStateChange}
        loading={loading}
        initialCompoundSearch={initialCompoundSearch}
        initialTargetId={initialTargetId}
        initialTargetNames={initialTargetNames}
        initialProtocolNames={initialProtocolNames}
        initialOutputFormat={initialOutputFormat}
      />

      {error && (
        <Alert severity="error" sx={{ mb: 3 }}>
          {error}
        </Alert>
      )}

      <AggregationTable
        data={data}
        loading={loading}
        aggregations={currentAggregations}
        concentrationDisplay={concentrationDisplay}
        onConcentrationDisplayChange={setConcentrationDisplay}
      />

      <Snackbar
        open={snackbarOpen}
        autoHideDuration={2000}
        onClose={() => setSnackbarOpen(false)}
        message={snackbarMessage}
      />
    </Container>
  );
}
