'use client';

import { useState, useCallback, useRef, Suspense } from 'react';
import { useSearchParams } from 'next/navigation';
import { Container, Typography, Box, Alert, CircularProgress } from '@mui/material';
import { TableChart } from '@mui/icons-material';
import { Breadcrumbs } from '@/components/compounds/Breadcrumbs';
import { routes } from '@/lib/compounds/routes';
import { PredicateBuilder } from '@/components/compounds/PredicateBuilder';
import { AggregationTable } from '@/components/compounds/AggregationTable';
import {
  Predicates,
  AggregationType,
  OutputFormat,
  AggregationResponse,
} from '@/types/compounds/aggregation';
import { fetchAggregation } from '@/lib/compounds/aggregation-api';

export default function AggregationPage() {
  return (
    <Suspense fallback={<Container maxWidth="xl" sx={{ py: 3 }}><CircularProgress /></Container>}>
      <AggregationPageContent />
    </Suspense>
  );
}

function AggregationPageContent() {
  const searchParams = useSearchParams();

  // Read initial values from URL params for deep linking
  const initialCompoundSearch = searchParams.get('compound') || undefined;
  const initialTargetId = searchParams.get('target') || undefined;

  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [data, setData] = useState<AggregationResponse | null>(null);
  const [currentAggregations, setCurrentAggregations] = useState<AggregationType[]>(['geomean', 'count']);

  // Track current request to avoid race conditions
  const requestIdRef = useRef(0);

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
    <Container maxWidth="xl" sx={{ py: 3 }}>
      <Breadcrumbs
        items={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Assays', href: routes.assays.list(), icon: 'assay' },
          { label: 'Data Aggregation', icon: 'aggregate' },
        ]}
      />

      <Box sx={{ mb: 3, display: 'flex', alignItems: 'center', gap: 2 }}>
        <TableChart sx={{ fontSize: 40, color: 'primary.main' }} />
        <Box>
          <Typography variant="h4" gutterBottom sx={{ mb: 0 }}>
            Data Aggregation
          </Typography>
          <Typography color="text.secondary">
            Aggregate KPI values across compounds and protocols
          </Typography>
        </Box>
      </Box>

      <PredicateBuilder
        onChange={handleChange}
        loading={loading}
        initialCompoundSearch={initialCompoundSearch}
        initialTargetId={initialTargetId}
      />

      {error && (
        <Alert severity="error" sx={{ mb: 3 }}>
          {error}
        </Alert>
      )}

      <AggregationTable data={data} loading={loading} aggregations={currentAggregations} />
    </Container>
  );
}
