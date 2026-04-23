'use client';

import { use, useCallback, useEffect, useMemo, useState } from 'react';
import { useRouter } from 'next/navigation';
import {
  Alert,
  Box,
  Breadcrumbs,
  Button,
  CircularProgress,
  Container,
  Link,
  Paper,
  Stack,
  Typography,
} from '@mui/material';
import { ArrowBack, Save } from '@mui/icons-material';
import { useCompoundsApi } from '@/lib/compounds/api';
import { useAggregation } from '@/lib/compounds/aggregation-api';
import { routes } from '@/lib/compounds/routes';
import { ScorecardEditor } from '@/components/compounds/ScorecardEditor';
import { CompoundSpider } from '@/components/compounds/CompoundSpider';
import type {
  Protocol,
  ScorecardConfig,
  Target,
} from '@/types/compounds/models';
import type { CompactAggregationResponse } from '@/types/compounds/aggregation';

interface PageProps {
  params: Promise<{ id: string }>;
}

export default function TargetScorecardPage({ params }: PageProps) {
  const { id } = use(params);
  const router = useRouter();
  const api = useCompoundsApi();

  const { data: target, isLoading: targetLoading, mutate: refreshTarget } =
    api.get<Target>(`targets/${id}/`);
  const { data: protocols, isLoading: protocolsLoading } =
    api.get<Protocol[]>('protocols/');

  // Sample data for the live-preview spider — a single-target aggregation
  // with all molecular properties (so Lipinski axes evaluate) and geomean
  // values (so protocol/ratio/worst_of axes evaluate).
  const { data: sampleData } = useAggregation({
    predicates: { targets: [id] },
    output_format: 'compact',
    aggregations: ['geomean', 'count'],
    include_properties: ['molecular_weight', 'clogp', 'hbd', 'hba'],
  });
  const sampleCompound = useMemo(() => {
    const compact = sampleData as CompactAggregationResponse | undefined;
    return compact?.data?.[0] ?? null;
  }, [sampleData]);

  const [draft, setDraft] = useState<ScorecardConfig>({ axes: [] });
  const [saving, setSaving] = useState(false);
  const [saveError, setSaveError] = useState<string | null>(null);
  const [dirty, setDirty] = useState(false);

  // Seed the draft from the target once it loads.
  useEffect(() => {
    if (target) {
      setDraft(target.scorecard_config ?? { axes: [] });
      setDirty(false);
    }
  }, [target]);

  const handleDraftChange = useCallback((next: ScorecardConfig) => {
    setDraft(next);
    setDirty(true);
    setSaveError(null);
  }, []);

  const handleSave = useCallback(async () => {
    setSaving(true);
    setSaveError(null);
    try {
      await api.patch(`targets/${id}/`, { scorecard_config: draft });
      await refreshTarget();
      setDirty(false);
    } catch (e) {
      setSaveError(e instanceof Error ? e.message : 'Failed to save');
    } finally {
      setSaving(false);
    }
  }, [api, id, draft, refreshTarget]);

  const loading = targetLoading || protocolsLoading;

  return (
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <Breadcrumbs sx={{ mb: 2 }}>
        <Link
          component="button"
          variant="body2"
          underline="hover"
          onClick={() => router.push(routes.registry.targets())}
        >
          Targets
        </Link>
        <Link
          component="button"
          variant="body2"
          underline="hover"
          onClick={() => router.push(routes.registry.target(id))}
        >
          {target?.name ?? '…'}
        </Link>
        <Typography variant="body2">Scorecard</Typography>
      </Breadcrumbs>

      <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 2 }}>
        <Button
          startIcon={<ArrowBack />}
          onClick={() => router.push(routes.registry.target(id))}
          variant="text"
        >
          Back to target
        </Button>
        <Box sx={{ flex: 1 }} />
        <Button
          variant="contained"
          startIcon={saving ? <CircularProgress size={16} /> : <Save />}
          onClick={handleSave}
          disabled={saving || !dirty}
        >
          {saving ? 'Saving…' : dirty ? 'Save changes' : 'Saved'}
        </Button>
      </Box>

      <Typography variant="h4" gutterBottom>
        Scorecard configuration
      </Typography>
      <Typography variant="body2" color="text.secondary" sx={{ mb: 3 }}>
        Define the spider-plot axes for{' '}
        <strong>{target?.name ?? 'this target'}</strong>. Each axis maps a
        compound&apos;s measurements to a single score, normalised against the
        axis&apos;s excellent / poor anchors.
      </Typography>

      {loading ? (
        <Stack direction="row" spacing={1} alignItems="center">
          <CircularProgress size={20} />
          <Typography variant="body2" color="text.secondary">
            Loading target and protocols…
          </Typography>
        </Stack>
      ) : !target ? (
        <Alert severity="error">Target not found.</Alert>
      ) : (
        <Box sx={{ display: 'flex', gap: 3, alignItems: 'flex-start', flexWrap: 'wrap' }}>
          <Paper sx={{ p: 3, flex: 1, minWidth: 360 }}>
            <ScorecardEditor
              value={draft}
              protocols={protocols ?? []}
              onChange={handleDraftChange}
              error={saveError}
              disabled={saving}
            />
          </Paper>

          <Paper sx={{ p: 3, minWidth: 400, position: 'sticky', top: 16 }}>
            <Typography variant="subtitle1" fontWeight={600} gutterBottom>
              Live preview
            </Typography>
            {sampleCompound ? (
              <>
                <Typography
                  variant="caption"
                  color="text.secondary"
                  sx={{ display: 'block', mb: 1 }}
                >
                  Sample compound: <strong>{sampleCompound.formatted_id}</strong>
                </Typography>
                <CompoundSpider config={draft} compound={sampleCompound} size="large" />
              </>
            ) : (
              <Typography variant="caption" color="text.secondary">
                No compound data available for this target yet — the preview
                will appear here once the target has at least one assay.
              </Typography>
            )}
          </Paper>
        </Box>
      )}
    </Container>
  );
}
