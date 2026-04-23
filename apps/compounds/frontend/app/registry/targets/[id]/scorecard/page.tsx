'use client';

import { use, useCallback, useEffect, useState } from 'react';
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
import { routes } from '@/lib/compounds/routes';
import { ScorecardEditor } from '@/components/compounds/ScorecardEditor';
import type {
  Protocol,
  ScorecardConfig,
  Target,
} from '@/types/compounds/models';

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
        <Paper sx={{ p: 3 }}>
          <ScorecardEditor
            value={draft}
            protocols={protocols ?? []}
            onChange={handleDraftChange}
            error={saveError}
            disabled={saving}
          />
        </Paper>
      )}
    </Container>
  );
}
