'use client';

import { use, useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { usePathname, useRouter, useSearchParams } from 'next/navigation';
import {
  Alert,
  Autocomplete,
  Box,
  Breadcrumbs,
  Button,
  CircularProgress,
  Container,
  IconButton,
  Link,
  Paper,
  Stack,
  TextField,
  Tooltip,
  Typography,
} from '@mui/material';
import { ArrowBack, Check, ContentCopy, Save } from '@mui/icons-material';
import html2canvas from 'html2canvas';
import { pinCaptureFonts } from '@/lib/compounds/html2canvas-fonts';
import { useCompoundsApi } from '@/lib/compounds/api';
import { useAggregation } from '@/lib/compounds/aggregation-api';
import { routes } from '@/lib/compounds/routes';
import { ScorecardEditor } from '@/components/compounds/ScorecardEditor';
import { CompoundSpider } from '@/components/compounds/CompoundSpider';
import { ScorecardValuesTable } from '@/components/compounds/ScorecardValuesTable';
import { ScorecardKeyTable } from '@/components/compounds/ScorecardKeyTable';
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
  const pathname = usePathname();
  const searchParams = useSearchParams();
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
  const sampleCompounds = useMemo(() => {
    const compact = sampleData as CompactAggregationResponse | undefined;
    return compact?.data ?? [];
  }, [sampleData]);

  // Optional URL-encoded compound choice — `?compound=NCL-00031070` pins
  // the live preview to a specific compound for shareable bookmarks.
  // When unset (or unmatched), default to the first compound returned by
  // the aggregation so something always renders.
  const pickedFormattedId = searchParams?.get('compound') ?? null;
  const sampleCompound = useMemo(() => {
    if (sampleCompounds.length === 0) return null;
    if (pickedFormattedId) {
      const match = sampleCompounds.find((c) => c.formatted_id === pickedFormattedId);
      if (match) return match;
    }
    return sampleCompounds[0];
  }, [sampleCompounds, pickedFormattedId]);

  const handlePickCompound = useCallback(
    (formattedId: string | null) => {
      const params = new URLSearchParams(searchParams?.toString() ?? '');
      if (formattedId) params.set('compound', formattedId);
      else params.delete('compound');
      const qs = params.toString();
      router.replace(`${pathname}${qs ? '?' + qs : ''}`, { scroll: false });
    },
    [pathname, router, searchParams],
  );

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
            {sampleCompounds.length > 0 ? (
              <>
                <Autocomplete
                  size="small"
                  options={sampleCompounds}
                  value={sampleCompound}
                  getOptionLabel={(c) => c.formatted_id ?? ''}
                  isOptionEqualToValue={(a, b) => a.compound_id === b.compound_id}
                  onChange={(_, picked) =>
                    handlePickCompound(picked?.formatted_id ?? null)
                  }
                  renderInput={(params) => (
                    <TextField
                      {...params}
                      label="Sample compound"
                      helperText="Pick which compound's data to render in the preview"
                    />
                  )}
                  sx={{ mb: 2 }}
                />
                {sampleCompound && (
                  <>
                    <CompoundSpider config={draft} compound={sampleCompound} size="large" />
                    <Box sx={{ mt: 2 }}>
                      <ScorecardValuesTable
                        config={draft}
                        compound={sampleCompound}
                        protocols={(sampleData as CompactAggregationResponse | undefined)?.protocols}
                        caption={`Values for ${sampleCompound.formatted_id}`}
                      />
                    </Box>
                  </>
                )}
              </>
            ) : (
              <Typography variant="caption" color="text.secondary">
                No compound data available for this target yet — the preview
                will appear here once the target has at least one assay.
              </Typography>
            )}
          </Paper>

          <ScorecardKeyPanel
            config={draft}
            protocols={(sampleData as CompactAggregationResponse | undefined)?.protocols}
          />
        </Box>
      )}
    </Container>
  );
}

// ---------------------------------------------------------------------------
// Compact "key" panel — explains each axis's formula + thresholds in
// presentation-ready form. Sits in the editor's right column under the
// live preview; copy-as-image button lets chemists drop it onto a
// PowerPoint slide alongside spider cards. Sized to fit ~5 spider cards
// across a slide.
// ---------------------------------------------------------------------------

function ScorecardKeyPanel({
  config,
  protocols,
}: {
  config: ScorecardConfig;
  protocols?: import('@/types/compounds/aggregation').ProtocolInfo[];
}) {
  const keyRef = useRef<HTMLDivElement>(null);
  const [copied, setCopied] = useState(false);

  const handleCopy = useCallback(() => {
    const node = keyRef.current;
    if (!node) return;
    // Safari and recent Chrome both invalidate the clipboard user-gesture
    // chain if the page awaits before calling navigator.clipboard.write.
    // html2canvas easily takes >100ms, so awaiting the canvas first
    // triggers NotAllowedError. Workaround: hand a Promise<Blob> directly
    // to ClipboardItem — the spec recognises this as a continuation of
    // the original click and the browsers honour it.
    const blobPromise = (async () => {
      // Wait for web fonts to finish loading — html2canvas otherwise
      // captures with system-font fallback metrics, which produces
      // mis-spaced text (e.g. "Lipinski compliance" overlapping).
      if (typeof document !== 'undefined' && document.fonts?.ready) {
        await document.fonts.ready;
      }
      const canvas = await html2canvas(node, {
        backgroundColor: '#ffffff',
        scale: 2, // higher DPI for crisp paste at slide scale
        onclone: pinCaptureFonts,
      });
      return new Promise<Blob>((resolve, reject) => {
        canvas.toBlob((blob) => {
          if (blob) resolve(blob);
          else reject(new Error('Could not convert canvas to PNG'));
        }, 'image/png');
      });
    })();

    navigator.clipboard
      .write([new ClipboardItem({ 'image/png': blobPromise })])
      .then(() => {
        setCopied(true);
        setTimeout(() => setCopied(false), 1800);
      })
      .catch((err) => {
        console.error('Failed to copy key to clipboard:', err);
      });
  }, []);

  if (config.axes.length === 0) return null;

  return (
    <Paper sx={{ p: 2, minWidth: 480, maxWidth: 800 }}>
      <Box sx={{ display: 'flex', alignItems: 'center', mb: 1 }}>
        <Typography variant="subtitle1" fontWeight={600} sx={{ flex: 1 }}>
          Scorecard key
        </Typography>
        <Tooltip title={copied ? 'Copied!' : 'Copy as image'} arrow>
          <IconButton onClick={handleCopy} size="small">
            {copied ? <Check fontSize="small" color="success" /> : <ContentCopy fontSize="small" />}
          </IconButton>
        </Tooltip>
      </Box>
      <Typography
        variant="caption"
        color="text.secondary"
        sx={{ display: 'block', mb: 1.5 }}
      >
        Drop this beside your spider cards on a slide. ★ excellent → ✗ poor.
      </Typography>
      <Box
        ref={keyRef}
        sx={{
          // White background and a hairline frame so the captured PNG
          // stands alone visually when pasted into a slide.
          bgcolor: '#fff',
          border: '1px solid',
          borderColor: 'divider',
          borderRadius: 1,
          p: 1.5,
        }}
      >
        <ScorecardKeyTable config={config} protocols={protocols} />
      </Box>
    </Paper>
  );
}

