'use client';

/**
 * Chip strip mounted above the inline scatter view: each chip is a
 * chemotype (substructure) from the curated catalog OR an ad-hoc
 * SMARTS the chemist typed in. Toggling a catalog chip adds / removes
 * the scaffold's canonical name from `?colour_by=`; submitting the
 * SMARTS input adds an entry prefixed `smarts:`. The page-level
 * effect picks up the URL change and recomputes per-point colouring
 * via the RDKit context.
 *
 * Default visibility: catalog scaffolds with at least one matching
 * compound in the current scatter. A toggle exposes the full catalog
 * when the chemist wants to flick on a scaffold currently absent.
 * Ad-hoc SMARTS chips are always shown alongside catalog chips.
 */

import { useCallback, useMemo, useState } from 'react';
import { Box, Chip, CircularProgress, IconButton, Stack, TextField, Tooltip, Typography } from '@mui/material';
import { Palette, AddCircleOutline, Close } from '@mui/icons-material';
import useSWR from 'swr';
import { usePathname, useRouter, useSearchParams } from 'next/navigation';
import { SUBSTRUCTURES_KEY, listSubstructures } from '@/lib/compounds/substructures-api';
import { SCATTER_CATEGORY_COLOURS } from '@/lib/compounds/scatter-palette';
import type { ScaffoldDefinition } from '@/lib/compounds/scaffold-matching';

const SMARTS_PREFIX = 'smarts:';

function truncate(s: string, n: number): string {
  return s.length > n ? `${s.slice(0, n - 1)}…` : s;
}

interface Props {
  /**
   * Map from scaffold name to the set of compound formatted_ids that
   * carry it. Computed at the page level via RDKit. The strip uses
   * this only for filtering chips by relevance and showing the
   * matching count.
   */
  matchesByScaffold: ReadonlyMap<string, ReadonlySet<string>>;
}

export function ScatterCategorisationStrip({ matchesByScaffold }: Props) {
  const router = useRouter();
  const pathname = usePathname();
  const searchParams = useSearchParams();

  const colourByEntries = useMemo(() => {
    return searchParams?.get('colour_by')?.split(',').map((s) => s.trim()).filter(Boolean) ?? [];
  }, [searchParams]);
  const adhocSmarts = useMemo(
    () => colourByEntries.filter((e) => e.startsWith(SMARTS_PREFIX)),
    [colourByEntries],
  );
  const colourByNames = useMemo(
    () => colourByEntries.filter((e) => !e.startsWith(SMARTS_PREFIX)),
    [colourByEntries],
  );

  const [showAll, setShowAll] = useState(false);
  const [smartsInput, setSmartsInput] = useState<string | null>(null);

  const { data: scaffolds, isLoading } = useSWR<ScaffoldDefinition[]>(
    SUBSTRUCTURES_KEY,
    () => listSubstructures(),
    { revalidateOnFocus: false },
  );

  const setColourBy = useCallback(
    (entries: string[]) => {
      const params = new URLSearchParams(searchParams?.toString() ?? '');
      if (entries.length === 0) params.delete('colour_by');
      else params.set('colour_by', entries.join(','));
      router.replace(`${pathname}?${params.toString()}`);
    },
    [router, pathname, searchParams],
  );

  const toggle = useCallback(
    (entry: string) => {
      setColourBy(
        colourByEntries.includes(entry)
          ? colourByEntries.filter((x) => x !== entry)
          : [...colourByEntries, entry],
      );
    },
    [colourByEntries, setColourBy],
  );

  const submitSmarts = useCallback(
    (raw: string) => {
      const smarts = raw.trim();
      if (!smarts) return;
      // SMARTS may contain commas (atom lists like [c,n]). Reject those
      // here — the URL contract uses ',' as the entry separator. The
      // chemist can rephrase to bracket-list form without the comma if
      // needed; nothing in the curated catalog needs it.
      if (smarts.includes(',')) {
        return;
      }
      const entry = `${SMARTS_PREFIX}${smarts}`;
      if (colourByEntries.includes(entry)) return;
      setColourBy([...colourByEntries, entry]);
      setSmartsInput(null);
    },
    [colourByEntries, setColourBy],
  );

  // Visible catalog chips: scaffolds with ≥1 matching compound, plus
  // any currently-active scaffold (so toggling never silently drops
  // one). Ad-hoc SMARTS chips are tracked separately and always shown.
  const { visible, relevantCount, totalCount } = useMemo(() => {
    if (!scaffolds) return { visible: [] as ScaffoldDefinition[], relevantCount: 0, totalCount: 0 };
    const relevant = scaffolds.filter((s) => (matchesByScaffold.get(s.name)?.size ?? 0) > 0);
    const pool = showAll ? scaffolds : relevant;
    const poolNames = new Set(pool.map((s) => s.name));
    const orphans = scaffolds.filter((s) => colourByNames.includes(s.name) && !poolNames.has(s.name));
    return { visible: [...pool, ...orphans], relevantCount: relevant.length, totalCount: scaffolds.length };
  }, [scaffolds, showAll, colourByNames, matchesByScaffold]);

  if (isLoading) {
    return (
      <Box sx={{ mb: 2, p: 1.5, display: 'flex', alignItems: 'center', gap: 1, color: 'text.secondary' }}>
        <CircularProgress size={14} />
        <Typography variant="caption">Loading chemotypes…</Typography>
      </Box>
    );
  }
  if (!scaffolds) return null;
  if (visible.length === 0 && colourByEntries.length === 0 && smartsInput === null) return null;

  // Index lookup for swatch colours — colourByEntries is the source
  // of truth for ordering (matches the dataset-split priority in the
  // chart), so we read from it for both catalog and ad-hoc chips.
  const swatchFor = (entry: string): string | null => {
    const i = colourByEntries.indexOf(entry);
    return i >= 0 ? SCATTER_CATEGORY_COLOURS[i % SCATTER_CATEGORY_COLOURS.length].border : null;
  };

  return (
    <Box sx={{ mb: 2, p: 1.5, border: 1, borderColor: 'divider', borderRadius: 1 }}>
      <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, mb: 1, flexWrap: 'wrap' }}>
        <Palette fontSize="small" color="action" />
        <Typography variant="caption" color="text.secondary">
          Colour points by chemotype
        </Typography>
        {totalCount > relevantCount && (
          <Typography
            component="button"
            onClick={() => setShowAll((s) => !s)}
            sx={{
              border: 'none', bgcolor: 'transparent', cursor: 'pointer',
              ml: 1, color: 'primary.main', fontSize: '0.72rem',
              '&:hover': { textDecoration: 'underline' },
            }}
          >
            {showAll
              ? `Relevant only (${relevantCount})`
              : `Show full catalog (${totalCount})`}
          </Typography>
        )}
        {colourByEntries.length > 0 && (
          <Typography
            component="button"
            onClick={() => setColourBy([])}
            sx={{
              border: 'none', bgcolor: 'transparent', cursor: 'pointer',
              ml: 'auto', color: 'primary.main', fontSize: '0.72rem',
              '&:hover': { textDecoration: 'underline' },
            }}
          >
            Clear all
          </Typography>
        )}
        <Typography
          component="a"
          href="/registry/scaffold-extensions"
          sx={{
            ml: colourByEntries.length > 0 ? 1 : 'auto',
            color: 'primary.main',
            fontSize: '0.72rem',
            textDecoration: 'none',
            '&:hover': { textDecoration: 'underline' },
          }}
        >
          Manage scaffolds
        </Typography>
      </Box>

      <Stack direction="row" spacing={1} sx={{ flexWrap: 'wrap', gap: 0.75, alignItems: 'center' }}>
        {visible.map((sc) => {
          const enabled = colourByEntries.includes(sc.name);
          const swatch = enabled ? swatchFor(sc.name) : null;
          const matchCount = matchesByScaffold.get(sc.name)?.size ?? 0;
          return (
            <Tooltip
              key={sc.name}
              title={`${matchCount} matching compound${matchCount === 1 ? '' : 's'}${sc.source ? ` · ${sc.source}` : ''}`}
            >
              <Chip
                label={`${sc.name}${matchCount > 0 ? ` (${matchCount})` : ''}`}
                size="small"
                onClick={() => toggle(sc.name)}
                variant={enabled ? 'filled' : 'outlined'}
                sx={{
                  ...(swatch && {
                    bgcolor: swatch,
                    color: '#fff',
                    '&:hover': { bgcolor: swatch, opacity: 0.85 },
                  }),
                }}
              />
            </Tooltip>
          );
        })}

        {adhocSmarts.map((entry) => {
          const smarts = entry.slice(SMARTS_PREFIX.length);
          const matchCount = matchesByScaffold.get(entry)?.size ?? 0;
          const swatch = swatchFor(entry);
          return (
            <Tooltip key={entry} title={`SMARTS: ${smarts} · ${matchCount} matching`}>
              <Chip
                label={`SMARTS: ${truncate(smarts, 24)}${matchCount > 0 ? ` (${matchCount})` : ''}`}
                size="small"
                onDelete={() => setColourBy(colourByEntries.filter((x) => x !== entry))}
                deleteIcon={<Close fontSize="small" />}
                variant="filled"
                sx={{
                  ...(swatch && {
                    bgcolor: swatch,
                    color: '#fff',
                    '& .MuiChip-deleteIcon': { color: 'rgba(255,255,255,0.85)' },
                  }),
                }}
              />
            </Tooltip>
          );
        })}

        {smartsInput === null ? (
          <Tooltip title="Add an ad-hoc SMARTS pattern">
            <IconButton size="small" onClick={() => setSmartsInput('')}>
              <AddCircleOutline fontSize="small" />
            </IconButton>
          </Tooltip>
        ) : (
          <Box sx={{ display: 'flex', alignItems: 'center', gap: 0.5 }}>
            <TextField
              size="small"
              autoFocus
              placeholder="SMARTS pattern (e.g. c1ccncc1)"
              value={smartsInput}
              onChange={(e) => setSmartsInput(e.target.value)}
              onKeyDown={(e) => {
                if (e.key === 'Enter') submitSmarts(smartsInput);
                if (e.key === 'Escape') setSmartsInput(null);
              }}
              sx={{ minWidth: 220, '& .MuiInputBase-input': { fontFamily: 'monospace', fontSize: '0.8rem' } }}
              helperText={smartsInput.includes(',') ? 'Commas not supported in URL — use bracket lists without commas' : ''}
              error={smartsInput.includes(',')}
            />
            <IconButton size="small" onClick={() => submitSmarts(smartsInput)} disabled={!smartsInput.trim() || smartsInput.includes(',')}>
              <AddCircleOutline fontSize="small" color="primary" />
            </IconButton>
            <IconButton size="small" onClick={() => setSmartsInput(null)}>
              <Close fontSize="small" />
            </IconButton>
          </Box>
        )}
      </Stack>

      {visible.length === 0 && adhocSmarts.length === 0 && (
        <Typography variant="caption" color="text.secondary" sx={{ display: 'block', mt: 1 }}>
          No catalog scaffolds match any compound here. Add a SMARTS pattern with the + button.
        </Typography>
      )}
    </Box>
  );
}
