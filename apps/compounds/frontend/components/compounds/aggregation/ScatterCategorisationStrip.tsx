'use client';

/**
 * Chip strip mounted above the inline scatter view: each chip is one of
 * the chemist's saved selections, toggling it adds/removes that selection
 * from the `?colour_by=` URL param. The page-level effect picks up the
 * URL change and refetches the named selections to repaint the points.
 *
 * Defaults to showing only saved selections — the ephemeral pool is
 * dominated by recent ad-hoc queries that aren't useful as colour
 * categories. A toggle exposes the full pool when needed (e.g. the
 * chemist wants to colour by "the query I ran two hours ago" without
 * having saved it first).
 *
 * Currently-active categories that fall outside the visible set are
 * still rendered (so toggling "saved only" doesn't silently drop chips
 * the chemist is currently using).
 */

import { useCallback, useMemo, useState } from 'react';
import { Box, Chip, Stack, Tooltip, Typography } from '@mui/material';
import { Palette } from '@mui/icons-material';
import useSWR from 'swr';
import { usePathname, useRouter, useSearchParams } from 'next/navigation';
import {
  SELECTIONS_LIST_KEY,
  SelectionSummary,
  listSelections,
} from '@/lib/compounds/selections-api';
import { SCATTER_CATEGORY_COLOURS } from '@/lib/compounds/scatter-palette';

export function ScatterCategorisationStrip() {
  const router = useRouter();
  const pathname = usePathname();
  const searchParams = useSearchParams();

  const colourByIds = useMemo(() => {
    return searchParams?.get('colour_by')?.split(',').map((s) => s.trim()).filter(Boolean) ?? [];
  }, [searchParams]);

  const [showAll, setShowAll] = useState(false);

  const { data: selections } = useSWR<SelectionSummary[]>(
    SELECTIONS_LIST_KEY,
    () => listSelections(),
    { revalidateOnFocus: false },
  );

  const setColourBy = useCallback(
    (ids: string[]) => {
      const params = new URLSearchParams(searchParams?.toString() ?? '');
      if (ids.length === 0) params.delete('colour_by');
      else params.set('colour_by', ids.join(','));
      router.replace(`${pathname}?${params.toString()}`);
    },
    [router, pathname, searchParams],
  );

  const toggle = useCallback(
    (id: string) => {
      setColourBy(
        colourByIds.includes(id)
          ? colourByIds.filter((x) => x !== id)
          : [...colourByIds, id],
      );
    },
    [colourByIds, setColourBy],
  );

  // Two derived lists: the visible pool (saved-only or all) plus any
  // currently-coloured-by chip whose selection isn't in the pool —
  // hiding active chips would feel like the toggle "lost" them.
  const { visible, savedCount, totalCount } = useMemo(() => {
    if (!selections) return { visible: [] as SelectionSummary[], savedCount: 0, totalCount: 0 };
    const saved = selections.filter((s) => s.is_saved);
    const pool = showAll ? selections : saved;
    const poolIds = new Set(pool.map((s) => s.id));
    const orphans = selections.filter((s) => colourByIds.includes(s.id) && !poolIds.has(s.id));
    return { visible: [...pool, ...orphans], savedCount: saved.length, totalCount: selections.length };
  }, [selections, showAll, colourByIds]);

  if (!selections) return null;
  if (visible.length === 0 && colourByIds.length === 0) return null;

  return (
    <Box sx={{ mb: 2, p: 1.5, border: 1, borderColor: 'divider', borderRadius: 1 }}>
      <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, mb: 1, flexWrap: 'wrap' }}>
        <Palette fontSize="small" color="action" />
        <Typography variant="caption" color="text.secondary">
          Colour points by selection membership
        </Typography>
        {totalCount > savedCount && (
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
              ? `Saved only (${savedCount})`
              : `Show all my queries (${totalCount})`}
          </Typography>
        )}
        {colourByIds.length > 0 && (
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
      </Box>
      {visible.length === 0 ? (
        <Typography variant="caption" color="text.secondary" sx={{ display: 'block' }}>
          No saved selections yet — save one from the Ask page or the
          session list to colour by it here.
        </Typography>
      ) : (
        <Stack direction="row" spacing={1} sx={{ flexWrap: 'wrap', gap: 0.75 }}>
          {visible.map((sel) => {
            const enabledIndex = colourByIds.indexOf(sel.id);
            const enabled = enabledIndex >= 0;
            const swatch = enabled
              ? SCATTER_CATEGORY_COLOURS[enabledIndex % SCATTER_CATEGORY_COLOURS.length].border
              : null;
            return (
              <Tooltip
                key={sel.id}
                title={sel.is_saved ? `Saved · ${sel.n_compounds} compounds` : `${sel.n_compounds} compounds`}
              >
                <Chip
                  label={sel.name || `${sel.n_compounds} compounds`}
                  size="small"
                  onClick={() => toggle(sel.id)}
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
        </Stack>
      )}
    </Box>
  );
}
