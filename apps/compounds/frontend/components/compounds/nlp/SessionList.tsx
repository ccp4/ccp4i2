'use client';

/**
 * Session-list panel — chips for the chemist's saved + recent NLP queries.
 *
 * Each row is a `Selection` row created server-side by the NLP view on
 * every successful compound query (or via direct POST). Clicking a chip
 * navigates to /assays/aggregate?selection=<uuid>; the bookmark icon
 * toggles `is_saved` (saved selections never expire); the bin deletes
 * (with a confirm dialog).
 */

import { useCallback, useState } from 'react';
import {
  Box,
  Chip,
  CircularProgress,
  IconButton,
  Paper,
  Stack,
  Tooltip,
  Typography,
} from '@mui/material';
import {
  BookmarkAdded,
  BookmarkBorder,
  DeleteOutline,
  OpenInNew,
} from '@mui/icons-material';
import useSWR from 'swr';
import { useRouter } from 'next/navigation';
import {
  SELECTIONS_LIST_KEY,
  SelectionSummary,
  deleteSelection,
  listSelections,
  patchSelection,
} from '@/lib/compounds/selections-api';
import { ConfirmDialog } from '@/components/compounds/ConfirmDialog';

const formatRelative = (iso: string, now: number): string => {
  const ms = now - new Date(iso).getTime();
  const mins = Math.round(ms / 60_000);
  if (mins < 1) return 'just now';
  if (mins < 60) return `${mins}m ago`;
  const hours = Math.round(mins / 60);
  if (hours < 24) return `${hours}h ago`;
  const days = Math.round(hours / 24);
  if (days < 7) return `${days}d ago`;
  return new Date(iso).toLocaleDateString();
};

export function SessionList() {
  const router = useRouter();
  const { data, error, isLoading, mutate } = useSWR<SelectionSummary[]>(
    SELECTIONS_LIST_KEY,
    () => listSelections(),
    { revalidateOnFocus: false },
  );

  const [confirmDelete, setConfirmDelete] = useState<SelectionSummary | null>(null);

  const handleOpen = useCallback((sel: SelectionSummary) => {
    router.push(`/assays/aggregate?selection=${sel.id}&format=cards`);
  }, [router]);

  const handleToggleSave = useCallback(
    async (sel: SelectionSummary) => {
      await patchSelection(sel.id, { is_saved: !sel.is_saved });
      mutate();
    },
    [mutate],
  );

  const handleConfirmDelete = useCallback(async () => {
    if (!confirmDelete) return;
    await deleteSelection(confirmDelete.id);
    setConfirmDelete(null);
    mutate();
  }, [confirmDelete, mutate]);

  if (isLoading) {
    return (
      <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, color: 'text.secondary' }}>
        <CircularProgress size={14} /> <Typography variant="caption">Loading session…</Typography>
      </Box>
    );
  }
  if (error) return null;
  const selections = data ?? [];
  if (selections.length === 0) return null;

  // Compute "now" once per render so the relative-time labels for all
  // chips agree and we don't pay the cost of `new Date()` per row.
  const now = Date.now();

  return (
    <>
      <Paper variant="outlined" sx={{ p: 2 }}>
        <Typography variant="overline" color="text.secondary" sx={{ display: 'block', mb: 1 }}>
          Your recent and saved queries
        </Typography>
        <Stack direction="row" spacing={1} sx={{ flexWrap: 'wrap', gap: 1 }}>
          {selections.map((sel) => (
            <Chip
              key={sel.id}
              icon={sel.is_saved ? <BookmarkAdded fontSize="small" /> : undefined}
              label={
                <Box sx={{ display: 'flex', alignItems: 'center', gap: 0.75 }}>
                  <Box component="span" sx={{ maxWidth: 360, overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}>
                    {sel.name || `${sel.n_compounds} compounds`}
                  </Box>
                  <Box component="span" sx={{ color: 'text.secondary', fontSize: '0.72rem' }}>
                    {sel.n_compounds}
                    {' · '}
                    {formatRelative(sel.created_at, now)}
                  </Box>
                  <Tooltip title={sel.is_saved ? 'Saved — click to unsave' : 'Save (never expire)'}>
                    <IconButton
                      size="small"
                      onClick={(e) => {
                        e.stopPropagation();
                        handleToggleSave(sel);
                      }}
                      sx={{ ml: 0.25, p: 0.25 }}
                    >
                      {sel.is_saved
                        ? <BookmarkAdded fontSize="inherit" />
                        : <BookmarkBorder fontSize="inherit" />}
                    </IconButton>
                  </Tooltip>
                  <Tooltip title="Open in aggregation view">
                    <IconButton
                      size="small"
                      onClick={(e) => {
                        e.stopPropagation();
                        handleOpen(sel);
                      }}
                      sx={{ p: 0.25 }}
                    >
                      <OpenInNew fontSize="inherit" />
                    </IconButton>
                  </Tooltip>
                  <Tooltip title="Delete">
                    <IconButton
                      size="small"
                      onClick={(e) => {
                        e.stopPropagation();
                        setConfirmDelete(sel);
                      }}
                      sx={{ p: 0.25 }}
                    >
                      <DeleteOutline fontSize="inherit" />
                    </IconButton>
                  </Tooltip>
                </Box>
              }
              variant={sel.is_saved ? 'filled' : 'outlined'}
              onClick={() => handleOpen(sel)}
              sx={{ height: 'auto', '& .MuiChip-label': { display: 'flex', py: 0.5 } }}
            />
          ))}
        </Stack>
      </Paper>

      <ConfirmDialog
        open={confirmDelete !== null}
        title="Delete selection?"
        message={
          <>
            Remove {confirmDelete?.name ? <strong>{confirmDelete.name}</strong> : 'this selection'} from your session list?
            This can&apos;t be undone.
          </>
        }
        onConfirm={handleConfirmDelete}
        onCancel={() => setConfirmDelete(null)}
      />
    </>
  );
}
