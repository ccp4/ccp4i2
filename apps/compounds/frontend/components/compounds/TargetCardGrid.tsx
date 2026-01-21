'use client';

import { useRef, useMemo, useState, useCallback, useEffect } from 'react';
import { Box, TextField, Typography, InputAdornment, Paper, LinearProgress } from '@mui/material';
import { Search as SearchIcon } from '@mui/icons-material';
import { useVirtualizer } from '@tanstack/react-virtual';
import { TargetCard } from './TargetCard';
import { Target } from '@/types/compounds/models';

interface TargetCardGridProps {
  targets: Target[] | undefined;
  loading?: boolean;
  onTargetClick?: (target: Target) => void;
  emptyMessage?: string;
}

// Card dimensions
const CARD_WIDTH = 280;
const CARD_HEIGHT = 220;
const CARD_GAP = 16;

export function TargetCardGrid({
  targets,
  loading,
  onTargetClick,
  emptyMessage = 'No targets found',
}: TargetCardGridProps) {
  const parentRef = useRef<HTMLDivElement>(null);
  const [searchQuery, setSearchQuery] = useState('');
  const [columns, setColumns] = useState(3);

  // Calculate number of columns based on container width
  const updateColumns = useCallback(() => {
    if (parentRef.current) {
      const containerWidth = parentRef.current.clientWidth;
      const newColumns = Math.max(1, Math.floor((containerWidth + CARD_GAP) / (CARD_WIDTH + CARD_GAP)));
      setColumns(newColumns);
    }
  }, []);

  // Update columns on mount and resize
  useEffect(() => {
    updateColumns();
    const resizeObserver = new ResizeObserver(updateColumns);
    if (parentRef.current) {
      resizeObserver.observe(parentRef.current);
    }
    return () => resizeObserver.disconnect();
  }, [updateColumns]);

  // Filter targets by search query
  const filteredTargets = useMemo(() => {
    if (!targets) return [];
    if (!searchQuery) return targets;

    const query = searchQuery.toLowerCase();
    return targets.filter((target) =>
      target.name.toLowerCase().includes(query)
    );
  }, [targets, searchQuery]);

  // Calculate number of rows
  const rowCount = Math.ceil(filteredTargets.length / columns);

  // Set up virtualizer for rows
  const rowVirtualizer = useVirtualizer({
    count: rowCount,
    getScrollElement: () => parentRef.current,
    estimateSize: () => CARD_HEIGHT + CARD_GAP,
    overscan: 2,
  });

  const virtualRows = rowVirtualizer.getVirtualItems();

  return (
    <Paper sx={{ width: '100%', overflow: 'hidden' }}>
      {/* Header with search */}
      <Box sx={{ p: 2, display: 'flex', alignItems: 'center', gap: 2 }}>
        <Typography variant="h6" sx={{ flexGrow: 1 }}>
          {filteredTargets.length} {filteredTargets.length === 1 ? 'target' : 'targets'}
        </Typography>
        <TextField
          size="small"
          placeholder="Search targets..."
          value={searchQuery}
          onChange={(e) => setSearchQuery(e.target.value)}
          InputProps={{
            startAdornment: (
              <InputAdornment position="start">
                <SearchIcon fontSize="small" />
              </InputAdornment>
            ),
          }}
          sx={{ minWidth: 200 }}
        />
      </Box>

      {/* Loading indicator */}
      {loading && <LinearProgress />}

      {/* Grid container */}
      <Box
        ref={parentRef}
        sx={{
          height: 600,
          overflow: 'auto',
          px: 2,
          pb: 2,
        }}
      >
        {filteredTargets.length === 0 ? (
          <Box sx={{ py: 8, textAlign: 'center' }}>
            <Typography color="text.secondary">
              {loading ? 'Loading...' : emptyMessage}
            </Typography>
          </Box>
        ) : (
          <Box
            sx={{
              height: rowVirtualizer.getTotalSize(),
              width: '100%',
              position: 'relative',
            }}
          >
            {virtualRows.map((virtualRow) => {
              const rowStartIndex = virtualRow.index * columns;
              const rowTargets = filteredTargets.slice(rowStartIndex, rowStartIndex + columns);

              return (
                <Box
                  key={virtualRow.key}
                  sx={{
                    position: 'absolute',
                    top: virtualRow.start,
                    left: 0,
                    width: '100%',
                    height: CARD_HEIGHT,
                    display: 'grid',
                    gridTemplateColumns: `repeat(${columns}, minmax(0, 1fr))`,
                    gap: `${CARD_GAP}px`,
                  }}
                >
                  {rowTargets.map((target) => (
                    <TargetCard
                      key={target.id}
                      target={target}
                      onClick={() => onTargetClick?.(target)}
                    />
                  ))}
                </Box>
              );
            })}
          </Box>
        )}
      </Box>
    </Paper>
  );
}
