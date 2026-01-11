"use client";

import { useRef, useMemo, ReactNode, useCallback } from "react";
import { Box, Typography } from "@mui/material";
import { useVirtualizer } from "@tanstack/react-virtual";

interface VirtualizedCardGridProps<T> {
  /** Data items to render */
  data: T[];
  /** Function to render each card */
  renderCard: (item: T, index: number) => ReactNode;
  /** Function to get unique key for each item */
  getItemKey: (item: T) => string | number;
  /** Number of columns at different breakpoints */
  columns?: {
    xs?: number;
    sm?: number;
    md?: number;
    lg?: number;
    xl?: number;
  };
  /** Estimated card height for virtualization */
  estimateCardHeight?: number;
  /** Gap between cards in pixels */
  gap?: number;
  /** Maximum height of the container */
  maxHeight?: number;
  /** Message when no data */
  emptyMessage?: string;
  /** Loading state */
  loading?: boolean;
}

/**
 * A virtualized grid component for rendering large lists of cards efficiently.
 * Only renders cards that are visible in the viewport plus a small overscan.
 */
export function VirtualizedCardGrid<T>({
  data,
  renderCard,
  getItemKey,
  columns = { xs: 1, sm: 2, md: 3, lg: 4, xl: 4 },
  estimateCardHeight = 280,
  gap = 24,
  maxHeight = 800,
  emptyMessage = "No items found",
  loading = false,
}: VirtualizedCardGridProps<T>) {
  const parentRef = useRef<HTMLDivElement>(null);

  // Calculate current number of columns based on container width
  // We'll use a simple approach - detect on mount and resize
  const getColumnsCount = useCallback(() => {
    if (typeof window === "undefined") return columns.md || 3;
    const width = window.innerWidth;
    if (width < 600) return columns.xs || 1;
    if (width < 900) return columns.sm || 2;
    if (width < 1200) return columns.md || 3;
    if (width < 1536) return columns.lg || 4;
    return columns.xl || 4;
  }, [columns]);

  // For simplicity, we'll use a fixed column count based on common breakpoints
  // In a more sophisticated implementation, we'd use ResizeObserver
  const columnsCount = getColumnsCount();

  // Group items into rows
  const rows = useMemo(() => {
    const result: T[][] = [];
    for (let i = 0; i < data.length; i += columnsCount) {
      result.push(data.slice(i, i + columnsCount));
    }
    return result;
  }, [data, columnsCount]);

  // Set up virtualizer for rows
  const rowVirtualizer = useVirtualizer({
    count: rows.length,
    getScrollElement: () => parentRef.current,
    estimateSize: () => estimateCardHeight + gap,
    overscan: 2, // Render 2 extra rows above/below for smoother scrolling
  });

  const virtualItems = rowVirtualizer.getVirtualItems();

  if (data.length === 0) {
    return (
      <Box sx={{ textAlign: "center", py: 8 }}>
        <Typography variant="h6" color="text.secondary" gutterBottom>
          {loading ? "Loading..." : emptyMessage}
        </Typography>
      </Box>
    );
  }

  return (
    <Box
      ref={parentRef}
      sx={{
        maxHeight,
        overflow: "auto",
        position: "relative",
      }}
    >
      {/* Container with total height for proper scrollbar */}
      <Box
        sx={{
          height: rowVirtualizer.getTotalSize(),
          width: "100%",
          position: "relative",
        }}
      >
        {/* Render only visible rows */}
        {virtualItems.map((virtualRow) => {
          const rowItems = rows[virtualRow.index];
          return (
            <Box
              key={virtualRow.key}
              data-index={virtualRow.index}
              ref={rowVirtualizer.measureElement}
              sx={{
                position: "absolute",
                top: 0,
                left: 0,
                width: "100%",
                transform: `translateY(${virtualRow.start}px)`,
                display: "grid",
                gridTemplateColumns: `repeat(${columnsCount}, 1fr)`,
                gap: `${gap}px`,
                pb: `${gap}px`,
              }}
            >
              {rowItems.map((item, colIndex) => (
                <Box key={getItemKey(item)}>
                  {renderCard(item, virtualRow.index * columnsCount + colIndex)}
                </Box>
              ))}
            </Box>
          );
        })}
      </Box>
    </Box>
  );
}
