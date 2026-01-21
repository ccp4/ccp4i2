'use client';

import { useRef, ReactNode, useCallback } from 'react';
import {
  Box,
  Paper,
  Typography,
  IconButton,
  Skeleton,
} from '@mui/material';
import { ChevronLeft, ChevronRight } from '@mui/icons-material';
import { useVirtualizer } from '@tanstack/react-virtual';

interface HorizontalCarouselProps<T> {
  /** Array of items to display */
  items: T[];
  /** Render function for each item */
  renderItem: (item: T, index: number) => ReactNode;
  /** Function to get unique key for each item */
  getItemKey: (item: T) => string | number;
  /** Optional title displayed above the carousel */
  title?: string;
  /** Width of each item in pixels (default: 200) */
  itemWidth?: number;
  /** Height of the carousel container in pixels (default: 280) */
  height?: number;
  /** Gap between items in pixels (default: 16) */
  gap?: number;
  /** Callback when clicking the carousel background (not on an item) */
  onBackgroundClick?: () => void;
  /** Whether data is loading */
  loading?: boolean;
  /** Message to show when there are no items */
  emptyMessage?: string;
  /** Optional action element to display in the header (e.g., an "Add" button) */
  headerAction?: ReactNode;
}

export function HorizontalCarousel<T>({
  items,
  renderItem,
  getItemKey,
  title,
  itemWidth = 200,
  height = 280,
  gap = 16,
  onBackgroundClick,
  loading,
  emptyMessage = 'No items',
  headerAction,
}: HorizontalCarouselProps<T>) {
  const scrollContainerRef = useRef<HTMLDivElement>(null);

  const virtualizer = useVirtualizer({
    count: items.length,
    getScrollElement: () => scrollContainerRef.current,
    estimateSize: () => itemWidth + gap,
    horizontal: true,
    overscan: 3,
  });

  const scrollLeft = useCallback(() => {
    if (scrollContainerRef.current) {
      scrollContainerRef.current.scrollBy({
        left: -(itemWidth + gap) * 3,
        behavior: 'smooth',
      });
    }
  }, [itemWidth, gap]);

  const scrollRight = useCallback(() => {
    if (scrollContainerRef.current) {
      scrollContainerRef.current.scrollBy({
        left: (itemWidth + gap) * 3,
        behavior: 'smooth',
      });
    }
  }, [itemWidth, gap]);

  const handleBackgroundClick = useCallback(
    (e: React.MouseEvent<HTMLDivElement>) => {
      // Only trigger if clicking directly on the scroll container or its inner wrapper
      if (
        e.target === e.currentTarget ||
        (e.target as HTMLElement).dataset.carouselBackground === 'true'
      ) {
        onBackgroundClick?.();
      }
    },
    [onBackgroundClick]
  );

  if (loading) {
    return (
      <Paper sx={{ p: 2, mb: 3 }}>
        {title && (
          <Typography variant="h6" sx={{ mb: 2 }}>
            {title}
          </Typography>
        )}
        <Box sx={{ display: 'flex', gap: `${gap}px`, overflow: 'hidden' }}>
          {[1, 2, 3, 4, 5].map((i) => (
            <Skeleton
              key={i}
              variant="rounded"
              width={itemWidth}
              height={height - 40}
            />
          ))}
        </Box>
      </Paper>
    );
  }

  return (
    <Paper sx={{ p: 2, mb: 3 }}>
      {/* Header with title, action, and scroll buttons */}
      <Box
        sx={{
          display: 'flex',
          alignItems: 'center',
          gap: 1,
          mb: 2,
        }}
      >
        {title && (
          <Typography variant="h6" sx={{ flex: 1 }}>
            {title}
          </Typography>
        )}
        {headerAction}
        {items.length > 0 && (
          <Box>
            <IconButton
              onClick={scrollLeft}
              size="small"
              aria-label="Scroll left"
            >
              <ChevronLeft />
            </IconButton>
            <IconButton
              onClick={scrollRight}
              size="small"
              aria-label="Scroll right"
            >
              <ChevronRight />
            </IconButton>
          </Box>
        )}
      </Box>

      {/* Carousel content */}
      {items.length === 0 ? (
        <Box
          sx={{
            height: height - 40,
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            bgcolor: 'grey.50',
            borderRadius: 1,
            cursor: onBackgroundClick ? 'pointer' : 'default',
          }}
          onClick={onBackgroundClick}
        >
          <Typography color="text.secondary">{emptyMessage}</Typography>
        </Box>
      ) : (
        <Box
          ref={scrollContainerRef}
          onClick={handleBackgroundClick}
          sx={{
            height,
            overflow: 'auto',
            scrollbarWidth: 'thin',
            cursor: onBackgroundClick ? 'pointer' : 'default',
            '&::-webkit-scrollbar': {
              height: 8,
            },
            '&::-webkit-scrollbar-track': {
              bgcolor: 'grey.100',
              borderRadius: 4,
            },
            '&::-webkit-scrollbar-thumb': {
              bgcolor: 'grey.400',
              borderRadius: 4,
              '&:hover': {
                bgcolor: 'grey.500',
              },
            },
          }}
        >
          <Box
            data-carousel-background="true"
            sx={{
              height: '100%',
              width: virtualizer.getTotalSize(),
              position: 'relative',
              minWidth: '100%',
            }}
          >
            {virtualizer.getVirtualItems().map((virtualItem) => {
              const item = items[virtualItem.index];
              return (
                <Box
                  key={getItemKey(item)}
                  sx={{
                    position: 'absolute',
                    top: 0,
                    left: virtualItem.start,
                    width: itemWidth,
                    height: '100%',
                    pr: `${gap}px`,
                    boxSizing: 'border-box',
                  }}
                >
                  <Box
                    sx={{
                      height: '100%',
                      cursor: 'pointer',
                    }}
                    onClick={(e) => e.stopPropagation()}
                  >
                    {renderItem(item, virtualItem.index)}
                  </Box>
                </Box>
              );
            })}
          </Box>
        </Box>
      )}
    </Paper>
  );
}
