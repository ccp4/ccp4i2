'use client';

import { useRef, useState, useLayoutEffect, useMemo } from 'react';
import {
  Box,
  Container,
  Paper,
  Typography,
  Skeleton,
  Collapse,
  Divider,
  IconButton,
  Tooltip,
  useTheme,
  useMediaQuery,
} from '@mui/material';
import { ExpandMore, ExpandLess } from '@mui/icons-material';
import { PageHeader } from './PageHeader';
import type {
  DetailPageLayoutProps,
  SummaryField,
} from '@/types/compounds/detail-page-layout';

const TRANSITION_DURATION = 250;
const BOTTOM_PADDING = 16;

/**
 * Renders a single summary field as "Label: Value"
 */
function SummaryFieldItem({ field }: { field: SummaryField }) {
  return (
    <Box
      sx={{
        display: 'flex',
        alignItems: 'center',
        gap: 0.5,
        minWidth: 0,
        overflow: 'hidden',
      }}
    >
      {field.icon}
      <Typography
        variant="body2"
        color="text.secondary"
        sx={{ fontWeight: 500, flexShrink: 0 }}
      >
        {field.label}:
      </Typography>
      <Typography
        variant="body2"
        sx={{
          overflow: 'hidden',
          textOverflow: 'ellipsis',
          whiteSpace: 'nowrap',
        }}
      >
        {field.value ?? '-'}
      </Typography>
    </Box>
  );
}

/**
 * Layout component for detail pages with a collapsible header.
 *
 * Features:
 * - Header collapses when scrolling down in the table, expands when scrolling up
 * - Collapsed state: compact summary with title + key fields
 * - Expanded state: full detail content (structure viewer, accordions, etc.)
 * - Table/content area fills remaining viewport height for virtualization
 * - Single scroll context (table scrolls, page does not)
 */
export function DetailPageLayout({
  breadcrumbs,
  summary,
  detailContent,
  loading = false,
  children,
  collapsedFieldCount = 3,
  defaultCollapsed = true,
}: DetailPageLayoutProps) {
  const theme = useTheme();
  const isMobile = useMediaQuery(theme.breakpoints.down('sm'));

  const headerRef = useRef<HTMLDivElement>(null);
  const [headerHeight, setHeaderHeight] = useState(0);
  const [isCollapsed, setIsCollapsed] = useState(defaultCollapsed);

  // Measure header height using ResizeObserver
  useLayoutEffect(() => {
    const headerEl = headerRef.current;
    if (!headerEl) return;

    const updateHeight = () => {
      const rect = headerEl.getBoundingClientRect();
      setHeaderHeight(rect.height);
    };

    // Initial measurement
    updateHeight();

    const resizeObserver = new ResizeObserver(updateHeight);
    resizeObserver.observe(headerEl);

    return () => {
      resizeObserver.disconnect();
    };
  }, []);

  // Calculate content area height
  const contentHeight = useMemo(() => {
    if (headerHeight === 0) return '50vh'; // Fallback while measuring
    return `calc(100vh - ${headerHeight}px - ${BOTTOM_PADDING}px)`;
  }, [headerHeight]);

  // Fields to show in collapsed inline summary
  const collapsedFields = useMemo(() => {
    const count = isMobile ? 1 : collapsedFieldCount;
    return summary.fields.slice(0, count);
  }, [summary.fields, collapsedFieldCount, isMobile]);

  return (
    <Box
      sx={{
        display: 'flex',
        flexDirection: 'column',
        height: '100vh',
        overflow: 'hidden',
        bgcolor: 'background.default',
      }}
    >
      {/* Header - not sticky, just fixed at top */}
      <Paper
        ref={headerRef}
        elevation={isCollapsed ? 2 : 0}
        sx={{
          flexShrink: 0,
          bgcolor: 'background.paper',
          borderRadius: 0,
          transition: `box-shadow ${TRANSITION_DURATION}ms ease-in-out`,
        }}
      >
        <Container maxWidth="xl">
          {/* Breadcrumbs row */}
          <Box sx={{ pt: 1, pb: isCollapsed ? 0 : 1 }}>
            <PageHeader breadcrumbs={breadcrumbs} hideActions={isCollapsed && !isMobile} />
          </Box>

          {/* Title and summary row */}
          <Box
            sx={{
              display: 'flex',
              alignItems: 'center',
              gap: 2,
              py: isCollapsed ? 1 : 1.5,
              transition: `padding ${TRANSITION_DURATION}ms ease-in-out`,
            }}
          >
            {/* Title icon */}
            {loading ? (
              <Skeleton
                variant="circular"
                width={isCollapsed ? 32 : 40}
                height={isCollapsed ? 32 : 40}
              />
            ) : (
              summary.titleIcon && (
                <Box
                  sx={{
                    fontSize: isCollapsed ? 32 : 40,
                    color: 'primary.main',
                    display: 'flex',
                    alignItems: 'center',
                    transition: `font-size ${TRANSITION_DURATION}ms ease-in-out`,
                  }}
                >
                  {summary.titleIcon}
                </Box>
              )
            )}

            {/* Title and inline summary */}
            <Box sx={{ flex: 1, minWidth: 0 }}>
              {loading ? (
                <Skeleton variant="text" width={200} height={32} />
              ) : (
                <Typography
                  variant={isCollapsed ? 'h6' : 'h5'}
                  sx={{
                    overflow: 'hidden',
                    textOverflow: 'ellipsis',
                    whiteSpace: 'nowrap',
                    transition: `all ${TRANSITION_DURATION}ms ease-in-out`,
                    lineHeight: 1.3,
                  }}
                >
                  {summary.title}
                </Typography>
              )}

              {/* Inline summary fields - shown when collapsed */}
              {isCollapsed && !loading && collapsedFields.length > 0 && (
                <Box
                  sx={{
                    display: 'flex',
                    gap: 2,
                    flexWrap: 'nowrap',
                    overflow: 'hidden',
                    mt: 0.25,
                  }}
                >
                  {collapsedFields.map((field, i) => (
                    <SummaryFieldItem key={i} field={field} />
                  ))}
                </Box>
              )}
            </Box>

            {/* Chips */}
            {!loading && summary.chips && (
              <Box sx={{ display: 'flex', gap: 0.5, flexShrink: 0 }}>
                {summary.chips}
              </Box>
            )}

            {/* Action buttons - always visible */}
            {!loading && summary.actions && (
              <Box sx={{ display: 'flex', gap: 1, flexShrink: 0 }}>
                {summary.actions}
              </Box>
            )}

            {/* Expand/Collapse toggle */}
            {!loading && detailContent && (
              <Tooltip title={isCollapsed ? 'Show details' : 'Hide details'}>
                <IconButton
                  onClick={() => setIsCollapsed(!isCollapsed)}
                  size="small"
                  sx={{ ml: 0.5 }}
                >
                  {isCollapsed ? <ExpandMore /> : <ExpandLess />}
                </IconButton>
              </Tooltip>
            )}
          </Box>

          {/* Expandable detail content */}
          <Collapse in={!isCollapsed} timeout={TRANSITION_DURATION}>
            {!loading && (
              <Box sx={{ pb: 2 }}>
                <Divider sx={{ mb: 2 }} />
                {detailContent}
              </Box>
            )}
          </Collapse>
        </Container>
      </Paper>

      {/* Content area (table) - fills remaining viewport, table handles its own scroll */}
      <Container
        maxWidth="xl"
        sx={{
          flex: 1,
          display: 'flex',
          flexDirection: 'column',
          py: 1,
          minHeight: 0,
          overflow: 'hidden',
        }}
      >
        <Box
          sx={{
            flex: 1,
            minHeight: 0,
            display: 'flex',
            flexDirection: 'column',
          }}
        >
          {children}
        </Box>
      </Container>
    </Box>
  );
}
