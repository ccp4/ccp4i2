'use client';

import { Breadcrumbs as MuiBreadcrumbs, Link, Typography, Box, useMediaQuery, useTheme } from '@mui/material';
import { Home, NavigateNext, Science, Medication, Inventory, Description, Assessment, TableChart, Search, Upload, LocalShipping, ChevronLeft } from '@mui/icons-material';
import NextLink from 'next/link';

export interface BreadcrumbItem {
  label: string;
  href?: string;
  icon?: 'home' | 'target' | 'compound' | 'batch' | 'qc' | 'protocol' | 'assay' | 'aggregate' | 'search' | 'import' | 'construct' | 'supplier';
}

interface BreadcrumbsProps {
  items: BreadcrumbItem[];
}

const iconMap = {
  home: <Home fontSize="small" />,
  target: <Science fontSize="small" />,
  compound: <Medication fontSize="small" />,
  batch: <Inventory fontSize="small" />,
  qc: <Description fontSize="small" />,
  protocol: <Description fontSize="small" />,
  assay: <Assessment fontSize="small" />,
  aggregate: <TableChart fontSize="small" />,
  search: <Search fontSize="small" />,
  import: <Upload fontSize="small" />,
  construct: <Science fontSize="small" />,
  supplier: <LocalShipping fontSize="small" />,
};

/**
 * Responsive breadcrumbs component.
 * - Desktop: Full breadcrumb trail with icons
 * - Mobile: Compact back link to parent ("< Parent Name")
 */
export function Breadcrumbs({ items }: BreadcrumbsProps) {
  const theme = useTheme();
  const isMobile = useMediaQuery(theme.breakpoints.down('sm'));

  // Mobile: Show back link to parent
  if (isMobile && items.length > 1) {
    // Find the last item with an href (the parent we can navigate back to)
    const parentIndex = items.slice(0, -1).findLastIndex(item => item.href);
    const parent = parentIndex >= 0 ? items[parentIndex] : null;

    if (parent?.href) {
      return (
        <Link
          component={NextLink}
          href={parent.href}
          underline="hover"
          sx={{
            display: 'inline-flex',
            alignItems: 'center',
            gap: 0.5,
            color: 'text.secondary',
            fontSize: '0.875rem',
            py: 0.5,
            '&:hover': { color: 'primary.main' },
          }}
        >
          <ChevronLeft fontSize="small" sx={{ ml: -0.5 }} />
          {parent.label}
        </Link>
      );
    }

    // No navigable parent - show current page name only
    return (
      <Typography
        variant="body2"
        color="text.secondary"
        sx={{ py: 0.5 }}
      >
        {items[items.length - 1].label}
      </Typography>
    );
  }

  // Desktop: Full breadcrumb trail
  return (
    <MuiBreadcrumbs
      separator={<NavigateNext fontSize="small" />}
      sx={{ py: 0.5 }}
    >
      {items.map((item, index) => {
        const isLast = index === items.length - 1;
        const icon = item.icon ? iconMap[item.icon] : null;

        if (isLast || !item.href) {
          return (
            <Box
              key={index}
              sx={{
                display: 'flex',
                alignItems: 'center',
                gap: 0.5,
                color: isLast ? 'text.primary' : 'text.secondary',
              }}
            >
              {icon}
              <Typography
                color={isLast ? 'text.primary' : 'text.secondary'}
                sx={{ fontWeight: isLast ? 600 : 400 }}
              >
                {item.label}
              </Typography>
            </Box>
          );
        }

        return (
          <Link
            key={index}
            component={NextLink}
            href={item.href}
            underline="hover"
            sx={{
              display: 'flex',
              alignItems: 'center',
              gap: 0.5,
              color: 'text.secondary',
              '&:hover': { color: 'primary.main' },
            }}
          >
            {icon}
            {item.label}
          </Link>
        );
      })}
    </MuiBreadcrumbs>
  );
}
