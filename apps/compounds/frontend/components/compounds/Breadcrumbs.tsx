'use client';

import { Breadcrumbs as MuiBreadcrumbs, Link, Typography, Box } from '@mui/material';
import { Home, NavigateNext, Science, Medication, Inventory, Description, Assessment, TableChart, Search, Upload } from '@mui/icons-material';
import NextLink from 'next/link';

export interface BreadcrumbItem {
  label: string;
  href?: string;
  icon?: 'home' | 'target' | 'compound' | 'batch' | 'qc' | 'protocol' | 'assay' | 'aggregate' | 'search' | 'import' | 'construct';
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
};

export function Breadcrumbs({ items }: BreadcrumbsProps) {
  return (
    <MuiBreadcrumbs
      separator={<NavigateNext fontSize="small" />}
      sx={{ mb: 2 }}
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
