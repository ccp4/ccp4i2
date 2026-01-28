'use client';

import { useState, useMemo, useRef, ReactNode } from 'react';
import {
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  TableSortLabel,
  Paper,
  TextField,
  Box,
  Typography,
  LinearProgress,
  InputAdornment,
  useMediaQuery,
  useTheme,
} from '@mui/material';
import { Search as SearchIcon } from '@mui/icons-material';
import { useVirtualizer } from '@tanstack/react-virtual';

export interface Column<T> {
  key: string;
  label: string;
  sortable?: boolean;
  searchable?: boolean;
  render?: (value: any, row: T) => ReactNode;
  width?: string | number;
  /** Hide this column on mobile devices (screen width < 600px) */
  hiddenOnMobile?: boolean;
}

interface DataTableProps<T> {
  data: T[] | undefined;
  columns: Column<T>[];
  loading?: boolean;
  onRowClick?: (row: T) => void;
  getRowKey: (row: T) => string | number;
  title?: string;
  emptyMessage?: string;
  /** Additional field names to include in search (fields not displayed as columns) */
  additionalSearchFields?: string[];
  /** Estimated row height for virtualization (auto-calculated based on comfortable setting if not specified) */
  estimateRowHeight?: number;
  /** Maximum height of the table container (default: 600) */
  maxHeight?: number;
  /** Optional action element to display in the header (e.g., an "Add" button) */
  headerAction?: ReactNode;
  /** Use comfortable padding for better mobile/touch experience (default: false) */
  comfortable?: boolean;
}

type Order = 'asc' | 'desc';

// Row height constants
const DENSE_ROW_HEIGHT = 53;
const COMFORTABLE_ROW_HEIGHT = 64;

export function DataTable<T extends Record<string, any>>({
  data,
  columns,
  loading,
  onRowClick,
  getRowKey,
  title,
  emptyMessage = 'No data found',
  additionalSearchFields = [],
  estimateRowHeight,
  maxHeight = 600,
  headerAction,
  comfortable = false,
}: DataTableProps<T>) {
  const theme = useTheme();
  const isMobile = useMediaQuery(theme.breakpoints.down('sm'));

  // Filter columns based on mobile visibility
  const visibleColumns = useMemo(() => {
    if (!isMobile) return columns;
    return columns.filter(col => !col.hiddenOnMobile);
  }, [columns, isMobile]);

  // Calculate row height based on comfortable setting if not explicitly provided
  const effectiveRowHeight = estimateRowHeight ?? (comfortable ? COMFORTABLE_ROW_HEIGHT : DENSE_ROW_HEIGHT);
  const cellPadding = comfortable ? 2 : undefined; // MUI spacing units (16px per unit)
  const [orderBy, setOrderBy] = useState<string | null>(null);
  const [order, setOrder] = useState<Order>('asc');
  const [searchQuery, setSearchQuery] = useState('');

  const parentRef = useRef<HTMLDivElement>(null);

  // Filter data by search query
  const filteredData = useMemo(() => {
    if (!data) return [];
    if (!searchQuery) return data;

    const query = searchQuery.toLowerCase();
    return data.filter((row) => {
      // Search in columns marked as searchable
      const matchesColumn = columns.some((col) => {
        if (!col.searchable) return false;
        const value = row[col.key];
        if (value == null) return false;
        return String(value).toLowerCase().includes(query);
      });

      if (matchesColumn) return true;

      // Search in additional fields not displayed as columns
      return additionalSearchFields.some((field) => {
        const value = row[field];
        if (value == null) return false;
        return String(value).toLowerCase().includes(query);
      });
    });
  }, [data, searchQuery, columns, additionalSearchFields]);

  // Sort data
  const sortedData = useMemo(() => {
    if (!orderBy) return filteredData;

    return [...filteredData].sort((a, b) => {
      const aVal = a[orderBy];
      const bVal = b[orderBy];

      if (aVal == null && bVal == null) return 0;
      if (aVal == null) return 1;
      if (bVal == null) return -1;

      let comparison = 0;
      if (typeof aVal === 'string' && typeof bVal === 'string') {
        comparison = aVal.localeCompare(bVal);
      } else if (typeof aVal === 'number' && typeof bVal === 'number') {
        comparison = aVal - bVal;
      } else {
        comparison = String(aVal).localeCompare(String(bVal));
      }

      return order === 'asc' ? comparison : -comparison;
    });
  }, [filteredData, orderBy, order]);

  // Set up virtualizer for windowed rendering
  const rowVirtualizer = useVirtualizer({
    count: sortedData.length,
    getScrollElement: () => parentRef.current,
    estimateSize: () => effectiveRowHeight,
    overscan: 5, // Render 5 extra rows above/below viewport for smoother scrolling
  });

  const handleSort = (column: string) => {
    const isAsc = orderBy === column && order === 'asc';
    setOrder(isAsc ? 'desc' : 'asc');
    setOrderBy(column);
  };

  // Note: search still works across all columns (including hidden ones)
  const hasSearchableColumns = columns.some((col) => col.searchable) || additionalSearchFields.length > 0;

  const virtualItems = rowVirtualizer.getVirtualItems();

  return (
    <Paper sx={{ width: '100%', overflow: 'hidden' }}>
      {/* Header with title and search */}
      <Box sx={{ p: 2, display: 'flex', alignItems: 'center', gap: 2 }}>
        {title && (
          <Typography variant="h6" sx={{ flexGrow: 1 }}>
            {title}
          </Typography>
        )}
        {!title && <Box sx={{ flexGrow: 1 }} />}
        {headerAction}
        {hasSearchableColumns && (
          <TextField
            size="small"
            placeholder="Search..."
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
        )}
        <Typography variant="body2" color="text.secondary">
          {sortedData.length} {sortedData.length === 1 ? 'row' : 'rows'}
        </Typography>
      </Box>

      {/* Loading indicator */}
      {loading && <LinearProgress />}

      {/* Table with virtualized scrolling */}
      <TableContainer
        ref={parentRef}
        sx={{
          maxHeight,
          overflow: 'auto',
        }}
      >
        <Table stickyHeader sx={{ tableLayout: 'fixed' }}>
          <TableHead>
            <TableRow>
              {visibleColumns.map((column) => (
                <TableCell
                  key={column.key}
                  sx={{ fontWeight: 600, width: column.width }}
                >
                  {column.sortable ? (
                    <TableSortLabel
                      active={orderBy === column.key}
                      direction={orderBy === column.key ? order : 'asc'}
                      onClick={() => handleSort(column.key)}
                    >
                      {column.label}
                    </TableSortLabel>
                  ) : (
                    column.label
                  )}
                </TableCell>
              ))}
            </TableRow>
          </TableHead>
          <TableBody>
            {sortedData.length === 0 ? (
              <TableRow>
                <TableCell colSpan={visibleColumns.length} align="center" sx={{ py: 4 }}>
                  <Typography color="text.secondary">
                    {loading ? 'Loading...' : emptyMessage}
                  </Typography>
                </TableCell>
              </TableRow>
            ) : (
              <>
                {/* Spacer for rows above viewport */}
                {virtualItems.length > 0 && virtualItems[0].start > 0 && (
                  <TableRow>
                    <TableCell
                      colSpan={visibleColumns.length}
                      sx={{
                        height: virtualItems[0].start,
                        padding: 0,
                        border: 'none',
                      }}
                    />
                  </TableRow>
                )}

                {/* Virtualized rows */}
                {virtualItems.map((virtualRow) => {
                  const row = sortedData[virtualRow.index];
                  return (
                    <TableRow
                      key={getRowKey(row)}
                      data-index={virtualRow.index}
                      ref={rowVirtualizer.measureElement}
                      hover
                      onClick={() => onRowClick?.(row)}
                      sx={{
                        cursor: onRowClick ? 'pointer' : 'default',
                        '&:hover': onRowClick
                          ? { bgcolor: 'action.hover' }
                          : undefined,
                      }}
                    >
                      {visibleColumns.map((column) => (
                        <TableCell key={column.key} sx={cellPadding ? { py: cellPadding } : undefined}>
                          {column.render
                            ? column.render(row[column.key], row)
                            : row[column.key] ?? '-'}
                        </TableCell>
                      ))}
                    </TableRow>
                  );
                })}

                {/* Spacer for rows below viewport */}
                {virtualItems.length > 0 && (
                  <TableRow>
                    <TableCell
                      colSpan={visibleColumns.length}
                      sx={{
                        height:
                          rowVirtualizer.getTotalSize() -
                          (virtualItems[virtualItems.length - 1]?.end ?? 0),
                        padding: 0,
                        border: 'none',
                      }}
                    />
                  </TableRow>
                )}
              </>
            )}
          </TableBody>
        </Table>
      </TableContainer>
    </Paper>
  );
}
