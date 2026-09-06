"use client";

import { useState, useMemo, useRef, ReactNode } from "react";
import {
  Box,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  TableRowProps,
  TableSortLabel,
  Typography,
} from "@mui/material";
import { useVirtualizer } from "@tanstack/react-virtual";

export interface Column<T> {
  key: string;
  label: string;
  sortable?: boolean;
  searchable?: boolean;
  render?: (value: any, row: T) => ReactNode;
  header?: ReactNode;
  width?: string | number;
}

interface DataTableProps<T> {
  data: T[];
  columns: Column<T>[];
  onRowClick?: (row: T) => void;
  getRowKey: (row: T) => string | number;
  /**
   * Extra props for a row — drag-and-drop handlers, a highlight — merged over
   * the table's own. An `sx` here is layered on top of the row's default.
   */
  getRowProps?: (row: T) => Partial<TableRowProps>;
  emptyMessage: string;
}

type Order = "asc" | "desc";

export function DataTable<T extends Record<string, any>>({
  data,
  columns,
  onRowClick,
  getRowKey,
  getRowProps,
  emptyMessage,
}: DataTableProps<T>) {
  const [orderBy, setOrderBy] = useState<string | null>(null);
  const [order, setOrder] = useState<Order>("asc");

  const parentRef = useRef<HTMLDivElement>(null);

  // Sort data
  const sortedData = useMemo(() => {
    if (!orderBy) return data;

    return [...data].sort((a, b) => {
      const aVal = a[orderBy];
      const bVal = b[orderBy];

      if (aVal == null && bVal == null) return 0;
      if (aVal == null) return 1;
      if (bVal == null) return -1;

      let comparison = 0;
      if (typeof aVal === "string" && typeof bVal === "string") {
        comparison = aVal.localeCompare(bVal);
      } else if (typeof aVal === "number" && typeof bVal === "number") {
        comparison = aVal - bVal;
      } else {
        comparison = String(aVal).localeCompare(String(bVal));
      }

      return order === "asc" ? comparison : -comparison;
    });
  }, [data, orderBy, order]);

  // Set up virtualizer for windowed rendering
  const rowVirtualizer = useVirtualizer({
    count: sortedData.length,
    getScrollElement: () => parentRef.current,
    estimateSize: () => 60,
    overscan: 5, // Render 5 extra rows above/below viewport for smoother scrolling
  });

  const handleSort = (column: string) => {
    const isAsc = orderBy === column && order === "asc";
    setOrder(isAsc ? "desc" : "asc");
    setOrderBy(column);
  };

  const virtualItems = rowVirtualizer.getVirtualItems();

  return (
    <Box
      sx={{
        width: "100%",
        overflow: "hidden",
        height: "100%",
        display: "flex",
        flexDirection: "column",
      }}
    >
      {/* Table with virtualized scrolling */}
      <TableContainer
        ref={parentRef}
        sx={{
          flex: 1,
          overflow: "auto",
        }}
      >
        <Table stickyHeader size="small" sx={{ tableLayout: "fixed" }}>
          <TableHead>
            <TableRow>
              {columns.map((column) => (
                <TableCell
                  key={column.key}
                  sx={{ fontWeight: 600, width: column.width }}
                >
                  {column.header ? (
                    column.header
                  ) : column.sortable ? (
                    <TableSortLabel
                      active={orderBy === column.key}
                      direction={orderBy === column.key ? order : "asc"}
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
                <TableCell
                  colSpan={columns.length}
                  align="center"
                  sx={{ py: 4 }}
                >
                  <Typography color="text.secondary">{emptyMessage}</Typography>
                </TableCell>
              </TableRow>
            ) : (
              <>
                {/* Spacer for rows above viewport */}
                {virtualItems.length > 0 && virtualItems[0].start > 0 && (
                  <TableRow>
                    <TableCell
                      colSpan={columns.length}
                      sx={{
                        height: virtualItems[0].start,
                        padding: 0,
                        border: "none",
                      }}
                    />
                  </TableRow>
                )}

                {/* Virtualized rows */}
                {virtualItems.map((virtualRow) => {
                  const row = sortedData[virtualRow.index];
                  const { sx: extraSx, ...extraProps } = getRowProps?.(row) ?? {};
                  return (
                    <TableRow
                      key={getRowKey(row)}
                      data-index={virtualRow.index}
                      ref={rowVirtualizer.measureElement}
                      hover
                      onClick={() => onRowClick?.(row)}
                      {...extraProps}
                      sx={[
                        {
                          cursor: onRowClick ? "pointer" : "default",
                          "&:hover": onRowClick
                            ? { bgcolor: "action.hover" }
                            : undefined,
                        },
                        ...(Array.isArray(extraSx) ? extraSx : [extraSx]),
                      ]}
                    >
                      {columns.map((column) => (
                        <TableCell key={column.key}>
                          {column.render
                            ? column.render(row[column.key], row)
                            : (row[column.key] ?? "-")}
                        </TableCell>
                      ))}
                    </TableRow>
                  );
                })}

                {/* Spacer for rows below viewport */}
                {virtualItems.length > 0 && (
                  <TableRow>
                    <TableCell
                      colSpan={columns.length}
                      sx={{
                        height:
                          rowVirtualizer.getTotalSize() -
                          (virtualItems[virtualItems.length - 1]?.end ?? 0),
                        padding: 0,
                        border: "none",
                      }}
                    />
                  </TableRow>
                )}
              </>
            )}
          </TableBody>
        </Table>
      </TableContainer>
    </Box>
  );
}
