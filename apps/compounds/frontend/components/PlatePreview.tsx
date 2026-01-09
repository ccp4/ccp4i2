'use client';

import { useMemo } from 'react';
import { Box, Typography, Paper } from '@mui/material';
import type { PlateLayout, PlateFormat } from '@/types/models';

/**
 * Plate dimensions for different formats
 */
const PLATE_DIMENSIONS: Record<PlateFormat, { rows: number; cols: number }> = {
  24: { rows: 4, cols: 6 },
  96: { rows: 8, cols: 12 },
  384: { rows: 16, cols: 24 },
  1536: { rows: 32, cols: 48 },
};

/**
 * Convert row letter to 0-indexed number
 */
function rowLetterToIndex(letter: string): number {
  return letter.charCodeAt(0) - 'A'.charCodeAt(0);
}

/**
 * Convert 0-indexed number to row letter
 */
function indexToRowLetter(index: number): string {
  return String.fromCharCode('A'.charCodeAt(0) + index);
}

/**
 * Well type for coloring
 */
type WellType = 'empty' | 'max_control' | 'min_control' | 'sample';

interface PlatePreviewProps {
  layout: Partial<PlateLayout>;
  width?: number;
  height?: number;
  showLabels?: boolean;
  showLegend?: boolean;
}

/**
 * Visual preview of a plate layout configuration.
 * Renders an SVG grid showing control wells, sample regions, and empty wells.
 */
export function PlatePreview({
  layout,
  width = 400,
  height = 280,
  showLabels = true,
  showLegend = true,
}: PlatePreviewProps) {
  const plateFormat = layout.plate_format || 96;
  const { rows: numRows, cols: numCols } = PLATE_DIMENSIONS[plateFormat];

  // Calculate well positions and types
  const wells = useMemo(() => {
    const wellData: { row: number; col: number; type: WellType }[] = [];

    // Check if using strip layout with embedded controls
    const isStripLayout = layout.controls?.placement === 'per_compound' && layout.strip_layout;

    for (let row = 0; row < numRows; row++) {
      for (let col = 0; col < numCols; col++) {
        let type: WellType = 'empty';
        const colNum = col + 1; // 1-indexed
        const rowLetter = indexToRowLetter(row);

        if (isStripLayout && layout.strip_layout) {
          // Strip layout: [min×N][data×M][max×N] repeated
          const { strip_width, min_wells, data_wells, max_wells } = layout.strip_layout;
          const posInStrip = (col % strip_width) + 1; // 1-indexed position within strip

          // Check if this row is in the sample region
          const inSampleRows = layout.sample_region
            ? rowLetter >= layout.sample_region.start_row && rowLetter <= layout.sample_region.end_row
            : true;

          if (inSampleRows) {
            if (posInStrip <= min_wells) {
              type = 'min_control';
            } else if (posInStrip <= min_wells + data_wells) {
              type = 'sample';
            } else if (posInStrip <= min_wells + data_wells + max_wells) {
              type = 'max_control';
            }
          }
        } else {
          // Standard layout: edge controls + sample region

          // Check if this well is a max control
          if (layout.controls?.max) {
            const { rows, columns } = layout.controls.max;
            if (rows?.includes(rowLetter) && columns?.includes(colNum)) {
              type = 'max_control';
            }
          }

          // Check if this well is a min control
          if (layout.controls?.min && type === 'empty') {
            const { rows, columns } = layout.controls.min;
            if (rows?.includes(rowLetter) && columns?.includes(colNum)) {
              type = 'min_control';
            }
          }

          // Check if this well is in the sample region
          if (layout.sample_region && type === 'empty') {
            const { start_column, end_column, start_row, end_row } = layout.sample_region;

            if (
              colNum >= start_column &&
              colNum <= end_column &&
              rowLetter >= start_row &&
              rowLetter <= end_row
            ) {
              type = 'sample';
            }
          }
        }

        wellData.push({ row, col, type });
      }
    }

    return wellData;
  }, [layout, numRows, numCols]);

  // Layout calculations
  const labelWidth = showLabels ? 20 : 0;
  const labelHeight = showLabels ? 16 : 0;
  const legendHeight = showLegend ? 30 : 0;

  const plateWidth = width - labelWidth - 10;
  const plateHeight = height - labelHeight - legendHeight - 10;

  const wellWidth = plateWidth / numCols;
  const wellHeight = plateHeight / numRows;
  const wellRadius = Math.min(wellWidth, wellHeight) * 0.4;

  // Colors for well types
  const wellColors: Record<WellType, string> = {
    empty: '#f5f5f5',
    max_control: '#2196f3',    // Blue
    min_control: '#f44336',    // Red
    sample: '#4caf50',         // Green
  };

  const wellBorders: Record<WellType, string> = {
    empty: '#e0e0e0',
    max_control: '#1976d2',
    min_control: '#d32f2f',
    sample: '#388e3c',
  };

  return (
    <Paper sx={{ p: 2, bgcolor: 'grey.50' }}>
      <svg width={width} height={height} viewBox={`0 0 ${width} ${height}`}>
        {/* Column labels */}
        {showLabels && Array.from({ length: numCols }, (_, i) => (
          <text
            key={`col-${i}`}
            x={labelWidth + wellWidth * i + wellWidth / 2}
            y={12}
            textAnchor="middle"
            fontSize={numCols > 24 ? 6 : 8}
            fill="#666"
          >
            {i + 1}
          </text>
        ))}

        {/* Row labels */}
        {showLabels && Array.from({ length: numRows }, (_, i) => (
          <text
            key={`row-${i}`}
            x={10}
            y={labelHeight + wellHeight * i + wellHeight / 2 + 3}
            textAnchor="middle"
            fontSize={numRows > 16 ? 6 : 8}
            fill="#666"
          >
            {indexToRowLetter(i)}
          </text>
        ))}

        {/* Wells */}
        {wells.map(({ row, col, type }) => (
          <circle
            key={`${row}-${col}`}
            cx={labelWidth + wellWidth * col + wellWidth / 2}
            cy={labelHeight + wellHeight * row + wellHeight / 2}
            r={wellRadius}
            fill={wellColors[type]}
            stroke={wellBorders[type]}
            strokeWidth={1}
          />
        ))}

        {/* Legend */}
        {showLegend && (
          <g transform={`translate(${labelWidth}, ${height - 25})`}>
            <circle cx={8} cy={8} r={6} fill={wellColors.max_control} stroke={wellBorders.max_control} />
            <text x={20} y={12} fontSize={9} fill="#333">Max</text>

            <circle cx={58} cy={8} r={6} fill={wellColors.min_control} stroke={wellBorders.min_control} />
            <text x={70} y={12} fontSize={9} fill="#333">Min</text>

            <circle cx={108} cy={8} r={6} fill={wellColors.sample} stroke={wellBorders.sample} />
            <text x={120} y={12} fontSize={9} fill="#333">Sample</text>

            <circle cx={178} cy={8} r={6} fill={wellColors.empty} stroke={wellBorders.empty} />
            <text x={190} y={12} fontSize={9} fill="#333">Empty</text>
          </g>
        )}
      </svg>

      {/* Summary stats */}
      <Box sx={{ mt: 1, display: 'flex', gap: 2, justifyContent: 'center' }}>
        <Typography variant="caption" color="text.secondary">
          {plateFormat}-well plate
        </Typography>
        <Typography variant="caption" color="text.secondary">
          {wells.filter(w => w.type === 'sample').length} sample wells
        </Typography>
        <Typography variant="caption" color="text.secondary">
          {wells.filter(w => w.type === 'max_control').length + wells.filter(w => w.type === 'min_control').length} control wells
        </Typography>
      </Box>
    </Paper>
  );
}

/**
 * Compact plate preview for use in tables or cards
 */
interface PlatePreviewThumbProps {
  layout: Partial<PlateLayout>;
  size?: number;
}

export function PlatePreviewThumb({ layout, size = 80 }: PlatePreviewThumbProps) {
  return (
    <PlatePreview
      layout={layout}
      width={size}
      height={size * 0.7}
      showLabels={false}
      showLegend={false}
    />
  );
}
