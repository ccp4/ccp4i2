/*
 * Copyright (C) 2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
'use client';

import { useMemo } from 'react';
import { Box, Typography, Paper } from '@mui/material';
import type { PlateLayout, PlateFormat } from '@/types/compounds/models';
import { computeWellMap, WellType, indexToRowLetter } from '@/lib/compounds/plate-extraction';
import {
  computeSeriesRanges,
  getSeriesColor,
  getSeriesBorderColor,
} from '@/lib/compounds/plateSeriesRanges';

/**
 * Plate dimensions for different formats
 */
const PLATE_DIMENSIONS: Record<PlateFormat, { rows: number; cols: number }> = {
  24: { rows: 4, cols: 6 },
  96: { rows: 8, cols: 12 },
  384: { rows: 16, cols: 24 },
  1536: { rows: 32, cols: 48 },
};

// Colors for well types
const WELL_COLORS: Record<WellType, string> = {
  empty: '#f5f5f5',
  max_control: '#2196f3',
  min_control: '#f44336',
  sample: '#4caf50',
};

const WELL_BORDERS: Record<WellType, string> = {
  empty: '#e0e0e0',
  max_control: '#1976d2',
  min_control: '#d32f2f',
  sample: '#388e3c',
};

interface PlatePreviewProps {
  layout: Partial<PlateLayout>;
  width?: number;
  height?: number;
  showLabels?: boolean;
  showLegend?: boolean;
  /** Show bounding boxes around each data series range */
  showSeriesRanges?: boolean;
}

/**
 * Visual preview of a plate layout configuration.
 * Renders an SVG grid showing control wells, sample regions, and empty wells.
 * Well classification is derived from the same logic used by the data extraction
 * code, ensuring the preview always matches what extraction would actually do.
 */
export function PlatePreview({
  layout,
  width = 400,
  height = 280,
  showLabels = true,
  showLegend = true,
  showSeriesRanges = false,
}: PlatePreviewProps) {
  const plateFormat = layout.plate_format || 96;
  const { rows: numRows, cols: numCols } = PLATE_DIMENSIONS[plateFormat];

  // Well classification — uses the same logic as plate-extraction.ts
  const wellMap = useMemo(() => computeWellMap(layout), [layout]);

  // Compute series ranges for visualization
  const seriesRanges = useMemo(() => {
    if (!showSeriesRanges) return [];
    return computeSeriesRanges(layout);
  }, [layout, showSeriesRanges]);

  // Layout calculations
  const labelWidth = showLabels ? 20 : 0;
  const labelHeight = showLabels ? 16 : 0;
  const legendHeight = showLegend ? 30 : 0;

  const plateWidth = width - labelWidth - 10;
  const plateHeight = height - labelHeight - legendHeight - 10;

  const wellWidth = plateWidth / numCols;
  const wellHeight = plateHeight / numRows;
  const wellRadius = Math.min(wellWidth, wellHeight) * 0.4;

  // Strip layout metadata for visual indicators
  const isStripLayout = layout.controls?.placement === 'per_compound' && layout.strip_layout;
  const stripWidth = layout.strip_layout?.strip_width || 0;
  const stripsPerRow = layout.strip_layout?.strips_per_row || 0;

  // Replicate grouping for alternating bands
  const replicateCount = layout.replicate?.count || 1;
  const replicatePattern = layout.replicate?.pattern || 'adjacent_rows';
  const startRowIdx = layout.sample_region?.start_row
    ? layout.sample_region.start_row.charCodeAt(0) - 'A'.charCodeAt(0)
    : 0;
  const endRowIdx = layout.sample_region?.end_row
    ? layout.sample_region.end_row.charCodeAt(0) - 'A'.charCodeAt(0)
    : numRows - 1;

  // Dilution direction
  const dilutionDir = layout.dilution?.direction || 'horizontal';

  // Count wells by type for summary
  const wellCounts = useMemo(() => {
    let sample = 0, control = 0;
    for (let r = 0; r < wellMap.length; r++) {
      for (let c = 0; c < wellMap[r].length; c++) {
        if (wellMap[r][c] === 'sample') sample++;
        else if (wellMap[r][c] === 'max_control' || wellMap[r][c] === 'min_control') control++;
      }
    }
    return { sample, control };
  }, [wellMap]);

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

        {/* Replicate group bands (alternating shading) */}
        {replicatePattern === 'adjacent_rows' && replicateCount > 1 && (() => {
          const bands: React.ReactNode[] = [];
          let groupIdx = 0;
          for (let row = startRowIdx; row <= endRowIdx; row += replicateCount) {
            if (groupIdx % 2 === 0) {
              const bandEnd = Math.min(row + replicateCount - 1, endRowIdx);
              const y = labelHeight + wellHeight * row;
              const bandHeight = wellHeight * (bandEnd - row + 1);
              bands.push(
                <rect
                  key={`band-${row}`}
                  x={labelWidth}
                  y={y}
                  width={plateWidth}
                  height={bandHeight}
                  fill="rgba(0, 0, 0, 0.04)"
                  rx={2}
                />
              );
            }
            groupIdx++;
          }
          return bands;
        })()}

        {/* Strip boundary lines */}
        {isStripLayout && stripsPerRow > 1 && Array.from({ length: stripsPerRow - 1 }, (_, i) => {
          const colBoundary = (i + 1) * stripWidth;
          if (colBoundary >= numCols) return null;
          const x = labelWidth + wellWidth * colBoundary;
          return (
            <line
              key={`strip-boundary-${i}`}
              x1={x}
              y1={labelHeight + wellHeight * startRowIdx}
              x2={x}
              y2={labelHeight + wellHeight * (endRowIdx + 1)}
              stroke="#9e9e9e"
              strokeWidth={1}
              strokeDasharray="4,3"
            />
          );
        })}

        {/* Series range boxes (rendered before wells so they appear behind) */}
        {showSeriesRanges && seriesRanges.map((range) => {
          const x = labelWidth + wellWidth * range.colStart;
          const y = labelHeight + wellHeight * range.rowStart;
          const boxWidth = wellWidth * (range.colEnd - range.colStart + 1);
          const boxHeight = wellHeight * (range.rowEnd - range.rowStart + 1);
          const padding = 2;

          return (
            <g key={`series-${range.seriesIndex}-${range.rowStart}-${range.colStart}`}>
              <rect
                x={x - padding}
                y={y - padding}
                width={boxWidth + padding * 2}
                height={boxHeight + padding * 2}
                fill={getSeriesColor(range.seriesIndex, 0.15)}
                stroke={getSeriesBorderColor(range.seriesIndex)}
                strokeWidth={1.5}
                rx={3}
                ry={3}
              />
              {/* Series label in corner */}
              <text
                x={x + 2}
                y={y + 8}
                fontSize={Math.min(wellHeight * 0.6, 9)}
                fontWeight="bold"
                fill={getSeriesBorderColor(range.seriesIndex)}
              >
                {range.label}
              </text>
            </g>
          );
        })}

        {/* Wells */}
        {wellMap.map((rowData, row) =>
          rowData.map((type, col) => (
            <circle
              key={`${row}-${col}`}
              cx={labelWidth + wellWidth * col + wellWidth / 2}
              cy={labelHeight + wellHeight * row + wellHeight / 2}
              r={wellRadius}
              fill={WELL_COLORS[type]}
              stroke={WELL_BORDERS[type]}
              strokeWidth={1}
            />
          ))
        )}

        {/* Dilution direction arrow */}
        {(() => {
          // Find the sample region for arrow placement
          const sampleStartCol = (layout.sample_region?.start_column || 1) - 1;
          const sampleEndCol = (layout.sample_region?.end_column || numCols) - 1;

          if (dilutionDir === 'horizontal') {
            // Arrow across columns below the plate
            const arrowY = labelHeight + wellHeight * (endRowIdx + 1) + 6;
            const arrowX1 = labelWidth + wellWidth * sampleStartCol + wellWidth / 2;
            const arrowX2 = labelWidth + wellWidth * sampleEndCol + wellWidth / 2;

            // For strip layouts, draw arrow per strip
            if (isStripLayout) {
              const minWells = layout.strip_layout?.min_wells || 0;
              const dataWells = layout.strip_layout?.data_wells || 0;
              return Array.from({ length: stripsPerRow }, (_, stripIdx) => {
                const stripStart = stripIdx * stripWidth;
                const dataStart = stripStart + minWells;
                const dataEnd = dataStart + dataWells - 1;
                const x1 = labelWidth + wellWidth * dataStart + wellWidth / 2;
                const x2 = labelWidth + wellWidth * dataEnd + wellWidth / 2;
                return (
                  <g key={`arrow-${stripIdx}`}>
                    <line x1={x1} y1={arrowY} x2={x2} y2={arrowY} stroke="#757575" strokeWidth={1.5} />
                    <polygon
                      points={`${x2},${arrowY - 3} ${x2 + 5},${arrowY} ${x2},${arrowY + 3}`}
                      fill="#757575"
                    />
                    <text x={(x1 + x2) / 2} y={arrowY + 10} textAnchor="middle" fontSize={7} fill="#757575">
                      dilution
                    </text>
                  </g>
                );
              });
            }

            return (
              <g key="arrow">
                <line x1={arrowX1} y1={arrowY} x2={arrowX2} y2={arrowY} stroke="#757575" strokeWidth={1.5} />
                <polygon
                  points={`${arrowX2},${arrowY - 3} ${arrowX2 + 5},${arrowY} ${arrowX2},${arrowY + 3}`}
                  fill="#757575"
                />
                <text x={(arrowX1 + arrowX2) / 2} y={arrowY + 10} textAnchor="middle" fontSize={7} fill="#757575">
                  dilution
                </text>
              </g>
            );
          } else {
            // Vertical dilution arrow to the right of the plate
            const arrowX = labelWidth + plateWidth + 4;
            const arrowY1 = labelHeight + wellHeight * startRowIdx + wellHeight / 2;
            const arrowY2 = labelHeight + wellHeight * endRowIdx + wellHeight / 2;
            return (
              <g key="arrow">
                <line x1={arrowX} y1={arrowY1} x2={arrowX} y2={arrowY2} stroke="#757575" strokeWidth={1.5} />
                <polygon
                  points={`${arrowX - 3},${arrowY2} ${arrowX},${arrowY2 + 5} ${arrowX + 3},${arrowY2}`}
                  fill="#757575"
                />
              </g>
            );
          }
        })()}

        {/* Legend */}
        {showLegend && (
          <g transform={`translate(${labelWidth}, ${height - 25})`}>
            <circle cx={8} cy={8} r={6} fill={WELL_COLORS.max_control} stroke={WELL_BORDERS.max_control} />
            <text x={20} y={12} fontSize={9} fill="#333">Max</text>

            <circle cx={58} cy={8} r={6} fill={WELL_COLORS.min_control} stroke={WELL_BORDERS.min_control} />
            <text x={70} y={12} fontSize={9} fill="#333">Min</text>

            <circle cx={108} cy={8} r={6} fill={WELL_COLORS.sample} stroke={WELL_BORDERS.sample} />
            <text x={120} y={12} fontSize={9} fill="#333">Sample</text>

            <circle cx={178} cy={8} r={6} fill={WELL_COLORS.empty} stroke={WELL_BORDERS.empty} />
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
          {wellCounts.sample} sample wells
        </Typography>
        <Typography variant="caption" color="text.secondary">
          {wellCounts.control} control wells
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
