/**
 * Utility for computing data series ranges from a PlateLayout.
 * Used by PlatePreview for visualization and by extraction logic.
 */

import type { PlateLayout, PlateFormat } from '@/types/compounds/models';

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
 * Represents the bounding box for a single data series on the plate.
 * All indices are 0-based for direct use in SVG rendering.
 */
export interface SeriesRange {
  /** Unique index for this series (0-based) */
  seriesIndex: number;
  /** Starting row (0-indexed) */
  rowStart: number;
  /** Ending row inclusive (0-indexed) */
  rowEnd: number;
  /** Starting column (0-indexed) */
  colStart: number;
  /** Ending column inclusive (0-indexed) */
  colEnd: number;
  /** Human-readable label (e.g., "A", "A-B", "A:1") */
  label: string;
  /** For strip layouts, which strip this is (0-indexed) */
  stripIndex?: number;
}

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
 * Compute the bounding boxes for all data series based on the plate layout.
 *
 * This function analyzes the PlateLayout configuration and returns an array
 * of SeriesRange objects, each representing the cell range that will be
 * extracted as one data series.
 *
 * For strip layouts (per_compound controls):
 *   - Each strip within each row becomes a separate series
 *   - The range includes min controls + data wells + max controls
 *
 * For standard edge control layouts:
 *   - Each row (or group of rows based on replicate pattern) becomes a series
 *   - The range spans from leftmost control column to rightmost control column
 */
export function computeSeriesRanges(layout: Partial<PlateLayout>): SeriesRange[] {
  const plateFormat = layout.plate_format || 96;
  const { rows: numRows, cols: numCols } = PLATE_DIMENSIONS[plateFormat];

  const ranges: SeriesRange[] = [];
  let seriesIndex = 0;

  // Determine row range from sample_region
  const startRowIdx = layout.sample_region?.start_row
    ? rowLetterToIndex(layout.sample_region.start_row)
    : 0;
  const endRowIdx = layout.sample_region?.end_row
    ? rowLetterToIndex(layout.sample_region.end_row)
    : numRows - 1;

  // Check if using strip layout with embedded controls
  if (layout.controls?.placement === 'per_compound' && layout.strip_layout) {
    // Strip layout: [min×N][data×M][max×N] repeated per row
    const strip = layout.strip_layout;
    const stripWidth = strip.strip_width;

    for (let plateRow = startRowIdx; plateRow <= endRowIdx; plateRow++) {
      for (let stripIdx = 0; stripIdx < strip.strips_per_row; stripIdx++) {
        const colStart = stripIdx * stripWidth;
        const colEnd = colStart + stripWidth - 1;

        ranges.push({
          seriesIndex,
          rowStart: plateRow,
          rowEnd: plateRow,
          colStart,
          colEnd,
          label: `${indexToRowLetter(plateRow)}:${stripIdx + 1}`,
          stripIndex: stripIdx,
        });
        seriesIndex++;
      }
    }
  } else {
    // Standard edge controls layout
    // Find the full column range including controls and sample region

    const sampleStartCol = (layout.sample_region?.start_column || 1) - 1; // to 0-indexed
    const sampleEndCol = (layout.sample_region?.end_column || numCols) - 1;

    // Find min and max of all involved columns
    const minControlCols = layout.controls?.min?.columns || [];
    const maxControlCols = layout.controls?.max?.columns || [];

    // Convert to 0-indexed
    const allCols = [
      sampleStartCol,
      sampleEndCol,
      ...minControlCols.map(c => c - 1),
      ...maxControlCols.map(c => c - 1),
    ];

    const colStart = Math.min(...allCols);
    const colEnd = Math.max(...allCols);

    // Handle replicate patterns to group rows
    const replicateCount = layout.replicate?.count || 1;
    const replicatePattern = layout.replicate?.pattern || 'adjacent_rows';

    if (replicatePattern === 'adjacent_rows' && replicateCount > 1) {
      // Adjacent rows: rows are grouped in sets of replicateCount
      // E.g., with count=2: A,B are one compound; C,D are another
      for (let plateRow = startRowIdx; plateRow <= endRowIdx; plateRow += replicateCount) {
        const rowEnd = Math.min(plateRow + replicateCount - 1, endRowIdx);
        const startLetter = indexToRowLetter(plateRow);
        const endLetter = indexToRowLetter(rowEnd);

        ranges.push({
          seriesIndex,
          rowStart: plateRow,
          rowEnd,
          colStart,
          colEnd,
          label: replicateCount > 1 && plateRow !== rowEnd
            ? `${startLetter}-${endLetter}`
            : startLetter,
        });
        seriesIndex++;
      }
    } else if (replicatePattern === 'grouped_rows' && replicateCount > 1) {
      // Grouped rows: first N rows are replicates of each other
      // E.g., with count=2: row A compound 1-8, row B is replicate of A
      // Each compound spans replicateCount rows at the same column positions
      // This is more complex - for now, show each row separately but with group indicator
      const compoundsPerGroup = Math.ceil((endRowIdx - startRowIdx + 1) / replicateCount);

      for (let compoundIdx = 0; compoundIdx < compoundsPerGroup; compoundIdx++) {
        const groupRows: number[] = [];
        for (let repIdx = 0; repIdx < replicateCount; repIdx++) {
          const row = startRowIdx + compoundIdx + repIdx * compoundsPerGroup;
          if (row <= endRowIdx) {
            groupRows.push(row);
          }
        }

        if (groupRows.length > 0) {
          // For grouped rows, show each replicate row separately but with same series index
          // This creates multiple ranges per compound, which is accurate
          for (const row of groupRows) {
            ranges.push({
              seriesIndex,
              rowStart: row,
              rowEnd: row,
              colStart,
              colEnd,
              label: `${indexToRowLetter(row)} (${seriesIndex + 1})`,
            });
          }
          seriesIndex++;
        }
      }
    } else {
      // No replicates or other patterns: one series per row
      for (let plateRow = startRowIdx; plateRow <= endRowIdx; plateRow++) {
        ranges.push({
          seriesIndex,
          rowStart: plateRow,
          rowEnd: plateRow,
          colStart,
          colEnd,
          label: indexToRowLetter(plateRow),
        });
        seriesIndex++;
      }
    }
  }

  return ranges;
}

/**
 * Generate a color for a series based on its index.
 * Uses a palette of distinct, visually pleasing colors.
 */
export function getSeriesColor(seriesIndex: number, alpha: number = 0.3): string {
  // Palette of distinct colors that work well with transparency
  const palette = [
    [255, 152, 0],   // Orange
    [156, 39, 176],  // Purple
    [0, 188, 212],   // Cyan
    [255, 87, 34],   // Deep Orange
    [63, 81, 181],   // Indigo
    [0, 150, 136],   // Teal
    [233, 30, 99],   // Pink
    [139, 195, 74],  // Light Green
    [121, 85, 72],   // Brown
    [255, 193, 7],   // Amber
  ];

  const color = palette[seriesIndex % palette.length];
  return `rgba(${color[0]}, ${color[1]}, ${color[2]}, ${alpha})`;
}

/**
 * Get a border color (more saturated) for a series
 */
export function getSeriesBorderColor(seriesIndex: number): string {
  return getSeriesColor(seriesIndex, 0.8);
}
