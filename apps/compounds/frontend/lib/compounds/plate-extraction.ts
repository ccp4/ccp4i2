/**
 * Pure plate data extraction logic.
 *
 * Extracts dose-response data series from a 2D spreadsheet grid
 * using a PlateLayout definition. The grid must use absolute Excel
 * coordinates: cells[r][c] maps to Excel row r, column c (0-indexed).
 */

import type { PlateLayout, PlateFormat } from '@/types/compounds/models';

/** Extracted data series from a single plate row (or strip). */
export interface ExtractedSeries {
  row: number;              // Plate row (0-indexed)
  rowLetter: string;
  stripIndex?: number;      // Strip index within the row (0-indexed)
  compoundGroupIndex?: number;
  compoundName: string | null;
  dataValues: (number | null)[];   // [min_ctrl, data..., max_ctrl]
  minControlValues: (number | null)[];
  maxControlValues: (number | null)[];
  minControl: number | null;       // Averaged
  maxControl: number | null;       // Averaged
  startColumn: number;
  endColumn: number;
  hasIssues: boolean;
  issues: string[];
}

export type CellGrid = (string | number | null)[][];

const PLATE_DIMENSIONS: Record<PlateFormat, { rows: number; cols: number }> = {
  24: { rows: 4, cols: 6 },
  96: { rows: 8, cols: 12 },
  384: { rows: 16, cols: 24 },
  1536: { rows: 32, cols: 48 },
};

export function indexToRowLetter(index: number): string {
  return String.fromCharCode('A'.charCodeAt(0) + index);
}

export function excelColumnToIndex(col: string): number {
  let result = 0;
  for (let i = 0; i < col.length; i++) {
    result = result * 26 + (col.charCodeAt(i) - 'A'.charCodeAt(0) + 1);
  }
  return result - 1; // 0-indexed
}

function averageNonNull(values: (number | null)[]): number | null {
  const valid = values.filter((v): v is number => v !== null);
  return valid.length > 0
    ? valid.reduce((a, b) => a + b, 0) / valid.length
    : null;
}

/**
 * Extract data series from a spreadsheet grid using a plate layout.
 *
 * The cells grid MUST use absolute Excel coordinates:
 *   cells[r][c] corresponds to Excel row r+1, column index c.
 *   i.e. cells[0][0] = cell A1, cells[0][1] = cell B1, etc.
 */
export function extractPlateData(
  cells: CellGrid,
  plateLayout: PlateLayout,
): ExtractedSeries[] {
  const format = plateLayout.plate_format || 384;
  const { rows: numRows, cols: numCols } = PLATE_DIMENSIONS[format];

  // Get origin offset from spreadsheet_origin (absolute Excel coords)
  const originCol = excelColumnToIndex(plateLayout.spreadsheet_origin?.column || 'A');
  const originRow = (plateLayout.spreadsheet_origin?.row || 1) - 1;

  const series: ExtractedSeries[] = [];

  // Determine row range based on sample_region
  const startRowIdx = plateLayout.sample_region?.start_row
    ? plateLayout.sample_region.start_row.charCodeAt(0) - 'A'.charCodeAt(0)
    : 0;
  const endRowIdx = plateLayout.sample_region?.end_row
    ? plateLayout.sample_region.end_row.charCodeAt(0) - 'A'.charCodeAt(0)
    : numRows - 1;

  if (plateLayout.controls?.placement === 'per_compound' && plateLayout.strip_layout) {
    extractStripLayout(
      cells, plateLayout, series,
      originCol, originRow, startRowIdx, endRowIdx,
    );
  } else {
    extractEdgeControlLayout(
      cells, plateLayout, series,
      originCol, originRow, startRowIdx, endRowIdx, numCols,
    );
  }

  return series;
}

function extractStripLayout(
  cells: CellGrid,
  plateLayout: PlateLayout,
  series: ExtractedSeries[],
  originCol: number,
  originRow: number,
  startRowIdx: number,
  endRowIdx: number,
): void {
  const strip = plateLayout.strip_layout!;
  const stripWidth = strip.strip_width || 12;
  const stripsPerRow = strip.strips_per_row || 2;
  const replicateCount = plateLayout.replicate?.count || 1;
  const isExplicitNaming = plateLayout.replicate?.pattern === 'explicit';

  const compoundsPerRow = isExplicitNaming
    ? stripsPerRow
    : Math.max(1, Math.floor(stripsPerRow / replicateCount));

  for (let plateRow = startRowIdx; plateRow <= endRowIdx; plateRow++) {
    const spreadsheetRow = originRow + plateRow;
    if (spreadsheetRow >= cells.length) continue;

    const rowData = cells[spreadsheetRow] || [];
    const rowLetter = indexToRowLetter(plateRow);

    for (let stripIdx = 0; stripIdx < stripsPerRow; stripIdx++) {
      const compoundGroupIdx = isExplicitNaming
        ? stripIdx
        : Math.floor(stripIdx / replicateCount);
      const stripStartCol = originCol + stripIdx * stripWidth;

      // Resolve compound name
      let compoundName: string | null = null;

      if (plateLayout.compound_source?.compound_name_row) {
        const nameRowOffset = plateRow - startRowIdx;
        const nameRow = plateLayout.compound_source.compound_name_row - 1 + nameRowOffset;
        const nameCol = stripStartCol;
        if (nameRow >= 0 && nameRow < cells.length) {
          const nameRowData = cells[nameRow] || [];
          if (nameCol < nameRowData.length) {
            const nameVal = nameRowData[nameCol];
            compoundName = nameVal ? String(nameVal) : null;
          }
        }
      } else if (plateLayout.compound_source?.type === 'row_header') {
        const nameColIdx = originCol - compoundsPerRow + compoundGroupIdx;
        if (nameColIdx >= 0 && nameColIdx < rowData.length) {
          const nameVal = rowData[nameColIdx];
          compoundName = nameVal ? String(nameVal) : null;
        }
      } else if (plateLayout.compound_source?.type === 'adjacent_column') {
        const totalStripWidth = stripWidth * stripsPerRow;
        const nameColIdx = originCol + totalStripWidth + compoundGroupIdx;
        if (nameColIdx < rowData.length) {
          const nameVal = rowData[nameColIdx];
          compoundName = nameVal ? String(nameVal) : null;
        }
      } else if (plateLayout.compound_source?.type === 'row_order') {
        if (compoundsPerRow > 1) {
          compoundName = `${rowLetter}${compoundGroupIdx + 1}`;
        } else {
          compoundName = rowLetter;
        }
      }

      const issues: string[] = [];

      // Extract min controls (start of strip)
      const minControlValues: (number | null)[] = [];
      for (let i = 0; i < strip.min_wells; i++) {
        const cellIdx = stripStartCol + i;
        const val = cellIdx < rowData.length ? rowData[cellIdx] : null;
        minControlValues.push(typeof val === 'number' ? val : null);
      }

      // Extract data values (middle of strip)
      const dataValues: (number | null)[] = [];
      for (let i = 0; i < strip.data_wells; i++) {
        const cellIdx = stripStartCol + strip.min_wells + i;
        const val = cellIdx < rowData.length ? rowData[cellIdx] : null;
        dataValues.push(typeof val === 'number' ? val : null);
      }

      // Extract max controls (end of strip)
      const maxControlValues: (number | null)[] = [];
      for (let i = 0; i < strip.max_wells; i++) {
        const cellIdx = stripStartCol + strip.min_wells + strip.data_wells + i;
        const val = cellIdx < rowData.length ? rowData[cellIdx] : null;
        maxControlValues.push(typeof val === 'number' ? val : null);
      }

      const averagedMinControl = averageNonNull(minControlValues);
      const averagedMaxControl = averageNonNull(maxControlValues);

      if (dataValues.every(v => v === null)) {
        issues.push('No numeric data found');
      }
      if (averagedMinControl === null) {
        issues.push('Missing min controls');
      }
      if (averagedMaxControl === null) {
        issues.push('Missing max controls');
      }

      series.push({
        row: plateRow,
        rowLetter,
        stripIndex: stripIdx,
        compoundGroupIndex: compoundGroupIdx,
        compoundName,
        dataValues: [averagedMinControl, ...dataValues, averagedMaxControl],
        minControlValues,
        maxControlValues,
        minControl: averagedMinControl,
        maxControl: averagedMaxControl,
        startColumn: stripStartCol - originCol + 1,
        endColumn: stripStartCol - originCol + stripWidth,
        hasIssues: issues.length > 0,
        issues,
      });
    }
  }
}

function extractEdgeControlLayout(
  cells: CellGrid,
  plateLayout: PlateLayout,
  series: ExtractedSeries[],
  originCol: number,
  originRow: number,
  startRowIdx: number,
  endRowIdx: number,
  numCols: number,
): void {
  const startCol = plateLayout.sample_region?.start_column || 1;
  const endCol = plateLayout.sample_region?.end_column || numCols;

  for (let plateRow = startRowIdx; plateRow <= endRowIdx; plateRow++) {
    const spreadsheetRow = originRow + plateRow;
    if (spreadsheetRow >= cells.length) continue;

    const rowData = cells[spreadsheetRow] || [];
    const issues: string[] = [];

    // Extract min control values
    const minControlValues: (number | null)[] = [];
    if (plateLayout.controls?.min?.columns) {
      for (const col of plateLayout.controls.min.columns) {
        const cellIdx = originCol + col - 1;
        const val = cellIdx < rowData.length ? rowData[cellIdx] : null;
        minControlValues.push(typeof val === 'number' ? val : null);
      }
    }

    // Extract max control values
    const maxControlValues: (number | null)[] = [];
    if (plateLayout.controls?.max?.columns) {
      for (const col of plateLayout.controls.max.columns) {
        const cellIdx = originCol + col - 1;
        const val = cellIdx < rowData.length ? rowData[cellIdx] : null;
        maxControlValues.push(typeof val === 'number' ? val : null);
      }
    }

    // Extract data values
    const dataValues: (number | null)[] = [];
    for (let col = startCol; col <= endCol; col++) {
      const cellIdx = originCol + col - 1;
      const val = cellIdx < rowData.length ? rowData[cellIdx] : null;
      dataValues.push(typeof val === 'number' ? val : null);
    }

    // Compound name
    let compoundName: string | null = null;
    if (plateLayout.compound_source?.type === 'row_header') {
      const nameColIdx = originCol - 1;
      if (nameColIdx >= 0 && nameColIdx < rowData.length) {
        const nameVal = rowData[nameColIdx];
        compoundName = nameVal ? String(nameVal) : null;
      }
    } else if (plateLayout.compound_source?.type === 'adjacent_column') {
      const nameColIdx = originCol + endCol;
      if (nameColIdx < rowData.length) {
        const nameVal = rowData[nameColIdx];
        compoundName = nameVal ? String(nameVal) : null;
      }
    }

    const averagedMinControl = averageNonNull(minControlValues);
    const averagedMaxControl = averageNonNull(maxControlValues);

    if (dataValues.every(v => v === null)) {
      issues.push('No numeric data found');
    }

    series.push({
      row: plateRow,
      rowLetter: indexToRowLetter(plateRow),
      compoundName,
      dataValues: [averagedMinControl, ...dataValues, averagedMaxControl],
      minControlValues,
      maxControlValues,
      minControl: averagedMinControl,
      maxControl: averagedMaxControl,
      startColumn: startCol,
      endColumn: endCol,
      hasIssues: issues.length > 0,
      issues,
    });
  }
}
