import { describe, it, expect } from 'vitest';
import * as path from 'path';
import {
  extractPlateData,
  excelColumnToIndex,
  indexToRowLetter,
} from '../plate-extraction';
import type { CellGrid } from '../plate-extraction';
import type { PlateLayout } from '@/types/compounds/models';

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/** Build a CellGrid from absolute Excel coordinates, padding with nulls. */
function buildGrid(
  entries: Record<string, string | number>,
): CellGrid {
  let maxRow = 0;
  let maxCol = 0;
  const parsed: { r: number; c: number; v: string | number }[] = [];

  for (const [addr, value] of Object.entries(entries)) {
    const match = addr.match(/^([A-Z]+)(\d+)$/);
    if (!match) throw new Error(`Bad cell address: ${addr}`);
    const c = excelColumnToIndex(match[1]);
    const r = Number(match[2]) - 1; // 0-indexed
    if (r > maxRow) maxRow = r;
    if (c > maxCol) maxCol = c;
    parsed.push({ r, c, v: value });
  }

  const grid: CellGrid = [];
  for (let row = 0; row <= maxRow; row++) {
    const rowData: (string | number | null)[] = [];
    for (let col = 0; col <= maxCol; col++) {
      rowData.push(null);
    }
    grid.push(rowData);
  }
  for (const { r, c, v } of parsed) {
    grid[r][c] = v;
  }
  return grid;
}

// ---------------------------------------------------------------------------
// Unit: helper functions
// ---------------------------------------------------------------------------

describe('excelColumnToIndex', () => {
  it('converts single letters', () => {
    expect(excelColumnToIndex('A')).toBe(0);
    expect(excelColumnToIndex('B')).toBe(1);
    expect(excelColumnToIndex('Z')).toBe(25);
  });
  it('converts multi-letter columns', () => {
    expect(excelColumnToIndex('AA')).toBe(26);
    expect(excelColumnToIndex('AB')).toBe(27);
  });
});

describe('indexToRowLetter', () => {
  it('converts indices to letters', () => {
    expect(indexToRowLetter(0)).toBe('A');
    expect(indexToRowLetter(9)).toBe('J');
    expect(indexToRowLetter(15)).toBe('P');
  });
});

// ---------------------------------------------------------------------------
// Unit: edge-control layout extraction
// ---------------------------------------------------------------------------

describe('extractPlateData — edge controls', () => {
  const layout: PlateLayout = {
    plate_format: 96,
    controls: {
      placement: 'edge_columns',
      max: { columns: [1], rows: ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'] },
      min: { columns: [12], rows: ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'] },
    },
    sample_region: {
      start_column: 2,
      end_column: 11,
      start_row: 'A',
      end_row: 'B',
    },
    dilution: { direction: 'horizontal', num_concentrations: 10 },
    replicate: { count: 1, pattern: 'adjacent_rows' },
    compound_source: { type: 'row_order' },
    spreadsheet_origin: { column: 'A', row: 1 },
  };

  it('extracts correct number of series', () => {
    // Build a small 2-row × 12-col grid at A1
    const entries: Record<string, number> = {};
    for (let r = 0; r < 2; r++) {
      for (let c = 0; c < 12; c++) {
        const addr = String.fromCharCode(65 + c) + (r + 1);
        entries[addr] = 100 + r * 12 + c;
      }
    }
    const grid = buildGrid(entries);
    const series = extractPlateData(grid, layout);

    expect(series).toHaveLength(2);
    expect(series[0].rowLetter).toBe('A');
    expect(series[1].rowLetter).toBe('B');
  });

  it('reads correct control and data values', () => {
    const entries: Record<string, number> = {};
    // Row A: max=500, data=10..19, min=999
    entries['A1'] = 500;
    for (let i = 0; i < 10; i++) {
      entries[String.fromCharCode(66 + i) + '1'] = 10 + i;
    }
    entries['L1'] = 999;
    const grid = buildGrid(entries);

    const series = extractPlateData(grid, { ...layout, sample_region: { ...layout.sample_region, end_row: 'A' } });
    expect(series).toHaveLength(1);
    const s = series[0];

    expect(s.maxControlValues).toEqual([500]);
    expect(s.maxControl).toBe(500);
    expect(s.minControlValues).toEqual([999]);
    expect(s.minControl).toBe(999);
    // dataValues = [min_avg, data0..data9, max_avg]
    expect(s.dataValues[0]).toBe(999);  // min at position 0
    expect(s.dataValues[s.dataValues.length - 1]).toBe(500); // max at end
  });

  it('handles non-zero spreadsheet_origin', () => {
    // Data starts at C5 (col=2, row=4 in 0-indexed)
    const entries: Record<string, number> = {};
    for (let c = 0; c < 12; c++) {
      const col = String.fromCharCode(67 + c); // C, D, E, ...
      entries[col + '5'] = 200 + c;
    }
    const grid = buildGrid(entries);

    const offsetLayout: PlateLayout = {
      ...layout,
      sample_region: { ...layout.sample_region, end_row: 'A' },
      spreadsheet_origin: { column: 'C', row: 5 },
    };
    const series = extractPlateData(grid, offsetLayout);

    expect(series).toHaveLength(1);
    expect(series[0].maxControlValues).toEqual([200]); // C5
    expect(series[0].minControlValues).toEqual([211]); // N5 (C+11)
  });

  it('works correctly when grid built from non-zero range.s.r', () => {
    // Simulate what used to happen: grid starts at row 4 (0-indexed)
    // but with our fix, cells[0..3] are null-filled.
    // Data at row 5 (Excel row 5 = index 4)
    const entries: Record<string, number> = {};
    for (let c = 0; c < 12; c++) {
      entries[String.fromCharCode(65 + c) + '5'] = 300 + c;
    }
    const grid = buildGrid(entries);
    // Verify the grid has rows 0-4 (row 5 = index 4)
    expect(grid.length).toBe(5);
    expect(grid[0].every(v => v === null)).toBe(true); // rows 1-4 are null
    expect(grid[4][0]).toBe(300); // row 5 has data

    const offsetLayout: PlateLayout = {
      ...layout,
      sample_region: { ...layout.sample_region, end_row: 'A' },
      spreadsheet_origin: { column: 'A', row: 5 },
    };
    const series = extractPlateData(grid, offsetLayout);

    expect(series).toHaveLength(1);
    expect(series[0].maxControlValues).toEqual([300]); // A5
  });
});

// ---------------------------------------------------------------------------
// Unit: strip layout extraction
// ---------------------------------------------------------------------------

describe('extractPlateData — strip layout', () => {
  const stripLayout: PlateLayout = {
    plate_format: 384,
    controls: {
      placement: 'per_compound',
      max: { columns: [], rows: [] },
      min: { columns: [], rows: [] },
    },
    sample_region: {
      start_column: 1,
      end_column: 24,
      start_row: 'A',
      end_row: 'B',
    },
    dilution: { direction: 'horizontal', num_concentrations: 10 },
    replicate: { count: 2, pattern: 'adjacent_rows' },
    compound_source: { type: 'adjacent_column' },
    strip_layout: {
      strip_width: 12,
      min_wells: 1,
      data_wells: 10,
      max_wells: 1,
      strips_per_row: 2,
    },
    spreadsheet_origin: { column: 'B', row: 1 },
  };

  it('produces 2 strips per row', () => {
    const entries: Record<string, string | number> = {};
    for (let r = 0; r < 2; r++) {
      for (let c = 1; c <= 24; c++) {
        const col = String.fromCharCode(65 + c); // B=1, C=2, ...
        entries[col + (r + 1)] = 1000 + r * 100 + c;
      }
      // Compound name in column Z (index 25 = 1 + 24)
      entries['Z' + (r + 1)] = `CPD-${r + 1}`;
    }
    const grid = buildGrid(entries);
    const series = extractPlateData(grid, stripLayout);

    // 2 rows × 2 strips = 4 series
    expect(series).toHaveLength(4);
    expect(series[0].stripIndex).toBe(0);
    expect(series[1].stripIndex).toBe(1);
    expect(series[2].stripIndex).toBe(0);
    expect(series[3].stripIndex).toBe(1);
  });

  it('reads compound names from adjacent_column', () => {
    const entries: Record<string, string | number> = {};
    // One row of data + compound name
    for (let c = 1; c <= 24; c++) {
      entries[String.fromCharCode(65 + c) + '1'] = c * 100;
    }
    entries['Z1'] = 'NCL-00000042';
    const grid = buildGrid(entries);

    const oneRowLayout: PlateLayout = {
      ...stripLayout,
      sample_region: { ...stripLayout.sample_region, end_row: 'A' },
    };
    const series = extractPlateData(grid, oneRowLayout);

    // Both strips should share the same compound name
    expect(series[0].compoundName).toBe('NCL-00000042');
    expect(series[1].compoundName).toBe('NCL-00000042');
  });

  it('reads correct strip cell values', () => {
    const entries: Record<string, number> = {};
    // Strip 1: B1(min=10), C1-L1(data=20..29), M1(max=30)
    entries['B1'] = 10;
    for (let i = 0; i < 10; i++) {
      entries[String.fromCharCode(67 + i) + '1'] = 20 + i; // C..L
    }
    entries['M1'] = 30;
    // Strip 2: N1(min=40), O1-X1(data=50..59), Y1(max=60)
    entries['N1'] = 40;
    for (let i = 0; i < 10; i++) {
      entries[String.fromCharCode(79 + i) + '1'] = 50 + i; // O..X
    }
    entries['Y1'] = 60;
    entries['Z1'] = 'TEST';

    const grid = buildGrid(entries);
    const oneRowLayout: PlateLayout = {
      ...stripLayout,
      sample_region: { ...stripLayout.sample_region, end_row: 'A' },
    };
    const series = extractPlateData(grid, oneRowLayout);

    // Strip 1
    expect(series[0].minControlValues).toEqual([10]);
    expect(series[0].maxControlValues).toEqual([30]);
    expect(series[0].dataValues).toEqual([10, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]);

    // Strip 2
    expect(series[1].minControlValues).toEqual([40]);
    expect(series[1].maxControlValues).toEqual([60]);
    expect(series[1].dataValues).toEqual([40, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60]);
  });

  it('flags rows with no data', () => {
    // Empty grid — all nulls
    const grid: CellGrid = Array.from({ length: 5 }, () =>
      Array.from({ length: 30 }, () => null)
    );
    const oneRowLayout: PlateLayout = {
      ...stripLayout,
      sample_region: { ...stripLayout.sample_region, end_row: 'A' },
    };
    const series = extractPlateData(grid, oneRowLayout);
    expect(series.length).toBeGreaterThan(0);
    expect(series[0].hasIssues).toBe(true);
    expect(series[0].issues).toContain('No numeric data found');
  });
});

// ---------------------------------------------------------------------------
// Integration: real HTRF plate file
// ---------------------------------------------------------------------------

describe('extractPlateData — HTRF_06032024_T1.xlsx', () => {
  // Read the real Excel file and convert section 3 (Ratio data) to a CellGrid.
  // Section 3 starts at row 47 (header), plate rows at rows 49-64.
  let cells: CellGrid;

  // Read the Excel fixture once
  const XLSX = require(
    path.resolve(__dirname, '../../../node_modules/xlsx')
  );
  const wb = XLSX.readFile(
    path.resolve(__dirname, 'fixtures/HTRF_06032024_T1.xlsx')
  );
  const ws = wb.Sheets[wb.SheetNames[0]];
  const range = XLSX.utils.decode_range(ws['!ref'] || 'A1');

  // Build absolute-coordinate grid (the same way the fixed FileDropZone does)
  cells = [];
  for (let row = 0; row <= range.e.r; row++) {
    const rowData: (string | number | null)[] = [];
    for (let col = 0; col <= range.e.c; col++) {
      const cellAddress = XLSX.utils.encode_cell({ r: row, c: col });
      const cell = ws[cellAddress];
      rowData.push(cell ? (cell.v as string | number | null) : null);
    }
    cells.push(rowData);
  }

  // Layout matching the HTRF plate:
  //   - 384-well plate, section 3 (Ratio) starts at row 49 for plate row A
  //   - Strip: [ctrl(1), 10 data, ctrl(1)] × 2 side-by-side
  //   - Compound names in column Z (adjacent_column)
  //   - spreadsheet_origin: B49 (plate row A = Excel row 49, column B)
  const htfrLayout: PlateLayout = {
    plate_format: 384,
    controls: {
      placement: 'per_compound',
      max: { columns: [], rows: [] },
      min: { columns: [], rows: [] },
    },
    sample_region: {
      start_column: 1,
      end_column: 24,
      start_row: 'A',
      end_row: 'P',
    },
    dilution: { direction: 'horizontal', num_concentrations: 10 },
    replicate: { count: 2, pattern: 'adjacent_rows' },
    compound_source: { type: 'adjacent_column' },
    strip_layout: {
      strip_width: 12,
      min_wells: 1,
      data_wells: 10,
      max_wells: 1,
      strips_per_row: 2,
    },
    spreadsheet_origin: { column: 'B', row: 49 },
  };

  it('reads all 16 plate rows × 2 strips = 32 series', () => {
    const series = extractPlateData(cells, htfrLayout);
    expect(series).toHaveLength(32); // 16 rows × 2 strips
  });

  it('first 9 rows (A-I) have no data, marked with issues', () => {
    const series = extractPlateData(cells, htfrLayout);
    // Rows A-I = indices 0-8, each has 2 strips → first 18 series
    const emptyRows = series.filter(s => s.row < 9);
    expect(emptyRows).toHaveLength(18); // 9 rows × 2 strips
    for (const s of emptyRows) {
      expect(s.hasIssues).toBe(true);
      expect(s.issues).toContain('No numeric data found');
    }
  });

  it('rows J-P (indices 9-15) have valid data', () => {
    const series = extractPlateData(cells, htfrLayout);
    const dataRows = series.filter(s => s.row >= 9 && !s.hasIssues);
    // 7 rows × 2 strips = 14 valid series
    expect(dataRows).toHaveLength(14);
  });

  it('extracts correct compound names from column Z', () => {
    const series = extractPlateData(cells, htfrLayout);
    // Row J (index 9), strip 0 should be NCL-00000001
    const rowJ = series.filter(s => s.row === 9);
    expect(rowJ).toHaveLength(2);
    expect(rowJ[0].compoundName).toBe('NCL-00000001');
    expect(rowJ[1].compoundName).toBe('NCL-00000001'); // same compound, replicate

    // Row P (index 15) should be NCL-00000007
    const rowP = series.filter(s => s.row === 15);
    expect(rowP[0].compoundName).toBe('NCL-00000007');
  });

  it('reads correct data values for row J strip 1 (NCL-00000001)', () => {
    const series = extractPlateData(cells, htfrLayout);
    const s = series.find(s => s.row === 9 && s.stripIndex === 0)!;

    // Strip 1: B58=min(5158.37), C58-L58=data, M58=max(2738.25)
    // (Note: code calls first well "min" and last "max" regardless of actual signal)
    expect(s.minControl).toBeCloseTo(5158.37, 1);
    expect(s.maxControl).toBeCloseTo(2738.25, 1);

    // 10 data points (C58-L58)
    // dataValues = [min_avg, data0..data9, max_avg]
    expect(s.dataValues).toHaveLength(12); // 1 + 10 + 1
    expect(s.dataValues[1]).toBeCloseTo(4957.53, 0);  // C58
    expect(s.dataValues[10]).toBeCloseTo(2682.81, 0); // L58
  });

  it('reads correct data values for row J strip 2 (replicate)', () => {
    const series = extractPlateData(cells, htfrLayout);
    const s = series.find(s => s.row === 9 && s.stripIndex === 1)!;

    // Strip 2: N58=min(5558.85), O58-X58=data, Y58=max(2646.11)
    expect(s.minControl).toBeCloseTo(5558.85, 1);
    expect(s.maxControl).toBeCloseTo(2646.11, 1);
    expect(s.dataValues).toHaveLength(12);
    expect(s.dataValues[1]).toBeCloseTo(4829.08, 0);  // O58
    expect(s.dataValues[10]).toBeCloseTo(2644.09, 0); // X58
  });

  it('empty rows A-I have null compound names', () => {
    const series = extractPlateData(cells, htfrLayout);
    const rowA = series.filter(s => s.row === 0);
    expect(rowA).toHaveLength(2);
    expect(rowA[0].compoundName).toBeNull();
  });
});
