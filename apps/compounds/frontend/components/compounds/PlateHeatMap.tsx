'use client';

import { useState, useEffect } from 'react';
import {
  Box,
  Typography,
  Tooltip,
  Divider,
  Dialog,
  DialogTitle,
  DialogContent,
  IconButton,
  Button,
  Skeleton,
} from '@mui/material';
import { Close, Palette } from '@mui/icons-material';
import type { PlateLayout, PlateFormat } from '@/types/compounds/models';
import { useRDKit } from '@/lib/compounds/rdkit-context';

/**
 * Compound information for a well position.
 */
export interface WellCompoundInfo {
  compoundId?: string;      // Formatted ID like "TMB-001"
  compoundName?: string;    // Human-readable name
  smiles?: string;          // SMILES string for molecule rendering
}

/**
 * Map from well position (e.g., "A1", "B12") to compound info.
 */
export type WellCompoundMapping = Record<string, WellCompoundInfo>;

const PLATE_DIMENSIONS: Record<PlateFormat, { rows: number; cols: number }> = {
  24: { rows: 4, cols: 6 },
  96: { rows: 8, cols: 12 },
  384: { rows: 16, cols: 24 },
  1536: { rows: 32, cols: 48 },
};

function excelColumnToIndex(col: string): number {
  let result = 0;
  for (let i = 0; i < col.length; i++) {
    result = result * 26 + (col.charCodeAt(i) - 'A'.charCodeAt(0) + 1);
  }
  return result - 1; // 0-indexed
}

/**
 * Small molecule thumbnail for tooltips.
 * Renders SMILES as an SVG using RDKit.
 */
function MoleculeTooltipThumb({ smiles, size = 80 }: { smiles: string; size?: number }) {
  const { rdkitModule, isLoading } = useRDKit();
  const [svgUrl, setSvgUrl] = useState<string | null>(null);

  useEffect(() => {
    if (!rdkitModule || !smiles) return;

    try {
      const mol = rdkitModule.get_mol(smiles);
      if (!mol) return;

      const svg = mol.get_svg(size, size);
      mol.delete();

      const blob = new Blob([svg], { type: 'image/svg+xml' });
      const url = URL.createObjectURL(blob);
      setSvgUrl(url);

      return () => URL.revokeObjectURL(url);
    } catch (err) {
      console.error('RDKit error in tooltip:', err);
    }
  }, [smiles, rdkitModule, size]);

  if (isLoading) {
    return <Skeleton variant="rectangular" width={size} height={size} />;
  }

  if (!svgUrl) {
    return null;
  }

  return (
    <Box
      sx={{
        width: size,
        height: size,
        bgcolor: 'white',
        borderRadius: 0.5,
        mt: 0.5,
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
      }}
    >
      <img
        src={svgUrl}
        alt=""
        style={{ width: size, height: size, objectFit: 'contain' }}
      />
    </Box>
  );
}

interface PlateHeatMapProps {
  cells: (string | number | null)[][];
  plateLayout: PlateLayout;
  /** Optional mapping from well positions (e.g., "A1") to compound info */
  compoundMapping?: WellCompoundMapping;
}

export function PlateHeatMap({ cells, plateLayout, compoundMapping }: PlateHeatMapProps) {
  const format = plateLayout.plate_format || 384;
  const { rows: numRows, cols: numCols } = PLATE_DIMENSIONS[format];

  // Get origin offset
  const originCol = excelColumnToIndex(plateLayout.spreadsheet_origin?.column || 'A');
  const originRow = (plateLayout.spreadsheet_origin?.row || 1) - 1;

  // Extract plate data and find min/max for color scaling
  const plateData: (number | null)[][] = [];
  let minVal = Infinity;
  let maxVal = -Infinity;

  for (let row = 0; row < numRows; row++) {
    const rowData: (number | null)[] = [];
    const spreadsheetRow = originRow + row;

    for (let col = 0; col < numCols; col++) {
      const spreadsheetCol = originCol + col;
      let value: number | null = null;

      if (spreadsheetRow < cells.length) {
        const cellRow = cells[spreadsheetRow];
        if (cellRow && spreadsheetCol < cellRow.length) {
          const cellVal = cellRow[spreadsheetCol];
          if (typeof cellVal === 'number') {
            value = cellVal;
            if (value < minVal) minVal = value;
            if (value > maxVal) maxVal = value;
          }
        }
      }
      rowData.push(value);
    }
    plateData.push(rowData);
  }

  // If no numeric data found, set defaults
  if (minVal === Infinity) minVal = 0;
  if (maxVal === -Infinity) maxVal = 100;

  // Color interpolation function (blue -> white -> red)
  const getColor = (value: number | null): string => {
    if (value === null) return '#f5f5f5';

    const range = maxVal - minVal;
    if (range === 0) return '#ffffff';

    const normalized = (value - minVal) / range; // 0 to 1

    // Blue (low) -> White (mid) -> Red (high)
    if (normalized < 0.5) {
      // Blue to white
      const t = normalized * 2;
      const r = Math.round(59 + (255 - 59) * t);
      const g = Math.round(130 + (255 - 130) * t);
      const b = Math.round(246 + (255 - 246) * t);
      return `rgb(${r}, ${g}, ${b})`;
    } else {
      // White to red
      const t = (normalized - 0.5) * 2;
      const r = 255;
      const g = Math.round(255 - (255 - 82) * t);
      const b = Math.round(255 - (255 - 82) * t);
      return `rgb(${r}, ${g}, ${b})`;
    }
  };

  // Determine cell type (control, data, etc.) for border styling
  const getCellType = (row: number, col: number): 'max' | 'min' | 'data' | 'outside' => {
    const col1Indexed = col + 1;
    const rowLetter = String.fromCharCode('A'.charCodeAt(0) + row);

    if (plateLayout.controls?.placement === 'per_compound' && plateLayout.strip_layout) {
      // Strip layout: controls are embedded in strips
      const strip = plateLayout.strip_layout;
      const posInRow = col % strip.strip_width;

      if (posInRow < strip.min_wells) return 'min';
      if (posInRow >= strip.min_wells + strip.data_wells) return 'max';
      return 'data';
    }

    // Edge controls
    if (plateLayout.controls?.max?.columns?.includes(col1Indexed) &&
        plateLayout.controls?.max?.rows?.includes(rowLetter)) {
      return 'max';
    }
    if (plateLayout.controls?.min?.columns?.includes(col1Indexed) &&
        plateLayout.controls?.min?.rows?.includes(rowLetter)) {
      return 'min';
    }

    // Check if in sample region
    const sampleRegion = plateLayout.sample_region;
    if (sampleRegion) {
      const startCol = sampleRegion.start_column;
      const endCol = sampleRegion.end_column;
      const startRowIdx = sampleRegion.start_row.charCodeAt(0) - 'A'.charCodeAt(0);
      const endRowIdx = sampleRegion.end_row.charCodeAt(0) - 'A'.charCodeAt(0);

      if (col1Indexed >= startCol && col1Indexed <= endCol &&
          row >= startRowIdx && row <= endRowIdx) {
        return 'data';
      }
    }

    return 'outside';
  };

  // Calculate cell size based on plate format
  const cellSize = format === 1536 ? 12 : format === 384 ? 18 : format === 96 ? 28 : 40;

  return (
    <Box>
      {/* Legend */}
      <Box sx={{ display: 'flex', gap: 3, mb: 2, alignItems: 'center', flexWrap: 'wrap' }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Box sx={{ width: 16, height: 16, bgcolor: '#3b82f6', borderRadius: 0.5 }} />
          <Typography variant="caption">Low ({minVal.toFixed(0)})</Typography>
        </Box>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Box sx={{ width: 16, height: 16, bgcolor: '#ffffff', border: '1px solid #ccc', borderRadius: 0.5 }} />
          <Typography variant="caption">Mid</Typography>
        </Box>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Box sx={{ width: 16, height: 16, bgcolor: '#ef5252', borderRadius: 0.5 }} />
          <Typography variant="caption">High ({maxVal.toFixed(0)})</Typography>
        </Box>
        <Divider orientation="vertical" flexItem />
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Box sx={{ width: 16, height: 16, border: '2px solid #1976d2', borderRadius: 0.5 }} />
          <Typography variant="caption">Max Control</Typography>
        </Box>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Box sx={{ width: 16, height: 16, border: '2px solid #d32f2f', borderRadius: 0.5 }} />
          <Typography variant="caption">Min Control</Typography>
        </Box>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Box sx={{ width: 16, height: 16, bgcolor: '#f5f5f5', borderRadius: 0.5 }} />
          <Typography variant="caption">No Data</Typography>
        </Box>
      </Box>

      {/* Plate grid */}
      <Box sx={{ overflowX: 'auto' }}>
        <Box sx={{ display: 'inline-block', minWidth: 'fit-content' }}>
          {/* Column headers */}
          <Box sx={{ display: 'flex', ml: `${cellSize + 4}px` }}>
            {Array.from({ length: numCols }, (_, col) => (
              <Box
                key={col}
                sx={{
                  width: cellSize,
                  height: 20,
                  display: 'flex',
                  alignItems: 'center',
                  justifyContent: 'center',
                  fontSize: format > 96 ? 8 : 10,
                  color: 'text.secondary',
                }}
              >
                {col + 1}
              </Box>
            ))}
          </Box>

          {/* Rows */}
          {plateData.map((rowData, row) => (
            <Box key={row} sx={{ display: 'flex', alignItems: 'center' }}>
              {/* Row label */}
              <Box
                sx={{
                  width: cellSize,
                  height: cellSize,
                  display: 'flex',
                  alignItems: 'center',
                  justifyContent: 'center',
                  fontSize: format > 96 ? 9 : 11,
                  fontWeight: 500,
                  color: 'text.secondary',
                  mr: 0.5,
                }}
              >
                {String.fromCharCode('A'.charCodeAt(0) + row)}
              </Box>

              {/* Cells */}
              {rowData.map((value, col) => {
                const cellType = getCellType(row, col);
                const wellPosition = `${String.fromCharCode('A'.charCodeAt(0) + row)}${col + 1}`;
                const compoundInfo = compoundMapping?.[wellPosition];

                return (
                  <Tooltip
                    key={col}
                    title={
                      <Box sx={{ minWidth: compoundInfo?.smiles ? 100 : 'auto' }}>
                        <Typography variant="caption" display="block" fontWeight="bold">
                          {wellPosition}
                        </Typography>
                        <Typography variant="caption" display="block">
                          Value: {value !== null ? value.toFixed(2) : 'N/A'}
                        </Typography>
                        <Typography variant="caption" display="block" sx={{ textTransform: 'capitalize' }}>
                          Type: {cellType}
                        </Typography>
                        {compoundInfo && (
                          <>
                            <Divider sx={{ my: 0.5, borderColor: 'rgba(255,255,255,0.3)' }} />
                            {compoundInfo.compoundName && (
                              <Typography variant="caption" display="block">
                                Name: {compoundInfo.compoundName}
                              </Typography>
                            )}
                            {compoundInfo.compoundId && (
                              <Typography variant="caption" display="block" fontFamily="monospace">
                                ID: {compoundInfo.compoundId}
                              </Typography>
                            )}
                            {compoundInfo.smiles && (
                              <MoleculeTooltipThumb smiles={compoundInfo.smiles} size={80} />
                            )}
                          </>
                        )}
                      </Box>
                    }
                    arrow
                    placement="top"
                  >
                    <Box
                      sx={{
                        width: cellSize,
                        height: cellSize,
                        bgcolor: getColor(value),
                        border: cellType === 'max'
                          ? '2px solid #1976d2'
                          : cellType === 'min'
                          ? '2px solid #d32f2f'
                          : '1px solid #e0e0e0',
                        borderRadius: 0.5,
                        cursor: 'pointer',
                        transition: 'transform 0.1s',
                        '&:hover': {
                          transform: 'scale(1.1)',
                          zIndex: 1,
                        },
                      }}
                    />
                  </Tooltip>
                );
              })}
            </Box>
          ))}
        </Box>
      </Box>
    </Box>
  );
}

/**
 * Heat map dialog component that can be used standalone
 */
interface PlateHeatMapDialogProps {
  open: boolean;
  onClose: () => void;
  cells: (string | number | null)[][];
  plateLayout: PlateLayout;
  /** Optional mapping from well positions (e.g., "A1") to compound info */
  compoundMapping?: WellCompoundMapping;
}

export function PlateHeatMapDialog({ open, onClose, cells, plateLayout, compoundMapping }: PlateHeatMapDialogProps) {
  return (
    <Dialog open={open} onClose={onClose} maxWidth="lg" fullWidth>
      <DialogTitle sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Palette color="primary" />
          <Typography variant="h6">Plate Heat Map</Typography>
        </Box>
        <IconButton onClick={onClose} size="small">
          <Close />
        </IconButton>
      </DialogTitle>
      <DialogContent>
        <PlateHeatMap cells={cells} plateLayout={plateLayout} compoundMapping={compoundMapping} />
      </DialogContent>
    </Dialog>
  );
}

/**
 * Button with integrated dialog for showing heat map
 */
interface PlateHeatMapButtonProps {
  cells: (string | number | null)[][];
  plateLayout: PlateLayout;
  /** Optional mapping from well positions (e.g., "A1") to compound info */
  compoundMapping?: WellCompoundMapping;
  variant?: 'text' | 'outlined' | 'contained';
  size?: 'small' | 'medium' | 'large';
}

export function PlateHeatMapButton({ cells, plateLayout, compoundMapping, variant = 'outlined', size = 'small' }: PlateHeatMapButtonProps) {
  const [open, setOpen] = useState(false);

  return (
    <>
      <Button
        size={size}
        variant={variant}
        startIcon={<Palette />}
        onClick={() => setOpen(true)}
      >
        Heat Map
      </Button>
      <PlateHeatMapDialog
        open={open}
        onClose={() => setOpen(false)}
        cells={cells}
        plateLayout={plateLayout}
        compoundMapping={compoundMapping}
      />
    </>
  );
}
