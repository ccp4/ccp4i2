'use client';

import { useState, useCallback } from 'react';
import {
  Box,
  Stepper,
  Step,
  StepLabel,
  StepContent,
  Button,
  Typography,
  FormControl,
  FormControlLabel,
  FormLabel,
  RadioGroup,
  Radio,
  Select,
  MenuItem,
  InputLabel,
  TextField,
  Chip,
  ToggleButton,
  ToggleButtonGroup,
  Grid2 as Grid,
  Paper,
  Alert,
} from '@mui/material';
import {
  GridOn,
  Science,
  LinearScale,
  ContentCopy,
  Source,
  TableChart,
} from '@mui/icons-material';
import { PlatePreview } from './PlatePreview';
import type {
  PlateLayout,
  PlateFormat,
  ReplicatePattern,
  CompoundSourceType,
  DilutionDirection,
  ControlPlacement,
  SpreadsheetOrigin,
} from '@/types/compounds/models';

/**
 * Plate dimensions for generating row/column options
 */
const PLATE_DIMENSIONS: Record<PlateFormat, { rows: number; cols: number }> = {
  24: { rows: 4, cols: 6 },
  96: { rows: 8, cols: 12 },
  384: { rows: 16, cols: 24 },
  1536: { rows: 32, cols: 48 },
};

function indexToRowLetter(index: number): string {
  return String.fromCharCode('A'.charCodeAt(0) + index);
}

function getRowLetters(numRows: number): string[] {
  return Array.from({ length: numRows }, (_, i) => indexToRowLetter(i));
}

function getColumnNumbers(numCols: number): number[] {
  return Array.from({ length: numCols }, (_, i) => i + 1);
}

/**
 * Default plate layout for a 384-well plate
 */
const DEFAULT_LAYOUT: PlateLayout = {
  plate_format: 384,
  controls: {
    placement: 'edge_columns',
    max: { columns: [1, 2], rows: ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P'] },
    min: { columns: [23, 24], rows: ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P'] },
  },
  sample_region: {
    start_column: 3,
    end_column: 22,
    start_row: 'A',
    end_row: 'P',
  },
  dilution: {
    direction: 'horizontal',
    num_concentrations: 20,
  },
  replicate: {
    count: 2,
    pattern: 'adjacent_rows',
  },
  compound_source: {
    type: 'row_order',
  },
  spreadsheet_origin: {
    column: 'A',
    row: 1,
  },
};

/**
 * Common plate layout templates
 */
const TEMPLATES: { name: string; description: string; layout: PlateLayout }[] = [
  {
    name: 'Standard 384-well (20-point)',
    description: '2 high + 2 low edge control columns, 20-point dilution series, 2 replicates',
    layout: DEFAULT_LAYOUT,
  },
  {
    name: 'Standard 96-well (8-point)',
    description: 'Edge controls, 8 concentrations horizontal, no replicates',
    layout: {
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
        end_row: 'H',
      },
      dilution: {
        direction: 'horizontal',
        num_concentrations: 8,
      },
      replicate: {
        count: 1,
        pattern: 'adjacent_rows',
      },
      compound_source: {
        type: 'row_order',
      },
    },
  },
  {
    name: '384-well Vertical Dilution',
    description: 'Top/bottom controls, vertical dilution, adjacent column replicates',
    layout: {
      plate_format: 384,
      controls: {
        placement: 'edge_rows',
        max: { columns: Array.from({ length: 24 }, (_, i) => i + 1), rows: ['A', 'B'] },
        min: { columns: Array.from({ length: 24 }, (_, i) => i + 1), rows: ['O', 'P'] },
      },
      sample_region: {
        start_column: 1,
        end_column: 24,
        start_row: 'C',
        end_row: 'N',
      },
      dilution: {
        direction: 'vertical',
        num_concentrations: 12,
      },
      replicate: {
        count: 2,
        pattern: 'adjacent_columns',
      },
      compound_source: {
        type: 'column_header',
      },
    },
  },
  {
    name: '384-well Strip (10-point, 2 reps)',
    description: 'Embedded controls per strip: [min×1][data×10][max×1] × 2 per row',
    layout: {
      plate_format: 384,
      controls: {
        placement: 'per_compound',
        max: { columns: [], rows: [] },  // Not used for strip layout
        min: { columns: [], rows: [] },
      },
      sample_region: {
        start_column: 1,
        end_column: 24,
        start_row: 'A',
        end_row: 'P',
      },
      dilution: {
        direction: 'horizontal',
        num_concentrations: 10,
      },
      replicate: {
        count: 2,
        pattern: 'adjacent_columns',  // Strips are adjacent along row
      },
      compound_source: {
        type: 'row_order',
      },
      strip_layout: {
        strip_width: 12,
        min_wells: 1,
        data_wells: 10,
        max_wells: 1,
        strips_per_row: 2,
      },
    },
  },
  {
    name: '1536-well Strip (10-point, 4 strips paired)',
    description: 'Embedded controls: [min×1][data×10][max×1] × 4 strips, paired replicates (1+2, 3+4)',
    layout: {
      plate_format: 1536,
      controls: {
        placement: 'per_compound',
        max: { columns: [], rows: [] },
        min: { columns: [], rows: [] },
      },
      sample_region: {
        start_column: 1,
        end_column: 48,
        start_row: 'A',
        end_row: 'P',  // First 16 rows; plate has 32 total
      },
      dilution: {
        direction: 'horizontal',
        num_concentrations: 10,
      },
      replicate: {
        count: 2,
        pattern: 'adjacent_columns',  // Strips 1+2 paired, 3+4 paired
      },
      compound_source: {
        type: 'row_order',
      },
      strip_layout: {
        strip_width: 12,
        min_wells: 1,
        data_wells: 10,
        max_wells: 1,
        strips_per_row: 4,
      },
    },
  },
];

const STEPS = [
  { label: 'Plate Format', icon: <GridOn /> },
  { label: 'Controls', icon: <Science /> },
  { label: 'Sample Region', icon: <LinearScale /> },
  { label: 'Replicates', icon: <ContentCopy /> },
  { label: 'Compound Source', icon: <Source /> },
  { label: 'Data Import', icon: <TableChart /> },
];

interface PlateLayoutEditorProps {
  value: Partial<PlateLayout>;
  onChange: (layout: PlateLayout) => void;
  showPreview?: boolean;
}

export function PlateLayoutEditor({
  value,
  onChange,
  showPreview = true,
}: PlateLayoutEditorProps) {
  const [activeStep, setActiveStep] = useState(0);
  const [layout, setLayout] = useState<PlateLayout>({
    ...DEFAULT_LAYOUT,
    ...value,
  } as PlateLayout);
  // Track the raw compound name row input (can be "+2" relative or "25" absolute)
  const [compoundNameRowInput, setCompoundNameRowInput] = useState<string>('');

  const plateFormat = layout.plate_format;
  const { rows: numRows, cols: numCols } = PLATE_DIMENSIONS[plateFormat];
  const rowLetters = getRowLetters(numRows);
  const columnNumbers = getColumnNumbers(numCols);

  const updateLayout = useCallback((updates: Partial<PlateLayout>) => {
    const newLayout = { ...layout, ...updates };
    setLayout(newLayout);
    onChange(newLayout);
  }, [layout, onChange]);

  // Calculate the Excel row at the bottom of the plate (for relative compound name row calculation)
  // Use the full plate height, not just the sample region, so names don't overlap with plate data
  const getPlateBottomRow = useCallback((): number => {
    const originRow = layout.spreadsheet_origin?.row || 1;
    const { rows: plateRows } = PLATE_DIMENSIONS[layout.plate_format];
    return originRow + plateRows - 1;  // 1-indexed Excel row of last plate row
  }, [layout.spreadsheet_origin?.row, layout.plate_format]);

  // Parse compound name row input: "+N" means relative to data bottom, "N" means absolute
  const handleCompoundNameRowChange = useCallback((inputValue: string) => {
    setCompoundNameRowInput(inputValue);

    if (!inputValue.trim()) {
      // Clear the setting
      updateLayout({
        compound_source: { ...layout.compound_source, compound_name_row: undefined },
      });
      return;
    }

    const trimmed = inputValue.trim();
    let absoluteRow: number;

    if (trimmed.startsWith('+')) {
      // Relative: add offset to data bottom row
      const offset = parseInt(trimmed.slice(1), 10);
      if (!isNaN(offset)) {
        absoluteRow = getPlateBottomRow() + offset;
      } else {
        return; // Invalid input, don't update
      }
    } else {
      // Absolute row number
      absoluteRow = parseInt(trimmed, 10);
      if (isNaN(absoluteRow) || absoluteRow < 1) {
        return; // Invalid input, don't update
      }
    }

    updateLayout({
      compound_source: { ...layout.compound_source, compound_name_row: absoluteRow },
    });
  }, [layout.compound_source, getPlateBottomRow, updateLayout]);

  const handlePlateFormatChange = (format: PlateFormat) => {
    // When format changes, reset regions to sensible defaults
    const { rows: newRows, cols: newCols } = PLATE_DIMENSIONS[format];
    const newRowLetters = getRowLetters(newRows);

    updateLayout({
      plate_format: format,
      controls: {
        placement: 'edge_columns',
        max: { columns: [1, 2], rows: newRowLetters },
        min: { columns: [newCols - 1, newCols], rows: newRowLetters },
      },
      sample_region: {
        start_column: 3,
        end_column: newCols - 2,
        start_row: 'A',
        end_row: indexToRowLetter(newRows - 1),
      },
    });
  };

  const handleTemplateSelect = (template: PlateLayout) => {
    setLayout(template);
    onChange(template);
    setActiveStep(0);
  };

  const handleNext = () => {
    setActiveStep((prev) => Math.min(prev + 1, STEPS.length - 1));
  };

  const handleBack = () => {
    setActiveStep((prev) => Math.max(prev - 1, 0));
  };

  const toggleControlColumn = (controlType: 'max' | 'min', col: number) => {
    const current = layout.controls[controlType].columns;
    const updated = current.includes(col)
      ? current.filter(c => c !== col)
      : [...current, col].sort((a, b) => a - b);

    updateLayout({
      controls: {
        ...layout.controls,
        [controlType]: { ...layout.controls[controlType], columns: updated },
      },
    });
  };

  const toggleControlRow = (controlType: 'max' | 'min', row: string) => {
    const current = layout.controls[controlType].rows;
    const updated = current.includes(row)
      ? current.filter(r => r !== row)
      : [...current, row].sort();

    updateLayout({
      controls: {
        ...layout.controls,
        [controlType]: { ...layout.controls[controlType], rows: updated },
      },
    });
  };

  return (
    <Grid container spacing={3}>
      {/* Editor panel */}
      <Grid size={{ xs: 12, md: showPreview ? 7 : 12 }}>
        {/* Template selector */}
        <Paper sx={{ p: 2, mb: 2 }}>
          <Typography variant="subtitle2" gutterBottom>
            Start from template
          </Typography>
          <Box sx={{ display: 'flex', gap: 1, flexWrap: 'wrap' }}>
            {TEMPLATES.map((template, idx) => (
              <Chip
                key={idx}
                label={template.name}
                onClick={() => handleTemplateSelect(template.layout)}
                variant={
                  JSON.stringify(layout) === JSON.stringify(template.layout)
                    ? 'filled'
                    : 'outlined'
                }
                color="primary"
                size="small"
              />
            ))}
          </Box>
        </Paper>

        {/* Stepper - nonLinear allows clicking on any step to edit */}
        <Stepper activeStep={activeStep} orientation="vertical" nonLinear>
          {/* Step 1: Plate Format */}
          <Step completed={activeStep > 0}>
            <StepLabel
              onClick={() => setActiveStep(0)}
              sx={{ cursor: 'pointer', '&:hover': { bgcolor: 'action.hover' }, borderRadius: 1, py: 0.5 }}
            >
              {STEPS[0].label}
            </StepLabel>
            <StepContent>
              <FormControl component="fieldset">
                <FormLabel>Select plate format</FormLabel>
                <RadioGroup
                  row
                  value={plateFormat}
                  onChange={(e) => handlePlateFormatChange(Number(e.target.value) as PlateFormat)}
                >
                  <FormControlLabel value={24} control={<Radio />} label="24-well" />
                  <FormControlLabel value={96} control={<Radio />} label="96-well" />
                  <FormControlLabel value={384} control={<Radio />} label="384-well" />
                  <FormControlLabel value={1536} control={<Radio />} label="1536-well" />
                </RadioGroup>
              </FormControl>
              <Box sx={{ mt: 2 }}>
                <Button variant="contained" onClick={handleNext}>
                  Continue
                </Button>
              </Box>
            </StepContent>
          </Step>

          {/* Step 2: Controls */}
          <Step completed={activeStep > 1}>
            <StepLabel
              onClick={() => setActiveStep(1)}
              sx={{ cursor: 'pointer', '&:hover': { bgcolor: 'action.hover' }, borderRadius: 1, py: 0.5 }}
            >
              {STEPS[1].label}
            </StepLabel>
            <StepContent>
              {/* Control placement type */}
              <FormControl component="fieldset" sx={{ mb: 3 }}>
                <FormLabel>Control Placement</FormLabel>
                <RadioGroup
                  row
                  value={layout.controls.placement || 'edge_columns'}
                  onChange={(e) => {
                    const placement = e.target.value as ControlPlacement;
                    if (placement === 'per_compound') {
                      // Switch to strip layout
                      updateLayout({
                        controls: { ...layout.controls, placement },
                        strip_layout: layout.strip_layout || {
                          strip_width: 12,
                          min_wells: 2,
                          data_wells: 8,
                          max_wells: 2,
                          strips_per_row: 2,
                        },
                      });
                    } else {
                      updateLayout({
                        controls: { ...layout.controls, placement },
                        strip_layout: undefined,
                      });
                    }
                  }}
                >
                  <FormControlLabel
                    value="edge_columns"
                    control={<Radio size="small" />}
                    label="Edge Columns"
                  />
                  <FormControlLabel
                    value="edge_rows"
                    control={<Radio size="small" />}
                    label="Edge Rows"
                  />
                  <FormControlLabel
                    value="per_compound"
                    control={<Radio size="small" />}
                    label="Per Compound Strip"
                  />
                </RadioGroup>
              </FormControl>

              {layout.controls.placement === 'per_compound' ? (
                /* Strip layout configuration */
                <Box sx={{ mb: 2 }}>
                  <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
                    Each compound has its own control wells in a repeating strip pattern:
                    <br />
                    <strong>[min×N] [data×M] [max×N]</strong> repeated across the row
                  </Typography>

                  <Grid container spacing={2}>
                    <Grid size={{ xs: 6, sm: 3 }}>
                      <TextField
                        label="Min Wells"
                        type="number"
                        size="small"
                        fullWidth
                        value={layout.strip_layout?.min_wells || 2}
                        onChange={(e) => updateLayout({
                          strip_layout: {
                            ...layout.strip_layout!,
                            min_wells: Number(e.target.value),
                            strip_width: Number(e.target.value) + (layout.strip_layout?.data_wells || 8) + (layout.strip_layout?.max_wells || 2),
                          },
                        })}
                        inputProps={{ min: 1, max: 4 }}
                      />
                    </Grid>
                    <Grid size={{ xs: 6, sm: 3 }}>
                      <TextField
                        label="Data Wells"
                        type="number"
                        size="small"
                        fullWidth
                        value={layout.strip_layout?.data_wells || 8}
                        onChange={(e) => updateLayout({
                          strip_layout: {
                            ...layout.strip_layout!,
                            data_wells: Number(e.target.value),
                            strip_width: (layout.strip_layout?.min_wells || 2) + Number(e.target.value) + (layout.strip_layout?.max_wells || 2),
                          },
                          dilution: {
                            ...layout.dilution,
                            num_concentrations: Number(e.target.value),
                          },
                        })}
                        inputProps={{ min: 4, max: 20 }}
                      />
                    </Grid>
                    <Grid size={{ xs: 6, sm: 3 }}>
                      <TextField
                        label="Max Wells"
                        type="number"
                        size="small"
                        fullWidth
                        value={layout.strip_layout?.max_wells || 2}
                        onChange={(e) => updateLayout({
                          strip_layout: {
                            ...layout.strip_layout!,
                            max_wells: Number(e.target.value),
                            strip_width: (layout.strip_layout?.min_wells || 2) + (layout.strip_layout?.data_wells || 8) + Number(e.target.value),
                          },
                        })}
                        inputProps={{ min: 1, max: 4 }}
                      />
                    </Grid>
                    <Grid size={{ xs: 6, sm: 3 }}>
                      <TextField
                        label="Strips/Row"
                        type="number"
                        size="small"
                        fullWidth
                        value={layout.strip_layout?.strips_per_row || 2}
                        onChange={(e) => updateLayout({
                          strip_layout: {
                            ...layout.strip_layout!,
                            strips_per_row: Number(e.target.value),
                          },
                          replicate: {
                            ...layout.replicate,
                            count: Number(e.target.value),
                          },
                        })}
                        inputProps={{ min: 1, max: 4 }}
                      />
                    </Grid>
                  </Grid>

                  <Alert severity="info" sx={{ mt: 2 }}>
                    Strip width: {layout.strip_layout?.strip_width || 12} columns
                    ({layout.strip_layout?.min_wells || 2} min + {layout.strip_layout?.data_wells || 8} data + {layout.strip_layout?.max_wells || 2} max)
                    × {layout.strip_layout?.strips_per_row || 2} replicates per row
                  </Alert>
                </Box>
              ) : (
                /* Standard edge controls configuration */
                <>
                  {/* Max controls */}
                  <Box sx={{ mb: 3 }}>
                    <Typography variant="subtitle2" sx={{ color: 'primary.main', mb: 1 }}>
                      Max Controls (100% signal)
                    </Typography>
                    <Box sx={{ mb: 1 }}>
                      <Typography variant="caption" color="text.secondary">Columns:</Typography>
                      <Box sx={{ display: 'flex', gap: 0.5, flexWrap: 'wrap', mt: 0.5 }}>
                        {columnNumbers.slice(0, Math.min(24, numCols)).map((col) => (
                          <Chip
                            key={col}
                            label={col}
                            size="small"
                            variant={layout.controls.max.columns.includes(col) ? 'filled' : 'outlined'}
                            color="primary"
                            onClick={() => toggleControlColumn('max', col)}
                            sx={{ minWidth: 32 }}
                          />
                        ))}
                      </Box>
                    </Box>
                    <Box>
                      <Typography variant="caption" color="text.secondary">Rows:</Typography>
                      <Box sx={{ display: 'flex', gap: 0.5, flexWrap: 'wrap', mt: 0.5 }}>
                        <Chip
                          label="All"
                          size="small"
                          variant={layout.controls.max.rows.length === numRows ? 'filled' : 'outlined'}
                          onClick={() => updateLayout({
                            controls: {
                              ...layout.controls,
                              max: { ...layout.controls.max, rows: rowLetters },
                            },
                          })}
                        />
                        {rowLetters.slice(0, 8).map((row) => (
                          <Chip
                            key={row}
                            label={row}
                            size="small"
                            variant={layout.controls.max.rows.includes(row) ? 'filled' : 'outlined'}
                            color="primary"
                            onClick={() => toggleControlRow('max', row)}
                          />
                        ))}
                      </Box>
                    </Box>
                  </Box>

                  {/* Min controls */}
                  <Box sx={{ mb: 2 }}>
                    <Typography variant="subtitle2" sx={{ color: 'error.main', mb: 1 }}>
                      Min Controls (0% signal)
                    </Typography>
                    <Box sx={{ mb: 1 }}>
                      <Typography variant="caption" color="text.secondary">Columns:</Typography>
                      <Box sx={{ display: 'flex', gap: 0.5, flexWrap: 'wrap', mt: 0.5 }}>
                        {columnNumbers.slice(0, Math.min(24, numCols)).map((col) => (
                          <Chip
                            key={col}
                            label={col}
                            size="small"
                            variant={layout.controls.min.columns.includes(col) ? 'filled' : 'outlined'}
                            color="error"
                            onClick={() => toggleControlColumn('min', col)}
                            sx={{ minWidth: 32 }}
                          />
                        ))}
                      </Box>
                    </Box>
                    <Box>
                      <Typography variant="caption" color="text.secondary">Rows:</Typography>
                      <Box sx={{ display: 'flex', gap: 0.5, flexWrap: 'wrap', mt: 0.5 }}>
                        <Chip
                          label="All"
                          size="small"
                          variant={layout.controls.min.rows.length === numRows ? 'filled' : 'outlined'}
                          onClick={() => updateLayout({
                            controls: {
                              ...layout.controls,
                              min: { ...layout.controls.min, rows: rowLetters },
                            },
                          })}
                        />
                        {rowLetters.slice(0, 8).map((row) => (
                          <Chip
                            key={row}
                            label={row}
                            size="small"
                            variant={layout.controls.min.rows.includes(row) ? 'filled' : 'outlined'}
                            color="error"
                            onClick={() => toggleControlRow('min', row)}
                          />
                        ))}
                      </Box>
                    </Box>
                  </Box>
                </>
              )}

              <Box sx={{ mt: 2 }}>
                <Button onClick={handleBack} sx={{ mr: 1 }}>Back</Button>
                <Button variant="contained" onClick={handleNext}>Continue</Button>
              </Box>
            </StepContent>
          </Step>

          {/* Step 3: Sample Region & Dilution */}
          <Step completed={activeStep > 2}>
            <StepLabel
              onClick={() => setActiveStep(2)}
              sx={{ cursor: 'pointer', '&:hover': { bgcolor: 'action.hover' }, borderRadius: 1, py: 0.5 }}
            >
              {STEPS[2].label}
            </StepLabel>
            <StepContent>
              {layout.controls.placement === 'per_compound' ? (
                /* Strip layout: columns and concentrations are defined by strip config */
                <Box>
                  <Alert severity="info" sx={{ mb: 2 }}>
                    For strip layouts, these values are defined by the strip configuration in Step 2:
                    <ul style={{ margin: '8px 0 0 0', paddingLeft: 20 }}>
                      <li><strong>{layout.strip_layout?.data_wells || 8} concentrations</strong> per compound</li>
                      <li><strong>Horizontal</strong> dilution direction (across columns within each strip)</li>
                    </ul>
                  </Alert>

                  <Typography variant="subtitle2" gutterBottom>
                    Active Rows
                  </Typography>
                  <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
                    Select which rows contain compound data (e.g., exclude rows used for other purposes).
                  </Typography>

                  <Grid container spacing={2}>
                    <Grid size={{ xs: 6 }}>
                      <FormControl fullWidth size="small">
                        <InputLabel>Start Row</InputLabel>
                        <Select
                          value={layout.sample_region.start_row}
                          label="Start Row"
                          onChange={(e) => updateLayout({
                            sample_region: { ...layout.sample_region, start_row: e.target.value },
                          })}
                        >
                          {rowLetters.map((row) => (
                            <MenuItem key={row} value={row}>{row}</MenuItem>
                          ))}
                        </Select>
                      </FormControl>
                    </Grid>
                    <Grid size={{ xs: 6 }}>
                      <FormControl fullWidth size="small">
                        <InputLabel>End Row</InputLabel>
                        <Select
                          value={layout.sample_region.end_row}
                          label="End Row"
                          onChange={(e) => updateLayout({
                            sample_region: { ...layout.sample_region, end_row: e.target.value },
                          })}
                        >
                          {rowLetters.map((row) => (
                            <MenuItem key={row} value={row}>{row}</MenuItem>
                          ))}
                        </Select>
                      </FormControl>
                    </Grid>
                  </Grid>
                </Box>
              ) : (
                /* Standard layout: full sample region configuration */
                <>
                  <Grid container spacing={2}>
                    <Grid size={{ xs: 6 }}>
                      <FormControl fullWidth size="small">
                        <InputLabel>Start Column</InputLabel>
                        <Select
                          value={layout.sample_region.start_column}
                          label="Start Column"
                          onChange={(e) => updateLayout({
                            sample_region: { ...layout.sample_region, start_column: Number(e.target.value) },
                          })}
                        >
                          {columnNumbers.map((col) => (
                            <MenuItem key={col} value={col}>{col}</MenuItem>
                          ))}
                        </Select>
                      </FormControl>
                    </Grid>
                    <Grid size={{ xs: 6 }}>
                      <FormControl fullWidth size="small">
                        <InputLabel>End Column</InputLabel>
                        <Select
                          value={layout.sample_region.end_column}
                          label="End Column"
                          onChange={(e) => updateLayout({
                            sample_region: { ...layout.sample_region, end_column: Number(e.target.value) },
                          })}
                        >
                          {columnNumbers.map((col) => (
                            <MenuItem key={col} value={col}>{col}</MenuItem>
                          ))}
                        </Select>
                      </FormControl>
                    </Grid>
                    <Grid size={{ xs: 6 }}>
                      <FormControl fullWidth size="small">
                        <InputLabel>Start Row</InputLabel>
                        <Select
                          value={layout.sample_region.start_row}
                          label="Start Row"
                          onChange={(e) => updateLayout({
                            sample_region: { ...layout.sample_region, start_row: e.target.value },
                          })}
                        >
                          {rowLetters.map((row) => (
                            <MenuItem key={row} value={row}>{row}</MenuItem>
                          ))}
                        </Select>
                      </FormControl>
                    </Grid>
                    <Grid size={{ xs: 6 }}>
                      <FormControl fullWidth size="small">
                        <InputLabel>End Row</InputLabel>
                        <Select
                          value={layout.sample_region.end_row}
                          label="End Row"
                          onChange={(e) => updateLayout({
                            sample_region: { ...layout.sample_region, end_row: e.target.value },
                          })}
                        >
                          {rowLetters.map((row) => (
                            <MenuItem key={row} value={row}>{row}</MenuItem>
                          ))}
                        </Select>
                      </FormControl>
                    </Grid>
                  </Grid>

                  <Box sx={{ mt: 3 }}>
                    <Typography variant="subtitle2" gutterBottom>Dilution Direction</Typography>
                    <ToggleButtonGroup
                      value={layout.dilution.direction}
                      exclusive
                      onChange={(_, value) => value && updateLayout({
                        dilution: { ...layout.dilution, direction: value as DilutionDirection },
                      })}
                      size="small"
                    >
                      <ToggleButton value="horizontal">
                        Horizontal (across columns)
                      </ToggleButton>
                      <ToggleButton value="vertical">
                        Vertical (down rows)
                      </ToggleButton>
                    </ToggleButtonGroup>
                  </Box>

                  <Box sx={{ mt: 2 }}>
                    <TextField
                      label="Number of Concentrations"
                      type="number"
                      size="small"
                      value={layout.dilution.num_concentrations}
                      onChange={(e) => updateLayout({
                        dilution: { ...layout.dilution, num_concentrations: Number(e.target.value) },
                      })}
                      inputProps={{ min: 1, max: 24 }}
                      sx={{ width: 200 }}
                    />
                  </Box>
                </>
              )}

              <Box sx={{ mt: 2 }}>
                <Button onClick={handleBack} sx={{ mr: 1 }}>Back</Button>
                <Button variant="contained" onClick={handleNext}>Continue</Button>
              </Box>
            </StepContent>
          </Step>

          {/* Step 4: Replicates */}
          <Step completed={activeStep > 3}>
            <StepLabel
              onClick={() => setActiveStep(3)}
              sx={{ cursor: 'pointer', '&:hover': { bgcolor: 'action.hover' }, borderRadius: 1, py: 0.5 }}
            >
              {STEPS[3].label}
            </StepLabel>
            <StepContent>
              {layout.controls.placement === 'per_compound' ? (
                /* Strip layout: choose between grouped replicates or explicit naming */
                <Box>
                  <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
                    You have configured <strong>{layout.strip_layout?.strips_per_row || 2} strips per row</strong>.
                    Choose how replicates are determined:
                  </Typography>

                  <FormControl component="fieldset">
                    <RadioGroup
                      value={layout.replicate.pattern === 'explicit' ? 'explicit' : 'grouped'}
                      onChange={(e) => {
                        if (e.target.value === 'explicit') {
                          updateLayout({
                            replicate: { ...layout.replicate, pattern: 'explicit', count: 1 },
                          });
                        } else {
                          updateLayout({
                            replicate: { ...layout.replicate, pattern: 'adjacent_columns', count: layout.strip_layout?.strips_per_row || 2 },
                          });
                        }
                      }}
                    >
                      <FormControlLabel
                        value="grouped"
                        control={<Radio size="small" />}
                        label={
                          <Box>
                            <Typography variant="body2">Grouped Replicates</Typography>
                            <Typography variant="caption" color="text.secondary">
                              Strips are grouped as replicates of the same compound
                            </Typography>
                          </Box>
                        }
                      />
                      <FormControlLabel
                        value="explicit"
                        control={<Radio size="small" />}
                        label={
                          <Box>
                            <Typography variant="body2">Explicit Naming</Typography>
                            <Typography variant="caption" color="text.secondary">
                              Each strip has its own compound name; replicates inferred from matching names
                            </Typography>
                          </Box>
                        }
                      />
                    </RadioGroup>
                  </FormControl>

                  {layout.replicate.pattern === 'explicit' ? (
                    <Alert severity="info" sx={{ mt: 2 }}>
                      <Typography variant="body2">
                        With {layout.strip_layout?.strips_per_row || 2} strips per row, you&apos;ll provide{' '}
                        <strong>{layout.strip_layout?.strips_per_row || 2} compound names per row</strong>.
                      </Typography>
                      <Typography variant="caption" color="text.secondary" sx={{ display: 'block', mt: 0.5 }}>
                        Strips with identical compound names will be treated as replicates during analysis.
                      </Typography>
                    </Alert>
                  ) : (
                    <>
                      <Box sx={{ mt: 2, mb: 1 }}>
                        <TextField
                          label="Replicates per Compound"
                          type="number"
                          size="small"
                          value={layout.replicate.count}
                          onChange={(e) => updateLayout({
                            replicate: { ...layout.replicate, count: Number(e.target.value) },
                          })}
                          inputProps={{ min: 1, max: layout.strip_layout?.strips_per_row || 4 }}
                          sx={{ width: 200 }}
                          helperText={`${Math.floor((layout.strip_layout?.strips_per_row || 2) / layout.replicate.count)} compound(s) × ${layout.replicate.count} replicate(s) per row`}
                        />
                      </Box>
                      <Typography variant="body2" color="text.secondary" sx={{ mt: 1 }}>
                        Example with {layout.replicate.count} replicates:
                      </Typography>
                      <Typography variant="body2" color="text.secondary" sx={{ mt: 0.5, fontFamily: 'monospace' }}>
                        • A1-A12 = Compound 1, replicate 1
                      </Typography>
                      <Typography variant="body2" color="text.secondary" sx={{ fontFamily: 'monospace' }}>
                        • A13-A24 = Compound 1, replicate 2
                      </Typography>
                    </>
                  )}
                </Box>
              ) : (
                /* Standard layout: configure replicate pattern */
                <>
                  <Box sx={{ mb: 2 }}>
                    <TextField
                      label="Replicate Count"
                      type="number"
                      size="small"
                      value={layout.replicate.count}
                      onChange={(e) => updateLayout({
                        replicate: { ...layout.replicate, count: Number(e.target.value) },
                      })}
                      inputProps={{ min: 1, max: 8 }}
                      sx={{ width: 200 }}
                    />
                  </Box>

                  <FormControl component="fieldset">
                    <FormLabel>Replicate Pattern</FormLabel>
                    <RadioGroup
                      value={layout.replicate.pattern}
                      onChange={(e) => updateLayout({
                        replicate: { ...layout.replicate, pattern: e.target.value as ReplicatePattern },
                      })}
                    >
                      <FormControlLabel
                        value="adjacent_rows"
                        control={<Radio size="small" />}
                        label={
                          <Box>
                            <Typography variant="body2">Adjacent Rows</Typography>
                            <Typography variant="caption" color="text.secondary">
                              A1, B1 are replicates of compound 1
                            </Typography>
                          </Box>
                        }
                      />
                      <FormControlLabel
                        value="adjacent_columns"
                        control={<Radio size="small" />}
                        label={
                          <Box>
                            <Typography variant="body2">Adjacent Columns</Typography>
                            <Typography variant="caption" color="text.secondary">
                              A1, A2 are replicates of compound 1
                            </Typography>
                          </Box>
                        }
                      />
                      <FormControlLabel
                        value="grouped_rows"
                        control={<Radio size="small" />}
                        label={
                          <Box>
                            <Typography variant="body2">Grouped Rows</Typography>
                            <Typography variant="caption" color="text.secondary">
                              Row A and row B are replicates of the same compounds
                            </Typography>
                          </Box>
                        }
                      />
                      <FormControlLabel
                        value="interleaved_rows"
                        control={<Radio size="small" />}
                        label={
                          <Box>
                            <Typography variant="body2">Interleaved Rows</Typography>
                            <Typography variant="caption" color="text.secondary">
                              Alternating rows for the same compound
                            </Typography>
                          </Box>
                        }
                      />
                    </RadioGroup>
                  </FormControl>
                </>
              )}

              <Box sx={{ mt: 2 }}>
                <Button onClick={handleBack} sx={{ mr: 1 }}>Back</Button>
                <Button variant="contained" onClick={handleNext}>Continue</Button>
              </Box>
            </StepContent>
          </Step>

          {/* Step 5: Compound Source */}
          <Step completed={activeStep > 4}>
            <StepLabel
              onClick={() => setActiveStep(4)}
              sx={{ cursor: 'pointer', '&:hover': { bgcolor: 'action.hover' }, borderRadius: 1, py: 0.5 }}
            >
              {STEPS[4].label}
            </StepLabel>
            <StepContent>
              <FormControl component="fieldset">
                <FormLabel>How are compounds identified?</FormLabel>
                <RadioGroup
                  value={layout.compound_source.type}
                  onChange={(e) => updateLayout({
                    compound_source: { ...layout.compound_source, type: e.target.value as CompoundSourceType },
                  })}
                >
                  <FormControlLabel
                    value="row_order"
                    control={<Radio size="small" />}
                    label={
                      <Box>
                        <Typography variant="body2">Row Order</Typography>
                        <Typography variant="caption" color="text.secondary">
                          Compounds assigned sequentially by row position
                        </Typography>
                      </Box>
                    }
                  />
                  <FormControlLabel
                    value="column_header"
                    control={<Radio size="small" />}
                    label={
                      <Box>
                        <Typography variant="body2">Column Header</Typography>
                        <Typography variant="caption" color="text.secondary">
                          Compound IDs read from column headers in import file
                        </Typography>
                      </Box>
                    }
                  />
                  <FormControlLabel
                    value="row_header"
                    control={<Radio size="small" />}
                    label={
                      <Box>
                        <Typography variant="body2">Row Header</Typography>
                        <Typography variant="caption" color="text.secondary">
                          Compound IDs read from row labels in import file
                        </Typography>
                      </Box>
                    }
                  />
                  <FormControlLabel
                    value="adjacent_column"
                    control={<Radio size="small" />}
                    label={
                      <Box>
                        <Typography variant="body2">Adjacent Column (Right of Data)</Typography>
                        <Typography variant="caption" color="text.secondary">
                          Compound names in the column immediately after the data region
                        </Typography>
                      </Box>
                    }
                  />
                  <FormControlLabel
                    value="plate_map_file"
                    control={<Radio size="small" />}
                    label={
                      <Box>
                        <Typography variant="body2">Plate Map File</Typography>
                        <Typography variant="caption" color="text.secondary">
                          Separate file defines compound positions
                        </Typography>
                      </Box>
                    }
                  />
                </RadioGroup>
              </FormControl>

              {layout.compound_source.type === 'column_header' && (
                <Box sx={{ mt: 2 }}>
                  <TextField
                    label="ID Column Name"
                    size="small"
                    value={layout.compound_source.id_column || ''}
                    onChange={(e) => updateLayout({
                      compound_source: { ...layout.compound_source, id_column: e.target.value },
                    })}
                    placeholder="e.g., Compound_ID"
                    sx={{ width: 250 }}
                  />
                </Box>
              )}

              {/* Compound name row specification for stripe layouts */}
              {layout.controls.placement === 'per_compound' && (
                <Box sx={{ mt: 3 }}>
                  <Typography variant="subtitle2" gutterBottom>
                    Compound Name Location
                  </Typography>
                  <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
                    Specify the Excel row where compound names start. Names are read from the
                    left-most control column of each stripe, one per data row.
                  </Typography>
                  <TextField
                    label="Compound Name Row"
                    size="small"
                    value={compoundNameRowInput}
                    onChange={(e) => handleCompoundNameRowChange(e.target.value)}
                    placeholder="e.g., 25 or +2"
                    helperText={
                      layout.compound_source.compound_name_row
                        ? `Stored as absolute row ${layout.compound_source.compound_name_row}`
                        : `Enter absolute row (e.g., 25) or relative offset (e.g., +2 = 2 rows below plate row ${getPlateBottomRow()})`
                    }
                    sx={{ width: 250 }}
                  />
                  {layout.compound_source.compound_name_row && (
                    <Alert severity="info" sx={{ mt: 1 }}>
                      <Typography variant="body2">
                        For a {layout.strip_layout?.strips_per_row || 2}-stripe layout, compound names will be read from:
                      </Typography>
                      <Typography variant="body2" sx={{ fontFamily: 'monospace', mt: 0.5 }}>
                        {Array.from({ length: layout.strip_layout?.strips_per_row || 2 }, (_, i) => {
                          const stripStartCol = 1 + i * (layout.strip_layout?.strip_width || 12);
                          return `Stripe ${i + 1}: Column ${stripStartCol}`;
                        }).join(', ')}
                      </Typography>
                      <Typography variant="body2" sx={{ mt: 0.5 }}>
                        Starting at row {layout.compound_source.compound_name_row}
                      </Typography>
                    </Alert>
                  )}
                </Box>
              )}

              <Box sx={{ mt: 2 }}>
                <Button onClick={handleBack} sx={{ mr: 1 }}>Back</Button>
                <Button variant="contained" onClick={handleNext}>Continue</Button>
              </Box>
            </StepContent>
          </Step>

          {/* Step 6: Data Import (Spreadsheet Origin) */}
          <Step completed={activeStep > 5}>
            <StepLabel
              onClick={() => setActiveStep(5)}
              sx={{ cursor: 'pointer', '&:hover': { bgcolor: 'action.hover' }, borderRadius: 1, py: 0.5 }}
            >
              {STEPS[5].label}
            </StepLabel>
            <StepContent>
              <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
                Configure where plate data starts in imported Excel files.
                This is typically where cell A1 of the plate appears in the spreadsheet.
              </Typography>

              <Grid container spacing={2}>
                <Grid size={{ xs: 6 }}>
                  <TextField
                    label="Column"
                    size="small"
                    fullWidth
                    value={layout.spreadsheet_origin?.column || 'A'}
                    onChange={(e) => {
                      const col = e.target.value.toUpperCase().replace(/[^A-Z]/g, '');
                      updateLayout({
                        spreadsheet_origin: {
                          column: col || 'A',
                          row: layout.spreadsheet_origin?.row || 1,
                        },
                      });
                    }}
                    placeholder="A"
                    helperText="Excel column letter (e.g., A, B, AA)"
                    inputProps={{ maxLength: 3, style: { textTransform: 'uppercase' } }}
                  />
                </Grid>
                <Grid size={{ xs: 6 }}>
                  <TextField
                    label="Row"
                    size="small"
                    type="number"
                    fullWidth
                    value={layout.spreadsheet_origin?.row || 1}
                    onChange={(e) => updateLayout({
                      spreadsheet_origin: {
                        column: layout.spreadsheet_origin?.column || 'A',
                        row: Math.max(1, Number(e.target.value)),
                      },
                    })}
                    helperText="Starting row number (1-indexed)"
                    inputProps={{ min: 1 }}
                  />
                </Grid>
              </Grid>

              <Alert severity="info" sx={{ mt: 2 }}>
                <Typography variant="body2">
                  The plate data starts at cell <strong>{layout.spreadsheet_origin?.column || 'A'}{layout.spreadsheet_origin?.row || 1}</strong>.
                </Typography>
                <Typography variant="caption" color="text.secondary">
                  This cell should contain the value for well A1 of your plate.
                </Typography>
              </Alert>

              <Box sx={{ mt: 2 }}>
                <Button onClick={handleBack} sx={{ mr: 1 }}>Back</Button>
                <Button variant="contained" color="success" onClick={() => onChange(layout)}>
                  Apply Layout
                </Button>
              </Box>
            </StepContent>
          </Step>
        </Stepper>
      </Grid>

      {/* Preview panel */}
      {showPreview && (
        <Grid size={{ xs: 12, md: 5 }}>
          <Box sx={{ position: 'sticky', top: 16 }}>
            <Typography variant="subtitle2" gutterBottom>
              Layout Preview
            </Typography>
            <PlatePreview layout={layout} width={350} height={250} showSeriesRanges />
          </Box>
        </Grid>
      )}
    </Grid>
  );
}
