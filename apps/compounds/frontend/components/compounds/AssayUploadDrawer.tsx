'use client';

import { useState, useCallback, useMemo } from 'react';
import { useRouter } from 'next/navigation';
import {
  Drawer,
  Box,
  Typography,
  Button,
  IconButton,
  Alert,
  CircularProgress,
  Autocomplete,
  TextField,
  Divider,
  Stepper,
  Step,
  StepLabel,
  Paper,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Chip,
  Tooltip,
  FormControlLabel,
  Switch,
} from '@mui/material';
import {
  Close,
  CloudUpload,
  Science,
  CheckCircle,
  ArrowForward,
  ArrowBack,
  GridOn,
  Warning,
  Error as ErrorIcon,
  Analytics,
  Palette,
} from '@mui/icons-material';
import { useCompoundsApi } from '@/lib/compounds/api';
import { PlateHeatMapDialog } from './PlateHeatMap';
import type { PlateLayout, PlateFormat } from '@/types/compounds/models';

interface Target {
  id: string;
  name: string;
}

interface Protocol {
  id: string;
  name: string;
  plate_layout?: string | null;  // FK ID
  plate_layout_config?: Partial<PlateLayout> | null;  // Denormalized config
}

interface AssayUploadDrawerProps {
  open: boolean;
  onClose: () => void;
  protocolId: string;
  protocolName: string;
  onAssayCreated?: (assayId: string) => void;
}

// Spreadsheet data as 2D grid (row-major)
interface SpreadsheetGrid {
  cells: (string | number | null)[][];
  fileName: string;
  sheetName: string;
}

// Extracted data series from plate layout
interface ExtractedSeries {
  row: number;  // Plate row (0-indexed)
  rowLetter: string;
  stripIndex?: number;  // Strip index within the row (0-indexed), undefined for non-strip layouts
  compoundGroupIndex?: number;  // Compound group index for paired strips (e.g., 0 for strips 1+2, 1 for strips 3+4)
  compoundName: string | null;
  dataValues: (number | null)[];  // Data values per concentration
  minControlValues: (number | null)[];  // Raw control values for display
  maxControlValues: (number | null)[];  // Raw control values for display
  minControl: number | null;  // Averaged min control for analysis
  maxControl: number | null;  // Averaged max control for analysis
  startColumn: number;
  endColumn: number;
  hasIssues: boolean;
  issues: string[];
}

const STEPS = ['Upload Excel', 'Review Extraction', 'Create Assay'];

const PLATE_DIMENSIONS: Record<PlateFormat, { rows: number; cols: number }> = {
  24: { rows: 4, cols: 6 },
  96: { rows: 8, cols: 12 },
  384: { rows: 16, cols: 24 },
  1536: { rows: 32, cols: 48 },
};

function indexToRowLetter(index: number): string {
  return String.fromCharCode('A'.charCodeAt(0) + index);
}

function excelColumnToIndex(col: string): number {
  let result = 0;
  for (let i = 0; i < col.length; i++) {
    result = result * 26 + (col.charCodeAt(i) - 'A'.charCodeAt(0) + 1);
  }
  return result - 1; // 0-indexed
}

export function AssayUploadDrawer({
  open,
  onClose,
  protocolId,
  protocolName,
  onAssayCreated,
}: AssayUploadDrawerProps) {
  const router = useRouter();
  const api = useCompoundsApi();

  const [activeStep, setActiveStep] = useState(0);
  const [spreadsheetGrid, setSpreadsheetGrid] = useState<SpreadsheetGrid | null>(null);
  const [selectedFile, setSelectedFile] = useState<File | null>(null);
  const [selectedTarget, setSelectedTarget] = useState<Target | null>(null);
  const [labbook, setLabbook] = useState('');
  const [page, setPage] = useState('');
  const [comments, setComments] = useState('');
  const [creating, setCreating] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [runAnalysis, setRunAnalysis] = useState(true);  // Run curve fitting by default
  const [showHeatMap, setShowHeatMap] = useState(false);  // Heat map overlay visibility

  // Fetch targets and protocol details
  const { data: targetsData } = api.get<Target[]>('targets/');
  const { data: protocolData } = api.get<Protocol>(`protocols/${protocolId}/`);
  const targets = targetsData || [];
  const plateLayout = protocolData?.plate_layout_config as PlateLayout | null;

  // Extract data series from spreadsheet using plate layout
  const extractedSeries = useMemo((): ExtractedSeries[] => {
    if (!spreadsheetGrid || !plateLayout) return [];

    const { cells } = spreadsheetGrid;
    const format = plateLayout.plate_format || 384;
    const { rows: numRows, cols: numCols } = PLATE_DIMENSIONS[format];

    // Get origin offset
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

    // Process based on control placement
    if (plateLayout.controls?.placement === 'per_compound' && plateLayout.strip_layout) {
      // Strip layout: each row has embedded controls per strip
      // Each strip becomes a separate data series (technical replicates)
      const strip = plateLayout.strip_layout;
      // Ensure strip dimensions have sensible defaults
      const stripWidth = strip.strip_width || 12;
      const stripsPerRow = strip.strips_per_row || 2;
      const replicateCount = plateLayout.replicate?.count || 1;
      const isExplicitNaming = plateLayout.replicate?.pattern === 'explicit';

      // How many distinct compounds per row:
      // - Explicit naming: each strip gets its own name (strips_per_row compounds)
      // - Grouped replicates: e.g., 4 strips / 2 replicates = 2 compounds
      const compoundsPerRow = isExplicitNaming
        ? stripsPerRow
        : Math.max(1, Math.floor(stripsPerRow / replicateCount));

      for (let plateRow = startRowIdx; plateRow <= endRowIdx; plateRow++) {
        const spreadsheetRow = originRow + plateRow;
        if (spreadsheetRow >= cells.length) continue;

        const rowData = cells[spreadsheetRow] || [];
        const rowLetter = indexToRowLetter(plateRow);

        // Process each strip as a separate data series
        for (let stripIdx = 0; stripIdx < stripsPerRow; stripIdx++) {
          // Determine which compound group this strip belongs to:
          // - Explicit naming: each strip is its own compound group
          // - Grouped replicates: strips are grouped (e.g., strips 0,1 -> compound 0)
          const compoundGroupIdx = isExplicitNaming
            ? stripIdx
            : Math.floor(stripIdx / replicateCount);
          const stripStartCol = originCol + stripIdx * stripWidth;

          // Get compound name based on source type and compound group
          let compoundName: string | null = null;

          // New naming schema: compound_name_row specifies where names start
          // Names are in the left-most control column of each stripe
          if (plateLayout.compound_source?.compound_name_row) {
            // compound_name_row is absolute Excel row (1-indexed)
            // For each data row, read from corresponding name row
            const nameRowOffset = plateRow - startRowIdx;  // How far into data region
            const nameRow = plateLayout.compound_source.compound_name_row - 1 + nameRowOffset;  // 0-indexed
            const nameCol = stripStartCol;  // Left-most control column of this stripe

            if (nameRow >= 0 && nameRow < cells.length) {
              const nameRowData = cells[nameRow] || [];
              if (nameCol < nameRowData.length) {
                const nameVal = nameRowData[nameCol];
                compoundName = nameVal ? String(nameVal) : null;
              }
            }
          } else if (plateLayout.compound_source?.type === 'row_header') {
            // Legacy: For row_header with multiple compounds per row, read from column before plate
            // Each compound group reads from a different column (offset by compoundGroupIdx)
            const nameColIdx = originCol - compoundsPerRow + compoundGroupIdx;
            if (nameColIdx >= 0 && nameColIdx < rowData.length) {
              const nameVal = rowData[nameColIdx];
              compoundName = nameVal ? String(nameVal) : null;
            }
          } else if (plateLayout.compound_source?.type === 'adjacent_column') {
            // For adjacent_column, read compound names from columns after ALL plate data
            const totalStripWidth = stripWidth * stripsPerRow;
            const nameColIdx = originCol + totalStripWidth + compoundGroupIdx;
            if (nameColIdx < rowData.length) {
              const nameVal = rowData[nameColIdx];
              compoundName = nameVal ? String(nameVal) : null;
            }
          } else if (plateLayout.compound_source?.type === 'row_order') {
            // For row_order, generate sequential names
            // - Explicit naming: A1, A2, A3, A4 for 4 strips (each strip gets unique name)
            // - Grouped replicates: A1, A1, A2, A2 for 4 strips with 2 replicates
            if (compoundsPerRow > 1) {
              compoundName = `${rowLetter}${compoundGroupIdx + 1}`;
            } else {
              compoundName = rowLetter;
            }
          }
          const issues: string[] = [];

          // Extract min controls for this strip
          const minControlValues: (number | null)[] = [];
          for (let i = 0; i < strip.min_wells; i++) {
            const cellIdx = stripStartCol + i;
            const val = cellIdx < rowData.length ? rowData[cellIdx] : null;
            minControlValues.push(typeof val === 'number' ? val : null);
          }

          // Extract data values for this strip
          const dataValues: (number | null)[] = [];
          for (let i = 0; i < strip.data_wells; i++) {
            const cellIdx = stripStartCol + strip.min_wells + i;
            const val = cellIdx < rowData.length ? rowData[cellIdx] : null;
            dataValues.push(typeof val === 'number' ? val : null);
          }

          // Extract max controls for this strip
          const maxControlValues: (number | null)[] = [];
          for (let i = 0; i < strip.max_wells; i++) {
            const cellIdx = stripStartCol + strip.min_wells + strip.data_wells + i;
            const val = cellIdx < rowData.length ? rowData[cellIdx] : null;
            maxControlValues.push(typeof val === 'number' ? val : null);
          }

          // Average min control values for this strip
          const validMinValues = minControlValues.filter((v): v is number => v !== null);
          const averagedMinControl = validMinValues.length > 0
            ? validMinValues.reduce((a, b) => a + b, 0) / validMinValues.length
            : null;

          // Average max control values for this strip
          const validMaxValues = maxControlValues.filter((v): v is number => v !== null);
          const averagedMaxControl = validMaxValues.length > 0
            ? validMaxValues.reduce((a, b) => a + b, 0) / validMaxValues.length
            : null;

          // Validate data
          if (dataValues.every(v => v === null)) {
            issues.push('No numeric data found');
          }
          if (averagedMinControl === null) {
            issues.push('Missing min controls');
          }
          if (averagedMaxControl === null) {
            issues.push('Missing max controls');
          }

          // Build final data array: [min_control, data1, data2, ..., dataN, max_control]
          const finalDataValues: (number | null)[] = [
            averagedMinControl,
            ...dataValues,
            averagedMaxControl,
          ];

          series.push({
            row: plateRow,
            rowLetter,
            stripIndex: stripIdx,
            compoundGroupIndex: compoundGroupIdx,
            compoundName,
            dataValues: finalDataValues,
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
    } else {
      // Standard edge controls layout
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

        // Get compound name based on compound_source type
        let compoundName: string | null = null;
        if (plateLayout.compound_source?.type === 'row_header') {
          // Compound name in column before plate data
          const nameColIdx = originCol - 1;
          if (nameColIdx >= 0 && nameColIdx < rowData.length) {
            const nameVal = rowData[nameColIdx];
            compoundName = nameVal ? String(nameVal) : null;
          }
        } else if (plateLayout.compound_source?.type === 'adjacent_column') {
          // Compound name in column immediately after the data region
          const nameColIdx = originCol + endCol;  // endCol is already relative, so this is the next column
          if (nameColIdx < rowData.length) {
            const nameVal = rowData[nameColIdx];
            compoundName = nameVal ? String(nameVal) : null;
          }
        }

        // Compute averaged control values
        const validMinValues = minControlValues.filter((v): v is number => v !== null);
        const averagedMinControl = validMinValues.length > 0
          ? validMinValues.reduce((a, b) => a + b, 0) / validMinValues.length
          : null;

        const validMaxValues = maxControlValues.filter((v): v is number => v !== null);
        const averagedMaxControl = validMaxValues.length > 0
          ? validMaxValues.reduce((a, b) => a + b, 0) / validMaxValues.length
          : null;

        // Validate
        if (dataValues.every(v => v === null)) {
          issues.push('No numeric data found');
        }

        // Build final data array: [min_control, data1, data2, ..., dataN, max_control]
        const finalDataValues: (number | null)[] = [
          averagedMinControl,
          ...dataValues,
          averagedMaxControl,
        ];

        series.push({
          row: plateRow,
          rowLetter: indexToRowLetter(plateRow),
          compoundName,
          dataValues: finalDataValues,  // Controls at first/last positions
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

    return series;
  }, [spreadsheetGrid, plateLayout]);

  const validSeriesCount = extractedSeries.filter(s => !s.hasIssues).length;
  const hasPlateLayout = plateLayout && plateLayout.plate_format;

  const handleFileSelected = useCallback((file: File, grid: SpreadsheetGrid) => {
    setSelectedFile(file);
    setSpreadsheetGrid(grid);
    setActiveStep(1);
    setError(null);
  }, []);

  const handleBack = () => {
    setActiveStep((prev) => Math.max(prev - 1, 0));
  };

  const handleNext = () => {
    setActiveStep((prev) => Math.min(prev + 1, STEPS.length - 1));
  };

  const handleCreateAssay = async () => {
    if (!selectedFile) {
      setError('No file selected');
      return;
    }

    setCreating(true);
    setError(null);

    try {
      // Step 1: Create the assay
      const formData = new FormData();
      formData.append('protocol', protocolId);
      formData.append('data_file', selectedFile);

      if (selectedTarget) {
        formData.append('target', selectedTarget.id);
      }
      if (labbook) {
        formData.append('labbook_number', labbook);
      }
      if (page) {
        formData.append('page_number', page);
      }
      if (comments) {
        formData.append('comments', comments);
      }

      const newAssay = await api.upload<{ id: string }>('assays/', formData);

      // Step 2: Create data series sequentially to avoid SQLite locking
      const validSeries = extractedSeries.filter(s => !s.hasIssues);
      for (const series of validSeries) {
        // extracted_data format: [min_control, data1, ..., dataN, max_control]
        const seriesData = {
          assay: newAssay.id,
          compound_name: series.compoundName,
          row: series.row,
          start_column: series.startColumn,
          end_column: series.endColumn,
          extracted_data: series.dataValues,  // Controls embedded at first/last positions
        };

        try {
          await api.post('data-series/', seriesData);
        } catch {
          const seriesLabel = series.stripIndex !== undefined
            ? `${series.rowLetter}-${series.stripIndex + 1}`
            : series.rowLetter;
          console.warn(`Failed to create data series for ${seriesLabel}`);
        }
      }

      // Step 3: Optionally trigger curve fitting analysis
      console.log('Analysis check:', { runAnalysis, validSeriesCount: validSeries.length });
      if (runAnalysis && validSeries.length > 0) {
        console.log('Triggering analysis for assay:', newAssay.id);
        try {
          const analysisResult = await api.post(`assays/${newAssay.id}/analyse_all/`, {});
          console.log('Analysis completed:', analysisResult);
        } catch (err) {
          console.warn('Analysis request failed, but assay was created');
          console.error('Analysis failed:', err);
        }
      } else {
        console.log('Skipping analysis:', { runAnalysis, validSeriesCount: validSeries.length });
      }

      // Reset state
      setActiveStep(0);
      setSpreadsheetGrid(null);
      setSelectedFile(null);
      setSelectedTarget(null);
      setLabbook('');
      setPage('');
      setComments('');

      // Notify parent and navigate
      onAssayCreated?.(newAssay.id);
      onClose();
      router.push(`/assays/${newAssay.id}`);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to create assay');
    } finally {
      setCreating(false);
    }
  };

  const handleClose = () => {
    if (!creating) {
      onClose();
    }
  };

  const renderStepContent = () => {
    switch (activeStep) {
      case 0:
        return (
          <Box sx={{ py: 2 }}>
            {!hasPlateLayout && (
              <Alert severity="warning" sx={{ mb: 2 }}>
                <Typography variant="body2">
                  No plate layout configured for this protocol. Please configure a plate layout first
                  to enable automatic data extraction.
                </Typography>
              </Alert>
            )}
            <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
              Upload an Excel file containing your plate reader data.
              The data will be extracted according to the protocol&apos;s plate layout.
            </Typography>
            <FileDropZone
              onFileSelected={handleFileSelected}
              selectedFile={selectedFile}
              spreadsheetGrid={spreadsheetGrid}
            />
          </Box>
        );

      case 1:
        return (
          <Box sx={{ py: 2 }}>
            <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
              Preview of data series extracted from your spreadsheet using the plate layout.
            </Typography>

            {!hasPlateLayout ? (
              <Alert severity="error" sx={{ mb: 2 }}>
                No plate layout configured. Cannot extract data series.
              </Alert>
            ) : (
              <>
                {/* Summary */}
                <Paper sx={{ p: 2, mb: 2, bgcolor: 'grey.50' }}>
                  <Box sx={{ display: 'flex', gap: 2, flexWrap: 'wrap', alignItems: 'center' }}>
                    <Chip
                      icon={<GridOn />}
                      label={`${plateLayout?.plate_format}-well plate`}
                      size="small"
                    />
                    <Chip
                      icon={<CheckCircle />}
                      label={`${validSeriesCount} valid series`}
                      size="small"
                      color="success"
                    />
                    {extractedSeries.length - validSeriesCount > 0 && (
                      <Chip
                        icon={<Warning />}
                        label={`${extractedSeries.length - validSeriesCount} with issues`}
                        size="small"
                        color="warning"
                      />
                    )}
                    <Button
                      size="small"
                      startIcon={<Palette />}
                      onClick={() => setShowHeatMap(true)}
                      variant="outlined"
                      sx={{ ml: 'auto' }}
                    >
                      Heat Map
                    </Button>
                  </Box>
                  <Typography variant="caption" color="text.secondary" sx={{ display: 'block', mt: 1 }}>
                    Origin: cell {plateLayout?.spreadsheet_origin?.column || 'A'}{plateLayout?.spreadsheet_origin?.row || 1}
                  </Typography>
                </Paper>

                {/* Extracted series table */}
                <TableContainer component={Paper} sx={{ maxHeight: 400 }}>
                  <Table size="small" stickyHeader>
                    <TableHead>
                      <TableRow>
                        <TableCell sx={{ fontWeight: 600 }}>Series</TableCell>
                        <TableCell sx={{ fontWeight: 600 }}>Compound</TableCell>
                        <TableCell sx={{ fontWeight: 600 }}>Data Points</TableCell>
                        <TableCell sx={{ fontWeight: 600 }}>Controls</TableCell>
                        <TableCell sx={{ fontWeight: 600 }}>Status</TableCell>
                      </TableRow>
                    </TableHead>
                    <TableBody>
                      {extractedSeries.map((series, idx) => (
                        <TableRow
                          key={`${series.row}-${series.stripIndex ?? 0}-${idx}`}
                          sx={{ bgcolor: series.hasIssues ? 'warning.50' : undefined }}
                        >
                          <TableCell>
                            <Typography fontFamily="monospace">
                              {series.stripIndex !== undefined
                                ? `${series.rowLetter}-${series.stripIndex + 1}`
                                : series.rowLetter}
                            </Typography>
                          </TableCell>
                          <TableCell>
                            {series.compoundName || (
                              <Typography color="text.secondary" variant="body2">
                                Row {series.row + 1}
                              </Typography>
                            )}
                          </TableCell>
                          <TableCell>
                            <Typography variant="body2" fontFamily="monospace">
                              {series.dataValues.filter(v => v !== null).length} values
                            </Typography>
                          </TableCell>
                          <TableCell>
                            <Box sx={{ display: 'flex', gap: 0.5 }}>
                              <Chip
                                size="small"
                                label={`min: ${series.minControlValues.filter(v => v !== null).length}`}
                                color={series.minControlValues.some(v => v !== null) ? 'info' : 'default'}
                              />
                              <Chip
                                size="small"
                                label={`max: ${series.maxControlValues.filter(v => v !== null).length}`}
                                color={series.maxControlValues.some(v => v !== null) ? 'primary' : 'default'}
                              />
                            </Box>
                          </TableCell>
                          <TableCell>
                            {series.hasIssues ? (
                              <Tooltip title={series.issues.join(', ')}>
                                <Chip
                                  size="small"
                                  icon={<Warning />}
                                  label="Issues"
                                  color="warning"
                                />
                              </Tooltip>
                            ) : (
                              <Chip
                                size="small"
                                icon={<CheckCircle />}
                                label="OK"
                                color="success"
                              />
                            )}
                          </TableCell>
                        </TableRow>
                      ))}
                    </TableBody>
                  </Table>
                </TableContainer>
              </>
            )}
          </Box>
        );

      case 2:
        return (
          <Box sx={{ py: 2 }}>
            <Typography variant="body2" color="text.secondary" sx={{ mb: 3 }}>
              Configure the assay metadata and create the assay with extracted data series.
            </Typography>

            {/* Target selection */}
            <Autocomplete
              options={targets}
              getOptionLabel={(option) => option.name}
              value={selectedTarget}
              onChange={(_, newValue) => setSelectedTarget(newValue)}
              renderInput={(params) => (
                <TextField
                  {...params}
                  label="Target"
                  size="small"
                  helperText="Select the target being tested (optional)"
                />
              )}
              sx={{ mb: 2 }}
            />

            {/* Lab book info */}
            <Box sx={{ display: 'flex', gap: 2, mb: 2 }}>
              <TextField
                label="Lab Book #"
                value={labbook}
                onChange={(e) => setLabbook(e.target.value)}
                size="small"
                type="number"
                sx={{ flex: 1 }}
              />
              <TextField
                label="Page #"
                value={page}
                onChange={(e) => setPage(e.target.value)}
                size="small"
                type="number"
                sx={{ flex: 1 }}
              />
            </Box>

            {/* Comments */}
            <TextField
              label="Comments"
              value={comments}
              onChange={(e) => setComments(e.target.value)}
              multiline
              rows={2}
              fullWidth
              size="small"
              sx={{ mb: 2 }}
            />

            {/* Run analysis option */}
            <Paper variant="outlined" sx={{ p: 2, mb: 2 }}>
              <FormControlLabel
                control={
                  <Switch
                    checked={runAnalysis}
                    onChange={(e) => setRunAnalysis(e.target.checked)}
                    color="primary"
                  />
                }
                label={
                  <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                    <Analytics color={runAnalysis ? 'primary' : 'disabled'} />
                    <Box>
                      <Typography variant="body2" fontWeight={500}>
                        Run curve fitting analysis
                      </Typography>
                      <Typography variant="caption" color="text.secondary">
                        Automatically fit dose-response curves using the protocol&apos;s analysis method
                      </Typography>
                    </Box>
                  </Box>
                }
                sx={{ m: 0 }}
              />
            </Paper>

            <Divider sx={{ my: 2 }} />

            {/* Summary */}
            <Paper sx={{ p: 2, bgcolor: 'grey.50' }}>
              <Typography variant="subtitle2" gutterBottom>
                Summary
              </Typography>
              <Box sx={{ display: 'grid', gridTemplateColumns: '1fr 2fr', gap: 1 }}>
                <Typography variant="body2" color="text.secondary">File:</Typography>
                <Typography variant="body2">{spreadsheetGrid?.fileName}</Typography>

                <Typography variant="body2" color="text.secondary">Protocol:</Typography>
                <Typography variant="body2">{protocolName}</Typography>

                <Typography variant="body2" color="text.secondary">Target:</Typography>
                <Typography variant="body2">{selectedTarget?.name || 'Not specified'}</Typography>

                <Typography variant="body2" color="text.secondary">Data Series:</Typography>
                <Typography variant="body2">{validSeriesCount} to be created</Typography>
              </Box>
            </Paper>

            {error && (
              <Alert severity="error" sx={{ mt: 2 }}>
                {error}
              </Alert>
            )}
          </Box>
        );

      default:
        return null;
    }
  };

  return (
    <>
    <Drawer
      anchor="right"
      open={open}
      onClose={handleClose}
      PaperProps={{
        sx: { width: { xs: '100%', sm: 560 } },
      }}
    >
      {/* Header */}
      <Box
        sx={{
          p: 2,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'space-between',
          borderBottom: 1,
          borderColor: 'divider',
        }}
      >
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <CloudUpload color="primary" />
          <Typography variant="h6">Add Assay</Typography>
        </Box>
        <IconButton onClick={handleClose} disabled={creating} size="small">
          <Close />
        </IconButton>
      </Box>

      {/* Stepper */}
      <Box sx={{ px: 2, py: 2, borderBottom: 1, borderColor: 'divider' }}>
        <Stepper activeStep={activeStep} alternativeLabel>
          {STEPS.map((label) => (
            <Step key={label}>
              <StepLabel>{label}</StepLabel>
            </Step>
          ))}
        </Stepper>
      </Box>

      {/* Content */}
      <Box sx={{ flex: 1, overflow: 'auto', px: 2 }}>
        {renderStepContent()}
      </Box>

      {/* Footer */}
      <Box
        sx={{
          p: 2,
          borderTop: 1,
          borderColor: 'divider',
          display: 'flex',
          justifyContent: 'space-between',
        }}
      >
        <Button
          onClick={handleBack}
          disabled={activeStep === 0 || creating}
          startIcon={<ArrowBack />}
        >
          Back
        </Button>

        {activeStep < STEPS.length - 1 ? (
          <Button
            variant="contained"
            onClick={handleNext}
            disabled={activeStep === 0 && !spreadsheetGrid}
            endIcon={<ArrowForward />}
          >
            Continue
          </Button>
        ) : (
          <Button
            variant="contained"
            onClick={handleCreateAssay}
            disabled={creating || !spreadsheetGrid || validSeriesCount === 0}
            startIcon={creating ? <CircularProgress size={16} /> : <CheckCircle />}
          >
            {creating ? 'Creating...' : `Create Assay (${validSeriesCount} series)`}
          </Button>
        )}
      </Box>
    </Drawer>

    {/* Heat Map Dialog */}
    {spreadsheetGrid && plateLayout && (
      <PlateHeatMapDialog
        open={showHeatMap}
        onClose={() => setShowHeatMap(false)}
        cells={spreadsheetGrid.cells}
        plateLayout={plateLayout}
      />
    )}
    </>
  );
}

// File drop zone component
interface FileDropZoneProps {
  onFileSelected: (file: File, grid: SpreadsheetGrid) => void;
  selectedFile: File | null;
  spreadsheetGrid: SpreadsheetGrid | null;
}

function FileDropZone({ onFileSelected, selectedFile, spreadsheetGrid }: FileDropZoneProps) {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [dragOver, setDragOver] = useState(false);

  const parseFile = async (file: File) => {
    setLoading(true);
    setError(null);

    try {
      const XLSX = await import('xlsx');
      const buffer = await file.arrayBuffer();
      const workbook = XLSX.read(buffer, { type: 'array' });

      const sheetName = workbook.SheetNames[0];
      const worksheet = workbook.Sheets[sheetName];

      if (!worksheet) {
        throw new Error('No worksheet found');
      }

      // Get sheet as 2D array (preserving all cells including empty ones)
      const range = XLSX.utils.decode_range(worksheet['!ref'] || 'A1');
      const cells: (string | number | null)[][] = [];

      for (let row = range.s.r; row <= range.e.r; row++) {
        const rowData: (string | number | null)[] = [];
        for (let col = range.s.c; col <= range.e.c; col++) {
          const cellAddress = XLSX.utils.encode_cell({ r: row, c: col });
          const cell = worksheet[cellAddress];
          if (cell) {
            rowData.push(cell.v as string | number | null);
          } else {
            rowData.push(null);
          }
        }
        cells.push(rowData);
      }

      const grid: SpreadsheetGrid = {
        cells,
        fileName: file.name,
        sheetName,
      };

      onFileSelected(file, grid);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to parse file');
    } finally {
      setLoading(false);
    }
  };

  const handleDrop = (e: React.DragEvent) => {
    e.preventDefault();
    setDragOver(false);

    const file = e.dataTransfer.files[0];
    if (file) {
      parseFile(file);
    }
  };

  const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (file) {
      parseFile(file);
    }
  };

  if (selectedFile && spreadsheetGrid) {
    return (
      <Box
        sx={{
          p: 3,
          border: '2px solid',
          borderColor: 'success.main',
          borderRadius: 1,
          bgcolor: 'success.50',
          textAlign: 'center',
        }}
      >
        <CheckCircle sx={{ fontSize: 48, color: 'success.main', mb: 1 }} />
        <Typography fontWeight={600}>{spreadsheetGrid.fileName}</Typography>
        <Typography variant="body2" color="text.secondary">
          {spreadsheetGrid.cells.length} rows × {spreadsheetGrid.cells[0]?.length || 0} columns
        </Typography>
        <Typography variant="caption" color="text.secondary">
          Sheet: {spreadsheetGrid.sheetName}
        </Typography>
        <Box sx={{ mt: 2 }}>
          <Button
            size="small"
            onClick={() => {
              const input = document.getElementById('file-input-replace') as HTMLInputElement;
              input?.click();
            }}
          >
            Change File
          </Button>
          <input
            id="file-input-replace"
            type="file"
            hidden
            accept=".xlsx,.xls,.csv"
            onChange={handleFileChange}
          />
        </Box>
      </Box>
    );
  }

  return (
    <Box
      onDragOver={(e) => {
        e.preventDefault();
        setDragOver(true);
      }}
      onDragLeave={() => setDragOver(false)}
      onDrop={handleDrop}
      sx={{
        p: 4,
        border: '2px dashed',
        borderColor: dragOver ? 'primary.main' : 'grey.300',
        borderRadius: 1,
        bgcolor: dragOver ? 'primary.50' : 'grey.50',
        textAlign: 'center',
        cursor: 'pointer',
        transition: 'all 0.2s',
        '&:hover': {
          borderColor: 'primary.main',
          bgcolor: 'primary.50',
        },
      }}
      onClick={() => {
        const input = document.getElementById('file-input') as HTMLInputElement;
        input?.click();
      }}
    >
      {loading ? (
        <>
          <CircularProgress size={48} sx={{ mb: 2 }} />
          <Typography>Processing file...</Typography>
        </>
      ) : (
        <>
          <CloudUpload sx={{ fontSize: 48, color: 'grey.400', mb: 1 }} />
          <Typography fontWeight={500} gutterBottom>
            Drop Excel file here
          </Typography>
          <Typography variant="body2" color="text.secondary">
            or click to browse
          </Typography>
        </>
      )}
      <input
        id="file-input"
        type="file"
        hidden
        accept=".xlsx,.xls,.csv"
        onChange={handleFileChange}
      />
      {error && (
        <Alert severity="error" sx={{ mt: 2 }}>
          {error}
        </Alert>
      )}
    </Box>
  );
}
