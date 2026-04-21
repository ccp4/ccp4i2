'use client';

import { useState, useCallback, useMemo, useRef } from 'react';
import ExcelJS from 'exceljs';
import {
  Box,
  Paper,
  Typography,
  Button,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  IconButton,
  Alert,
  CircularProgress,
  Chip,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  SelectChangeEvent,
  Tooltip,
} from '@mui/material';
import {
  CloudUpload,
  Close,
  CheckCircle,
  Warning,
  TableChart,
  ArrowForward,
} from '@mui/icons-material';

export interface SpreadsheetRow {
  [key: string]: string | number | null;
}

/**
 * Hyperlink data extracted from Excel cells.
 * Maps column header -> hyperlink URL for each row.
 */
export interface SpreadsheetHyperlinks {
  [key: string]: string | undefined;
}

export interface SpreadsheetData {
  headers: string[];
  rows: SpreadsheetRow[];
  /** Hyperlinks extracted from Excel cells, indexed by row */
  hyperlinks?: SpreadsheetHyperlinks[];
  fileName: string;
  sheetName: string;
  /** Original file object for upload */
  originalFile?: File;
}

export interface FieldMapping {
  field: string;
  label: string;
  required?: boolean;
  type?: 'string' | 'number' | 'select';
}

interface SpreadsheetUploadProps {
  onDataLoaded?: (data: SpreadsheetData) => void;
  onMappingComplete?: (mapping: Record<string, string>, data: SpreadsheetData) => void;
  fieldMappings?: FieldMapping[];
  title?: string;
  previewRows?: number;
  showColumnMapping?: boolean;
}

function parseCSV(text: string): string[][] {
  const rows: string[][] = [];
  let current = '';
  let inQuotes = false;
  let row: string[] = [];
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (inQuotes) {
      if (ch === '"' && text[i + 1] === '"') { current += '"'; i++; }
      else if (ch === '"') { inQuotes = false; }
      else { current += ch; }
    } else {
      if (ch === '"') { inQuotes = true; }
      else if (ch === ',') { row.push(current.trim()); current = ''; }
      else if (ch === '\n' || (ch === '\r' && text[i + 1] === '\n')) {
        row.push(current.trim()); current = '';
        if (row.some(c => c !== '')) rows.push(row);
        row = [];
        if (ch === '\r') i++;
      } else { current += ch; }
    }
  }
  if (current || row.length > 0) {
    row.push(current.trim());
    if (row.some(c => c !== '')) rows.push(row);
  }
  return rows;
}

export function SpreadsheetUpload({
  onDataLoaded,
  onMappingComplete,
  fieldMappings = [],
  title = 'Upload Spreadsheet',
  previewRows = 10,
  showColumnMapping = true,
}: SpreadsheetUploadProps) {
  const [file, setFile] = useState<File | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [data, setData] = useState<SpreadsheetData | null>(null);
  const [columnMapping, setColumnMapping] = useState<Record<string, string>>({});
  const [sheets, setSheets] = useState<string[]>([]);
  const [selectedSheet, setSelectedSheet] = useState<string>('');
  const inputRef = useRef<HTMLInputElement>(null);

  const parseFile = useCallback(async (file: File, sheetName?: string) => {
    setLoading(true);
    setError(null);

    try {
      const buffer = await file.arrayBuffer();
      let sheetNames: string[];
      let targetSheet: string;
      let headers: string[];
      let jsonData: SpreadsheetRow[];
      let hyperlinkData: SpreadsheetHyperlinks[] = [];

      if (file.name.endsWith('.csv')) {
        // Parse CSV as text
        const text = new TextDecoder().decode(buffer);
        const parsed = parseCSV(text);
        if (parsed.length < 2) throw new Error('Spreadsheet appears to be empty');
        sheetNames = ['Sheet1'];
        targetSheet = 'Sheet1';
        headers = parsed[0];
        jsonData = parsed.slice(1).map(row => {
          const obj: SpreadsheetRow = {};
          headers.forEach((h, i) => { obj[h] = row[i] || null; });
          return obj;
        });
      } else {
        // Parse Excel with ExcelJS
        const wb = new ExcelJS.Workbook();
        await wb.xlsx.load(buffer);

        sheetNames = wb.worksheets.map(ws => ws.name);
        targetSheet = sheetName || sheetNames[0];
        const worksheet = wb.getWorksheet(targetSheet);
        if (!worksheet) throw new Error(`Sheet "${targetSheet}" not found`);

        // Read header row
        headers = [];
        const headerRow = worksheet.getRow(1);
        headerRow.eachCell({ includeEmpty: false }, (cell) => {
          headers.push(String(cell.value ?? ''));
        });

        // Read data rows, converting all values to strings for consistency
        // Also extract hyperlinks from cells (used for ELN references, document URLs)
        jsonData = [];
        hyperlinkData = [];
        for (let i = 2; i <= worksheet.rowCount; i++) {
          const row = worksheet.getRow(i);
          const obj: SpreadsheetRow = {};
          const rowHyperlinks: SpreadsheetHyperlinks = {};
          let hasData = false;
          headers.forEach((h, idx) => {
            const cell = row.getCell(idx + 1);
            if (cell.value != null) {
              // cell.text extracts display text for hyperlinks ({text, hyperlink}),
              // rich text ({richText}), formula results, etc. String(cell.value)
              // would yield "[object Object]" for these.
              obj[h] = cell.text;
              hasData = true;
            } else {
              obj[h] = null;
            }
            // Extract hyperlink if present
            if (cell.hyperlink) {
              // ExcelJS stores hyperlink in cell.hyperlink or cell.value.hyperlink
              const hyperlink = typeof cell.hyperlink === 'string'
                ? cell.hyperlink
                : (cell.hyperlink as any)?.target || (cell.hyperlink as any)?.address;
              if (hyperlink) {
                rowHyperlinks[h] = hyperlink;
              }
            }
          });
          if (hasData) {
            jsonData.push(obj);
            hyperlinkData.push(rowHyperlinks);
          }
        }
      }

      setSheets(sheetNames);
      setSelectedSheet(targetSheet);

      if (jsonData.length === 0) {
        throw new Error('Spreadsheet appears to be empty');
      }

      // For CSV files, hyperlinkData will be empty (no hyperlinks in CSV)
      // For Excel files, hyperlinkData contains extracted hyperlinks per row
      const spreadsheetData: SpreadsheetData = {
        headers,
        rows: jsonData,
        hyperlinks: hyperlinkData.length > 0 ? hyperlinkData : jsonData.map(() => ({})),
        fileName: file.name,
        sheetName: targetSheet,
        originalFile: file,
      };

      setData(spreadsheetData);

      // Auto-map columns if field mappings provided
      if (fieldMappings.length > 0) {
        const normalize = (s: string) =>
          s
            .toLowerCase()
            .replace(/\([^)]*\)/g, '')
            .replace(/[_\s/]+/g, ' ')
            .trim();

        const normHeaders = headers.map((h) => ({ orig: h, norm: normalize(h) }));
        const fieldNormsFor = (f: FieldMapping) =>
          [normalize(f.label), normalize(f.field)].filter((s) => s.length > 0);

        const autoMapping: Record<string, string> = {};
        const usedColumns = new Set<string>();

        // Pass 1: exact match against normalized label or field name.
        fieldMappings.forEach((field) => {
          const fieldNorms = fieldNormsFor(field);
          const match = normHeaders.find(
            (nh) => !usedColumns.has(nh.orig) && fieldNorms.includes(nh.norm)
          );
          if (match) {
            autoMapping[field.field] = match.orig;
            usedColumns.add(match.orig);
          }
        });

        // Pass 2: partial match, but only when exactly one unmapped field could
        // claim the column. This prevents short/ambiguous headers like "Batch"
        // from being silently assigned to "batch_amount" (or "batch_salt_code").
        const partialMatches = (headerNorm: string, fieldNorms: string[]) =>
          fieldNorms.some(
            (fn) =>
              fn.length > 0 &&
              headerNorm.length > 0 &&
              (headerNorm.includes(fn) || fn.includes(headerNorm))
          );

        fieldMappings.forEach((field) => {
          if (autoMapping[field.field]) return;
          const fieldNorms = fieldNormsFor(field);
          const candidate = normHeaders.find((nh) => {
            if (usedColumns.has(nh.orig)) return false;
            if (!partialMatches(nh.norm, fieldNorms)) return false;
            const competitors = fieldMappings.filter(
              (other) =>
                other.field !== field.field &&
                !autoMapping[other.field] &&
                partialMatches(nh.norm, fieldNormsFor(other))
            );
            return competitors.length === 0;
          });
          if (candidate) {
            autoMapping[field.field] = candidate.orig;
            usedColumns.add(candidate.orig);
          }
        });

        setColumnMapping(autoMapping);
      }

      onDataLoaded?.(spreadsheetData);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to parse file');
    } finally {
      setLoading(false);
    }
  }, [fieldMappings, onDataLoaded]);

  const handleFileChange = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) => {
      const selectedFile = event.target.files?.[0];
      // Reset so the same file can be picked again after an error.
      event.target.value = '';
      if (!selectedFile) return;

      // Validate file type
      const validTypes = [
        'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
        'application/vnd.ms-excel',
        'text/csv',
        'application/csv',
      ];
      const isValid =
        validTypes.includes(selectedFile.type) ||
        selectedFile.name.endsWith('.xlsx') ||
        selectedFile.name.endsWith('.xls') ||
        selectedFile.name.endsWith('.csv');

      if (!isValid) {
        setError('Please upload an Excel (.xlsx, .xls) or CSV file');
        return;
      }

      setFile(selectedFile);
      parseFile(selectedFile);
    },
    [parseFile]
  );

  const handleSheetChange = useCallback(
    (event: SelectChangeEvent) => {
      const newSheet = event.target.value;
      setSelectedSheet(newSheet);
      if (file) {
        parseFile(file, newSheet);
      }
    },
    [file, parseFile]
  );

  const handleMappingChange = useCallback(
    (field: string, column: string) => {
      setColumnMapping((prev) => ({
        ...prev,
        [field]: column,
      }));
    },
    []
  );

  const handleClearFile = useCallback(() => {
    setFile(null);
    setData(null);
    setError(null);
    setColumnMapping({});
    setSheets([]);
    setSelectedSheet('');
  }, []);

  const handleConfirmMapping = useCallback(() => {
    if (data && onMappingComplete) {
      onMappingComplete(columnMapping, data);
    }
  }, [data, columnMapping, onMappingComplete]);

  // Check if all required fields are mapped
  const requiredFieldsMapped = useMemo(() => {
    const requiredFields = fieldMappings.filter((f) => f.required);
    return requiredFields.every((f) => columnMapping[f.field]);
  }, [fieldMappings, columnMapping]);

  // Preview data with limited rows
  const previewData = useMemo(() => {
    if (!data) return null;
    return {
      ...data,
      rows: data.rows.slice(0, previewRows),
    };
  }, [data, previewRows]);

  return (
    <Paper sx={{ p: 3 }}>
      <Typography variant="h6" gutterBottom sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
        <TableChart color="primary" />
        {title}
      </Typography>

      {!data ? (
        // File upload area. Using an explicit ref + onClick (rather than
        // <Box component="label">) avoids a first-click-does-nothing symptom
        // seen under Next.js hydration, where the implicit label↔input
        // association isn't yet wired on the first render.
        <Box
          role="button"
          tabIndex={0}
          onClick={() => inputRef.current?.click()}
          onKeyDown={(e) => {
            if (e.key === 'Enter' || e.key === ' ') {
              e.preventDefault();
              inputRef.current?.click();
            }
          }}
          sx={{
            border: '2px dashed',
            borderColor: error ? 'error.main' : 'divider',
            borderRadius: 2,
            p: 4,
            textAlign: 'center',
            bgcolor: 'grey.50',
            cursor: 'pointer',
            '&:hover': { borderColor: 'primary.main', bgcolor: 'grey.100' },
          }}
        >
          <input
            ref={inputRef}
            type="file"
            hidden
            accept=".xlsx,.xls,.csv"
            onChange={handleFileChange}
          />
          {loading ? (
            <CircularProgress />
          ) : (
            <>
              <CloudUpload sx={{ fontSize: 48, color: 'primary.main', mb: 2 }} />
              <Typography variant="body1" gutterBottom>
                Drop Excel or CSV file here, or click to browse
              </Typography>
              <Typography variant="caption" color="text.secondary">
                Supported formats: .xlsx, .xls, .csv
              </Typography>
            </>
          )}
        </Box>
      ) : (
        // File loaded - show preview
        <Box>
          {/* File info header */}
          <Box
            sx={{
              display: 'flex',
              alignItems: 'center',
              gap: 2,
              mb: 2,
              p: 1.5,
              bgcolor: 'success.50',
              borderRadius: 1,
            }}
          >
            <CheckCircle color="success" />
            <Box sx={{ flex: 1 }}>
              <Typography variant="body2" fontWeight={600}>
                {data.fileName}
              </Typography>
              <Typography variant="caption" color="text.secondary">
                {data.rows.length} rows, {data.headers.length} columns
              </Typography>
            </Box>

            {/* Sheet selector if multiple sheets */}
            {sheets.length > 1 && (
              <FormControl size="small" sx={{ minWidth: 150 }}>
                <InputLabel>Sheet</InputLabel>
                <Select
                  value={selectedSheet}
                  label="Sheet"
                  onChange={handleSheetChange}
                >
                  {sheets.map((sheet) => (
                    <MenuItem key={sheet} value={sheet}>
                      {sheet}
                    </MenuItem>
                  ))}
                </Select>
              </FormControl>
            )}

            <IconButton onClick={handleClearFile} size="small">
              <Close />
            </IconButton>
          </Box>

          {/* Column mapping section */}
          {showColumnMapping && fieldMappings.length > 0 && (
            <Box sx={{ mb: 3 }}>
              <Typography variant="subtitle2" gutterBottom>
                Map Columns to Fields
              </Typography>
              <Box
                sx={{
                  display: 'grid',
                  gridTemplateColumns: 'repeat(auto-fill, minmax(280px, 1fr))',
                  gap: 2,
                }}
              >
                {fieldMappings.map((field) => (
                  <FormControl key={field.field} size="small" fullWidth>
                    <InputLabel>
                      {field.label}
                      {field.required && ' *'}
                    </InputLabel>
                    <Select
                      value={columnMapping[field.field] || ''}
                      label={`${field.label}${field.required ? ' *' : ''}`}
                      onChange={(e) =>
                        handleMappingChange(field.field, e.target.value)
                      }
                    >
                      <MenuItem value="">
                        <em>-- Not mapped --</em>
                      </MenuItem>
                      {data.headers.map((header) => (
                        <MenuItem key={header} value={header}>
                          {header}
                        </MenuItem>
                      ))}
                    </Select>
                  </FormControl>
                ))}
              </Box>

              {!requiredFieldsMapped && (
                <Alert severity="warning" sx={{ mt: 2 }}>
                  Please map all required fields before proceeding.
                </Alert>
              )}
            </Box>
          )}

          {/* Data preview table */}
          <Typography variant="subtitle2" gutterBottom>
            Preview ({Math.min(previewRows, data.rows.length)} of {data.rows.length} rows)
          </Typography>
          <TableContainer
            component={Paper}
            variant="outlined"
            sx={{ maxHeight: 400, mb: 2 }}
          >
            <Table size="small" stickyHeader>
              <TableHead>
                <TableRow>
                  <TableCell sx={{ fontWeight: 600, bgcolor: 'grey.100' }}>#</TableCell>
                  {previewData?.headers.map((header) => {
                    // Check if this column is mapped
                    const mappedField = Object.entries(columnMapping).find(
                      ([, col]) => col === header
                    );
                    return (
                      <TableCell
                        key={header}
                        sx={{
                          fontWeight: 600,
                          bgcolor: mappedField ? 'primary.50' : 'grey.100',
                          minWidth: 100,
                        }}
                      >
                        <Box>
                          {header}
                          {mappedField && (
                            <Chip
                              label={mappedField[0]}
                              size="small"
                              color="primary"
                              sx={{ ml: 1 }}
                            />
                          )}
                        </Box>
                      </TableCell>
                    );
                  })}
                </TableRow>
              </TableHead>
              <TableBody>
                {previewData?.rows.map((row, idx) => {
                  const rowHyperlinks = previewData.hyperlinks?.[idx];
                  return (
                    <TableRow key={idx} hover>
                      <TableCell sx={{ color: 'text.secondary' }}>{idx + 1}</TableCell>
                      {previewData.headers.map((header) => {
                        const value = row[header];
                        const url = rowHyperlinks?.[header];
                        if (value === null || value === undefined) {
                          return <TableCell key={header}>-</TableCell>;
                        }
                        const display = String(value);
                        return (
                          <TableCell key={header}>
                            {url ? (
                              <a
                                href={url}
                                target="_blank"
                                rel="noopener noreferrer"
                                title={url}
                                style={{ color: 'inherit' }}
                              >
                                {display}
                              </a>
                            ) : (
                              display
                            )}
                          </TableCell>
                        );
                      })}
                    </TableRow>
                  );
                })}
              </TableBody>
            </Table>
          </TableContainer>

          {/* Action buttons */}
          {showColumnMapping && onMappingComplete && (
            <Box sx={{ display: 'flex', justifyContent: 'flex-end', gap: 2 }}>
              <Button variant="outlined" onClick={handleClearFile}>
                Cancel
              </Button>
              <Button
                variant="contained"
                disabled={!requiredFieldsMapped}
                onClick={handleConfirmMapping}
                endIcon={<ArrowForward />}
              >
                Continue with Mapping
              </Button>
            </Box>
          )}
        </Box>
      )}

      {error && (
        <Alert severity="error" sx={{ mt: 2 }}>
          {error}
        </Alert>
      )}
    </Paper>
  );
}

/**
 * Standalone preview component for already-loaded data
 */
interface SpreadsheetPreviewProps {
  data: SpreadsheetData;
  maxRows?: number;
  highlightColumns?: string[];
}

export function SpreadsheetPreview({
  data,
  maxRows = 20,
  highlightColumns = [],
}: SpreadsheetPreviewProps) {
  const previewRows = data.rows.slice(0, maxRows);

  return (
    <TableContainer component={Paper} variant="outlined" sx={{ maxHeight: 500 }}>
      <Table size="small" stickyHeader>
        <TableHead>
          <TableRow>
            <TableCell sx={{ fontWeight: 600, bgcolor: 'grey.100' }}>#</TableCell>
            {data.headers.map((header) => (
              <TableCell
                key={header}
                sx={{
                  fontWeight: 600,
                  bgcolor: highlightColumns.includes(header) ? 'primary.50' : 'grey.100',
                  minWidth: 80,
                }}
              >
                {header}
              </TableCell>
            ))}
          </TableRow>
        </TableHead>
        <TableBody>
          {previewRows.map((row, idx) => (
            <TableRow key={idx} hover>
              <TableCell sx={{ color: 'text.secondary' }}>{idx + 1}</TableCell>
              {data.headers.map((header) => (
                <TableCell
                  key={header}
                  sx={{
                    bgcolor: highlightColumns.includes(header) ? 'primary.50' : undefined,
                  }}
                >
                  {row[header] !== null && row[header] !== undefined
                    ? String(row[header])
                    : '-'}
                </TableCell>
              ))}
            </TableRow>
          ))}
        </TableBody>
      </Table>
      {data.rows.length > maxRows && (
        <Box sx={{ p: 1, textAlign: 'center', bgcolor: 'grey.50' }}>
          <Typography variant="caption" color="text.secondary">
            Showing {maxRows} of {data.rows.length} rows
          </Typography>
        </Box>
      )}
    </TableContainer>
  );
}
