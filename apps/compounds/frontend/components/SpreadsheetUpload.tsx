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

import { useState, useCallback, useMemo } from 'react';
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

export interface SpreadsheetData {
  headers: string[];
  rows: SpreadsheetRow[];
  fileName: string;
  sheetName: string;
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

  const parseFile = useCallback(async (file: File, sheetName?: string) => {
    setLoading(true);
    setError(null);

    try {
      const buffer = await file.arrayBuffer();
      let sheetNames: string[];
      let targetSheet: string;
      let headers: string[];
      let jsonData: SpreadsheetRow[];

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
        jsonData = [];
        for (let i = 2; i <= worksheet.rowCount; i++) {
          const row = worksheet.getRow(i);
          const obj: SpreadsheetRow = {};
          let hasData = false;
          headers.forEach((h, idx) => {
            const cell = row.getCell(idx + 1);
            if (cell.value != null) {
              obj[h] = String(cell.value);
              hasData = true;
            } else {
              obj[h] = null;
            }
          });
          if (hasData) jsonData.push(obj);
        }
      }

      setSheets(sheetNames);
      setSelectedSheet(targetSheet);

      if (jsonData.length === 0) {
        throw new Error('Spreadsheet appears to be empty');
      }

      const spreadsheetData: SpreadsheetData = {
        headers,
        rows: jsonData,
        fileName: file.name,
        sheetName: targetSheet,
      };

      setData(spreadsheetData);

      // Auto-map columns if field mappings provided
      if (fieldMappings.length > 0) {
        const autoMapping: Record<string, string> = {};
        fieldMappings.forEach((field) => {
          // Try exact match first
          const exactMatch = headers.find(
            (h) => h.toLowerCase() === field.field.toLowerCase()
          );
          if (exactMatch) {
            autoMapping[field.field] = exactMatch;
            return;
          }
          // Try partial match
          const partialMatch = headers.find(
            (h) =>
              h.toLowerCase().includes(field.field.toLowerCase()) ||
              field.field.toLowerCase().includes(h.toLowerCase())
          );
          if (partialMatch) {
            autoMapping[field.field] = partialMatch;
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
        // File upload area
        <Box
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
          component="label"
        >
          <input
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
                {previewData?.rows.map((row, idx) => (
                  <TableRow key={idx} hover>
                    <TableCell sx={{ color: 'text.secondary' }}>{idx + 1}</TableCell>
                    {previewData.headers.map((header) => (
                      <TableCell key={header}>
                        {row[header] !== null && row[header] !== undefined
                          ? String(row[header])
                          : '-'}
                      </TableCell>
                    ))}
                  </TableRow>
                ))}
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
