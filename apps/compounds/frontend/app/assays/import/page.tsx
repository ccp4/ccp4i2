'use client';

import { useState, useCallback, useMemo } from 'react';
import { useRouter } from 'next/navigation';
import Link from 'next/link';
import {
  Container,
  Paper,
  Typography,
  Box,
  Button,
  Alert,
  CircularProgress,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Chip,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  Autocomplete,
  TextField,
  IconButton,
  Tooltip,
  Grid2 as Grid,
  Divider,
  Checkbox,
  FormControlLabel,
  Card,
  CardContent,
} from '@mui/material';
import {
  Upload,
  ArrowBack,
  Science,
  TableChart,
  Settings,
  CheckCircle,
  Warning,
  Visibility,
  VisibilityOff,
  Highlight,
} from '@mui/icons-material';
import { Breadcrumbs } from '@/components/compounds/Breadcrumbs';
import { SpreadsheetUpload, SpreadsheetData, SpreadsheetPreview, SpreadsheetRow } from '@/components/compounds/SpreadsheetUpload';
import { useCompoundsApi } from '@/lib/compounds/api';
import { routes } from '@/lib/compounds/routes';

interface Protocol {
  id: string;
  name: string;
  analysis_method: string;
}

interface Target {
  id: string;
  name: string;
}

interface ColumnConfig {
  type: 'compound_name' | 'data' | 'control_high' | 'control_low' | 'ignore';
  label?: string;
}

export default function ImportAssayPage() {
  const router = useRouter();
  const api = useCompoundsApi();

  // Data state
  const [spreadsheetData, setSpreadsheetData] = useState<SpreadsheetData | null>(null);
  const [selectedProtocol, setSelectedProtocol] = useState<string | null>(null);
  const [selectedTarget, setSelectedTarget] = useState<string | null>(null);
  const [columnConfig, setColumnConfig] = useState<Record<string, ColumnConfig>>({});
  const [showRawPreview, setShowRawPreview] = useState(true);
  const [highlightedRow, setHighlightedRow] = useState<number | null>(null);

  // Fetch protocols and targets
  const { data: protocolsData } = api.get<Protocol[]>('protocols/');
  const { data: targetsData } = api.get<Target[]>('targets/');

  const protocols = protocolsData || [];
  const targets = targetsData || [];

  // Handle spreadsheet data loaded
  const handleDataLoaded = useCallback((data: SpreadsheetData) => {
    setSpreadsheetData(data);

    // Auto-detect column types based on header names
    const autoConfig: Record<string, ColumnConfig> = {};
    data.headers.forEach((header) => {
      const lowerHeader = header.toLowerCase();
      if (
        lowerHeader.includes('compound') ||
        lowerHeader.includes('name') ||
        lowerHeader.includes('id') ||
        lowerHeader === 'row' ||
        lowerHeader === 'sample'
      ) {
        autoConfig[header] = { type: 'compound_name', label: header };
      } else if (
        lowerHeader.includes('high') ||
        lowerHeader.includes('pos') ||
        lowerHeader === 'h'
      ) {
        autoConfig[header] = { type: 'control_high', label: header };
      } else if (
        lowerHeader.includes('low') ||
        lowerHeader.includes('neg') ||
        lowerHeader === 'l'
      ) {
        autoConfig[header] = { type: 'control_low', label: header };
      } else {
        // Assume numeric columns are data
        const firstValue = data.rows[0]?.[header];
        if (firstValue !== null && !isNaN(Number(firstValue))) {
          autoConfig[header] = { type: 'data', label: header };
        } else {
          autoConfig[header] = { type: 'ignore', label: header };
        }
      }
    });
    setColumnConfig(autoConfig);
  }, []);

  // Handle column configuration change
  const handleColumnTypeChange = useCallback(
    (column: string, type: ColumnConfig['type']) => {
      setColumnConfig((prev) => ({
        ...prev,
        [column]: { ...prev[column], type },
      }));
    },
    []
  );

  // Get columns by type
  const columnsByType = useMemo(() => {
    const result: Record<string, string[]> = {
      compound_name: [],
      data: [],
      control_high: [],
      control_low: [],
      ignore: [],
    };

    for (const [column, config] of Object.entries(columnConfig)) {
      if (result[config.type]) {
        result[config.type].push(column);
      }
    }

    return result;
  }, [columnConfig]);

  // Summary stats
  const stats = useMemo(() => {
    if (!spreadsheetData) return null;

    return {
      totalRows: spreadsheetData.rows.length,
      dataColumns: columnsByType.data.length,
      hasCompoundColumn: columnsByType.compound_name.length > 0,
      hasHighControl: columnsByType.control_high.length > 0,
      hasLowControl: columnsByType.control_low.length > 0,
    };
  }, [spreadsheetData, columnsByType]);

  // Handle file clear
  const handleClearFile = useCallback(() => {
    setSpreadsheetData(null);
    setColumnConfig({});
    setSelectedProtocol(null);
    setSelectedTarget(null);
    setHighlightedRow(null);
  }, []);

  return (
    <Container maxWidth="xl" sx={{ py: 4 }}>
      <Breadcrumbs
        items={[
          { label: 'Assays', href: routes.assays.list() },
          { label: 'Import Data' },
        ]}
      />

      <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 3 }}>
        <Button
          component={Link}
          href={routes.assays.list()}
          startIcon={<ArrowBack />}
          size="small"
        >
          Back to Assays
        </Button>
        <Typography variant="h4" sx={{ flex: 1 }}>
          Import Assay Data
        </Typography>
      </Box>

      {!spreadsheetData ? (
        // Step 1: File upload
        <SpreadsheetUpload
          title="Upload Assay Data File"
          onDataLoaded={handleDataLoaded}
          showColumnMapping={false}
          previewRows={10}
        />
      ) : (
        // Step 2: Configure and preview
        <Grid container spacing={3}>
          {/* Left panel: Configuration */}
          <Grid size={{ xs: 12, md: 4 }}>
            <Paper sx={{ p: 3, position: 'sticky', top: 16 }}>
              <Typography variant="h6" gutterBottom sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                <Settings color="primary" />
                Configuration
              </Typography>

              {/* File info */}
              <Box sx={{ mb: 3, p: 2, bgcolor: 'success.50', borderRadius: 1 }}>
                <Typography variant="body2" fontWeight={600}>
                  {spreadsheetData.fileName}
                </Typography>
                <Typography variant="caption" color="text.secondary">
                  {spreadsheetData.rows.length} rows, {spreadsheetData.headers.length} columns
                  {spreadsheetData.sheetName && ` (Sheet: ${spreadsheetData.sheetName})`}
                </Typography>
                <Button size="small" onClick={handleClearFile} sx={{ mt: 1 }}>
                  Change File
                </Button>
              </Box>

              {/* Protocol selection */}
              <Autocomplete
                options={protocols}
                getOptionLabel={(option) => option.name}
                value={protocols.find((p) => p.id === selectedProtocol) || null}
                onChange={(_, newValue) => setSelectedProtocol(newValue?.id || null)}
                renderInput={(params) => (
                  <TextField
                    {...params}
                    label="Protocol"
                    margin="normal"
                    fullWidth
                    helperText="Select the assay protocol for this data"
                  />
                )}
              />

              {/* Target selection */}
              <Autocomplete
                options={targets}
                getOptionLabel={(option) => option.name}
                value={targets.find((t) => t.id === selectedTarget) || null}
                onChange={(_, newValue) => setSelectedTarget(newValue?.id || null)}
                renderInput={(params) => (
                  <TextField
                    {...params}
                    label="Target"
                    margin="normal"
                    fullWidth
                    helperText="Select the target being tested"
                  />
                )}
              />

              <Divider sx={{ my: 2 }} />

              {/* Column configuration */}
              <Typography variant="subtitle2" gutterBottom>
                Column Types
              </Typography>
              <Typography variant="caption" color="text.secondary" paragraph>
                Configure how each column should be interpreted
              </Typography>

              <Box sx={{ maxHeight: 300, overflow: 'auto' }}>
                {spreadsheetData.headers.map((header) => (
                  <FormControl key={header} fullWidth size="small" sx={{ mb: 1 }}>
                    <InputLabel>{header}</InputLabel>
                    <Select
                      value={columnConfig[header]?.type || 'ignore'}
                      label={header}
                      onChange={(e) =>
                        handleColumnTypeChange(header, e.target.value as ColumnConfig['type'])
                      }
                    >
                      <MenuItem value="compound_name">
                        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                          <Chip size="small" label="ID" color="primary" />
                          Compound Name/ID
                        </Box>
                      </MenuItem>
                      <MenuItem value="data">
                        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                          <Chip size="small" label="Data" color="success" />
                          Response Data
                        </Box>
                      </MenuItem>
                      <MenuItem value="control_high">
                        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                          <Chip size="small" label="H" color="warning" />
                          High Control
                        </Box>
                      </MenuItem>
                      <MenuItem value="control_low">
                        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                          <Chip size="small" label="L" color="info" />
                          Low Control
                        </Box>
                      </MenuItem>
                      <MenuItem value="ignore">
                        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                          <Chip size="small" label="—" />
                          Ignore
                        </Box>
                      </MenuItem>
                    </Select>
                  </FormControl>
                ))}
              </Box>

              <Divider sx={{ my: 2 }} />

              {/* Summary */}
              {stats && (
                <Card variant="outlined" sx={{ mb: 2 }}>
                  <CardContent sx={{ py: 1.5 }}>
                    <Typography variant="subtitle2" gutterBottom>
                      Configuration Summary
                    </Typography>
                    <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: 0.5 }}>
                      {stats.hasCompoundColumn ? (
                        <Chip
                          size="small"
                          icon={<CheckCircle />}
                          label="Compound column"
                          color="success"
                        />
                      ) : (
                        <Chip
                          size="small"
                          icon={<Warning />}
                          label="No compound column"
                          color="warning"
                        />
                      )}
                      <Chip
                        size="small"
                        label={`${stats.dataColumns} data columns`}
                        color={stats.dataColumns > 0 ? 'success' : 'error'}
                      />
                      {stats.hasHighControl && (
                        <Chip size="small" label="High control" color="info" />
                      )}
                      {stats.hasLowControl && (
                        <Chip size="small" label="Low control" color="info" />
                      )}
                    </Box>
                  </CardContent>
                </Card>
              )}

              <Alert severity="info" sx={{ mb: 2 }}>
                <Typography variant="caption">
                  Further configuration (compound matching, control selection, data extraction)
                  will be available in the next step.
                </Typography>
              </Alert>

              {/* Actions */}
              <Button
                variant="contained"
                fullWidth
                disabled={!selectedProtocol || columnsByType.data.length === 0}
                startIcon={<Science />}
              >
                Continue to Data Extraction
              </Button>
            </Paper>
          </Grid>

          {/* Right panel: Preview */}
          <Grid size={{ xs: 12, md: 8 }}>
            <Paper sx={{ p: 3 }}>
              <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 2 }}>
                <Typography variant="h6" sx={{ flex: 1, display: 'flex', alignItems: 'center', gap: 1 }}>
                  <TableChart color="primary" />
                  Data Preview
                </Typography>
                <FormControlLabel
                  control={
                    <Checkbox
                      checked={showRawPreview}
                      onChange={(e) => setShowRawPreview(e.target.checked)}
                    />
                  }
                  label="Show all columns"
                />
              </Box>

              <TableContainer sx={{ maxHeight: 600 }}>
                <Table size="small" stickyHeader>
                  <TableHead>
                    <TableRow>
                      <TableCell sx={{ fontWeight: 600, bgcolor: 'grey.100', width: 50 }}>
                        #
                      </TableCell>
                      {spreadsheetData.headers.map((header) => {
                        const config = columnConfig[header];
                        if (!showRawPreview && config?.type === 'ignore') {
                          return null;
                        }

                        let bgColor = 'grey.100';
                        let chipColor: any = 'default';
                        let chipLabel = '';

                        switch (config?.type) {
                          case 'compound_name':
                            bgColor = 'primary.50';
                            chipColor = 'primary';
                            chipLabel = 'ID';
                            break;
                          case 'data':
                            bgColor = 'success.50';
                            chipColor = 'success';
                            chipLabel = 'Data';
                            break;
                          case 'control_high':
                            bgColor = 'warning.50';
                            chipColor = 'warning';
                            chipLabel = 'H';
                            break;
                          case 'control_low':
                            bgColor = 'info.50';
                            chipColor = 'info';
                            chipLabel = 'L';
                            break;
                        }

                        return (
                          <TableCell
                            key={header}
                            sx={{ fontWeight: 600, bgcolor: bgColor, minWidth: 80 }}
                          >
                            <Box sx={{ display: 'flex', alignItems: 'center', gap: 0.5 }}>
                              {chipLabel && (
                                <Chip label={chipLabel} size="small" color={chipColor} />
                              )}
                              {header}
                            </Box>
                          </TableCell>
                        );
                      })}
                    </TableRow>
                  </TableHead>
                  <TableBody>
                    {spreadsheetData.rows.slice(0, 50).map((row, idx) => (
                      <TableRow
                        key={idx}
                        hover
                        selected={highlightedRow === idx}
                        onClick={() => setHighlightedRow(idx === highlightedRow ? null : idx)}
                        sx={{ cursor: 'pointer' }}
                      >
                        <TableCell sx={{ color: 'text.secondary' }}>{idx + 1}</TableCell>
                        {spreadsheetData.headers.map((header) => {
                          const config = columnConfig[header];
                          if (!showRawPreview && config?.type === 'ignore') {
                            return null;
                          }

                          let bgColor: string | undefined;
                          switch (config?.type) {
                            case 'compound_name':
                              bgColor = 'primary.50';
                              break;
                            case 'data':
                              bgColor = 'success.50';
                              break;
                            case 'control_high':
                              bgColor = 'warning.50';
                              break;
                            case 'control_low':
                              bgColor = 'info.50';
                              break;
                          }

                          const value = row[header];
                          const displayValue =
                            value !== null && value !== undefined
                              ? String(value)
                              : '-';

                          return (
                            <TableCell
                              key={header}
                              sx={{
                                bgcolor: bgColor,
                                fontFamily:
                                  config?.type === 'data' ? 'monospace' : undefined,
                              }}
                            >
                              {displayValue}
                            </TableCell>
                          );
                        })}
                      </TableRow>
                    ))}
                  </TableBody>
                </Table>
              </TableContainer>

              {spreadsheetData.rows.length > 50 && (
                <Typography
                  variant="caption"
                  color="text.secondary"
                  sx={{ display: 'block', mt: 1, textAlign: 'center' }}
                >
                  Showing first 50 of {spreadsheetData.rows.length} rows
                </Typography>
              )}
            </Paper>
          </Grid>
        </Grid>
      )}
    </Container>
  );
}
