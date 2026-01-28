'use client';

import { useState, useCallback, useMemo, Suspense } from 'react';
import { useRouter, useSearchParams } from 'next/navigation';
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
  Grid2 as Grid,
  Divider,
  Card,
  CardContent,
  Skeleton,
} from '@mui/material';
import {
  Upload,
  ArrowBack,
  Science,
  TableChart,
  Settings,
  CheckCircle,
  Warning,
  Error as ErrorIcon,
} from '@mui/icons-material';
import { PageHeader } from '@/components/compounds/PageHeader';
import { SpreadsheetUpload, SpreadsheetData } from '@/components/compounds/SpreadsheetUpload';
import { useCompoundsApi, apiPost, apiUpload } from '@/lib/compounds/api';
import { useCompoundConfig } from '@/lib/compounds/config';
import { routes } from '@/lib/compounds/routes';

interface Protocol {
  id: string;
  name: string;
  import_type: string;
}

interface Target {
  id: string;
  name: string;
}

// Loading fallback for Suspense
function ImportPageSkeleton() {
  return (
    <Container maxWidth="xl" sx={{ py: 4 }}>
      <Skeleton variant="text" width={300} height={40} />
      <Skeleton variant="text" width={200} sx={{ mb: 3 }} />
      <Skeleton variant="rectangular" height={200} />
    </Container>
  );
}

// Main page component wrapped in Suspense
export default function ImportTableOfValuesPage() {
  return (
    <Suspense fallback={<ImportPageSkeleton />}>
      <ImportTableOfValuesContent />
    </Suspense>
  );
}

function ImportTableOfValuesContent() {
  const router = useRouter();
  const searchParams = useSearchParams();
  const api = useCompoundsApi();
  const { config: compoundConfig } = useCompoundConfig();

  // Get protocol from URL query param - if provided, it's locked
  const lockedProtocolId = searchParams.get('protocol');

  // Data state
  const [spreadsheetData, setSpreadsheetData] = useState<SpreadsheetData | null>(null);
  const [selectedProtocol, setSelectedProtocol] = useState<string | null>(lockedProtocolId);
  const [selectedTarget, setSelectedTarget] = useState<string | null>(null);
  const [compoundColumn, setCompoundColumn] = useState<string>('');
  const [kpiColumn, setKpiColumn] = useState<string>('');
  const [imageColumn, setImageColumn] = useState<string>('');
  const [kpiUnitOverride, setKpiUnitOverride] = useState<string>('');

  // Import state
  const [isImporting, setIsImporting] = useState(false);
  const [importError, setImportError] = useState<string | null>(null);

  // Fetch the locked protocol details if provided via URL
  const { data: lockedProtocol } = api.get<Protocol>(
    lockedProtocolId ? `protocols/${lockedProtocolId}/` : null
  );

  // Fetch protocols (filtered to table_of_values only) and targets
  const { data: protocolsData } = api.get<Protocol[]>(
    !lockedProtocolId ? 'protocols/?import_type=table_of_values' : null
  );
  const { data: targetsData } = api.get<Target[]>('targets/');

  const protocols = protocolsData || [];
  const targets = targetsData || [];

  // Protocol is locked if provided via URL
  const isProtocolLocked = !!lockedProtocolId;

  // Handle spreadsheet data loaded
  const handleDataLoaded = useCallback((data: SpreadsheetData) => {
    setSpreadsheetData(data);
    setImportError(null);

    // Auto-detect columns based on common names
    const headers = data.headers;
    const lowerHeaders = headers.map(h => h.toLowerCase());

    // Auto-detect compound column
    const compoundIdx = lowerHeaders.findIndex(h =>
      h.includes('compound') || h.includes('name') || h === 'id' || h === 'sample'
    );
    if (compoundIdx >= 0) {
      setCompoundColumn(headers[compoundIdx]);
    }

    // Auto-detect KPI column
    const kpiIdx = lowerHeaders.findIndex(h => h === 'kpi');
    if (kpiIdx >= 0) {
      setKpiColumn(headers[kpiIdx]);
    }

    // Auto-detect image column
    const imageIdx = lowerHeaders.findIndex(h =>
      h.includes('image') || h.includes('file') || h.includes('plot')
    );
    if (imageIdx >= 0) {
      setImageColumn(headers[imageIdx]);
    }
  }, []);

  // Parse unit from KPI field name (mirrors backend kpi_utils.py logic)
  const parseUnitFromFieldName = useCallback((fieldName: string): string | null => {
    if (!fieldName) return null;
    // Matches units in parentheses or square brackets
    const match = fieldName.match(/[\(\[]([nμµu]M|mM|pM|M|min|s|h|%|[μµu]L\/min\/mg|mL\/min\/kg|1e-6\s*cm\/s|cm\/s)[\)\]]/i);
    if (match) {
      // Normalize unicode mu to ASCII
      let unit = match[1];
      unit = unit.replace(/[μµ]/g, 'u');
      return unit;
    }
    return null;
  }, []);

  // Validate KPI column - all rows must have the same value
  const kpiValidation = useMemo(() => {
    if (!spreadsheetData || !kpiColumn) {
      return { valid: false, message: 'Select a KPI column', kpiValue: null, inferredUnit: null };
    }

    const values = new Set(
      spreadsheetData.rows
        .map(row => row[kpiColumn])
        .filter(v => v !== null && v !== undefined && v !== '')
    );

    if (values.size === 0) {
      return { valid: false, message: 'KPI column is empty', kpiValue: null, inferredUnit: null };
    }

    if (values.size > 1) {
      return {
        valid: false,
        message: `KPI column has multiple values: ${Array.from(values).join(', ')}`,
        kpiValue: null,
        inferredUnit: null,
      };
    }

    const kpiValue = Array.from(values)[0] as string;

    // Check if the KPI value is a column name
    if (!spreadsheetData.headers.includes(kpiValue)) {
      return {
        valid: false,
        message: `KPI value "${kpiValue}" is not a column name in the data`,
        kpiValue,
        inferredUnit: null,
      };
    }

    // Try to infer unit from KPI field name
    const inferredUnit = parseUnitFromFieldName(kpiValue);

    return { valid: true, message: `KPI: ${kpiValue}`, kpiValue, inferredUnit };
  }, [spreadsheetData, kpiColumn, parseUnitFromFieldName]);

  // Handle file clear
  const handleClearFile = useCallback(() => {
    setSpreadsheetData(null);
    setCompoundColumn('');
    setKpiColumn('');
    setImageColumn('');
    setKpiUnitOverride('');
    setImportError(null);
  }, []);

  // Handle import
  const handleImport = async () => {
    if (!spreadsheetData || !selectedProtocol || !compoundColumn || !kpiValidation.valid) {
      return;
    }

    setIsImporting(true);
    setImportError(null);

    try {
      // First, create the assay with the uploaded file
      const formData = new FormData();
      formData.append('protocol', selectedProtocol);
      if (selectedTarget) {
        formData.append('target', selectedTarget);
      }
      formData.append('comments', `Imported from ${spreadsheetData.fileName}`);
      // Include the original spreadsheet as the data_file
      if (spreadsheetData.originalFile) {
        formData.append('data_file', spreadsheetData.originalFile);
      }

      const assayResponse = await apiUpload<{ id: string }>('assays/', formData);

      const assayId = assayResponse.id;

      // Then import the table of values data
      // Use override unit or inferred unit
      const effectiveUnit = kpiUnitOverride || kpiValidation.inferredUnit;

      const importResponse = await apiPost<{
        status: string;
        created: number;
        errors_count: number;
        errors?: { row: number; error: string }[];
      }>(`assays/${assayId}/import_table_of_values/`, {
        compound_column: compoundColumn,
        kpi_column: kpiColumn,
        image_column: imageColumn || undefined,
        kpi_unit: effectiveUnit || undefined,
        data: spreadsheetData.rows,
      });

      if (importResponse.status === 'completed') {
        // Navigate to the new assay
        router.push(routes.assays.detail(assayId));
      } else {
        setImportError(
          `Import completed with ${importResponse.errors_count} errors. ` +
          `Created ${importResponse.created} data series.`
        );
      }
    } catch (error: any) {
      console.error('Import failed:', error);
      setImportError(error.message || 'Import failed. Please try again.');
    } finally {
      setIsImporting(false);
    }
  };

  // Check if we can import
  const canImport =
    spreadsheetData &&
    selectedProtocol &&
    compoundColumn &&
    kpiValidation.valid &&
    !isImporting;

  return (
    <Container maxWidth="xl" sx={{ py: 4 }}>
      <PageHeader
        breadcrumbs={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          ...(lockedProtocol
            ? [
                { label: 'Protocols', href: routes.assays.protocols() },
                { label: lockedProtocol.name, href: routes.assays.protocol(lockedProtocol.id) },
              ]
            : [{ label: 'Assays', href: routes.assays.list() }]),
          { label: 'Import Table of Values' },
        ]}
      />

      <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 3 }}>
        <Button
          component={Link}
          href={lockedProtocol ? routes.assays.protocol(lockedProtocol.id) : routes.assays.list()}
          startIcon={<ArrowBack />}
          size="small"
        >
          {lockedProtocol ? `Back to ${lockedProtocol.name}` : 'Back to Assays'}
        </Button>
        <Typography variant="h4" sx={{ flex: 1 }}>
          Import Table of Values
        </Typography>
      </Box>

      <Alert severity="info" sx={{ mb: 3 }}>
        Import pre-analyzed data from external analysis tools. The spreadsheet should have:
        <ul style={{ margin: '8px 0 0 0', paddingLeft: 20 }}>
          <li><strong>Compound column</strong> - Compound identifiers (e.g., {compoundConfig.compound_id_prefix}-00042)</li>
          <li><strong>KPI column</strong> - Contains the name of the column with the primary metric (e.g., all rows have "EC50")</li>
          <li><strong>Data columns</strong> - EC50, IC50, Hill, etc. as referenced by the KPI column</li>
          <li><strong>Image File column</strong> (optional) - Filenames for plot images to upload later</li>
        </ul>
      </Alert>

      {!spreadsheetData ? (
        // Step 1: File upload
        <SpreadsheetUpload
          title="Upload Table of Values File"
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
              {isProtocolLocked ? (
                // Locked protocol - show read-only display
                <Box sx={{ mt: 2, mb: 1 }}>
                  <Typography variant="caption" color="text.secondary">
                    Protocol
                  </Typography>
                  <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, mt: 0.5 }}>
                    <Chip
                      label={lockedProtocol?.name || 'Loading...'}
                      color="primary"
                      variant="outlined"
                    />
                    <Typography variant="caption" color="text.secondary">
                      (from protocol page)
                    </Typography>
                  </Box>
                </Box>
              ) : (
                // Unlocked - show selection dropdown
                <>
                  <Autocomplete
                    options={protocols}
                    getOptionLabel={(option) => option.name}
                    value={protocols.find((p) => p.id === selectedProtocol) || null}
                    onChange={(_, newValue) => setSelectedProtocol(newValue?.id || null)}
                    renderInput={(params) => (
                      <TextField
                        {...params}
                        label="Protocol *"
                        margin="normal"
                        fullWidth
                        helperText="Select a Table of Values protocol"
                        error={!selectedProtocol}
                      />
                    )}
                  />

                  {protocols.length === 0 && (
                    <Alert severity="warning" sx={{ mt: 1, mb: 2 }}>
                      No Table of Values protocols found.{' '}
                      <Link href={routes.assays.protocols()}>Create one first</Link>.
                    </Alert>
                  )}
                </>
              )}

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
                    helperText="Optional: target being tested"
                  />
                )}
              />

              <Divider sx={{ my: 2 }} />

              {/* Column configuration */}
              <Typography variant="subtitle2" gutterBottom>
                Column Mapping
              </Typography>

              {/* Compound column */}
              <FormControl fullWidth margin="normal" error={!compoundColumn}>
                <InputLabel>Compound Column *</InputLabel>
                <Select
                  value={compoundColumn}
                  label="Compound Column *"
                  onChange={(e) => setCompoundColumn(e.target.value)}
                >
                  {spreadsheetData.headers.map((header) => (
                    <MenuItem key={header} value={header}>
                      {header}
                    </MenuItem>
                  ))}
                </Select>
              </FormControl>

              {/* KPI column */}
              <FormControl fullWidth margin="normal" error={!kpiValidation.valid}>
                <InputLabel>KPI Column *</InputLabel>
                <Select
                  value={kpiColumn}
                  label="KPI Column *"
                  onChange={(e) => setKpiColumn(e.target.value)}
                >
                  {spreadsheetData.headers.map((header) => (
                    <MenuItem key={header} value={header}>
                      {header}
                    </MenuItem>
                  ))}
                </Select>
              </FormControl>

              {/* KPI validation status */}
              {kpiColumn && (
                <Box sx={{ mt: 1, mb: 2 }}>
                  {kpiValidation.valid ? (
                    <Chip
                      icon={<CheckCircle />}
                      label={kpiValidation.message}
                      color="success"
                      size="small"
                    />
                  ) : (
                    <Chip
                      icon={<Warning />}
                      label={kpiValidation.message}
                      color="error"
                      size="small"
                    />
                  )}
                </Box>
              )}

              {/* KPI Unit (only shown when KPI is valid) */}
              {kpiValidation.valid && (
                <TextField
                  fullWidth
                  label="KPI Unit"
                  value={kpiUnitOverride}
                  onChange={(e) => setKpiUnitOverride(e.target.value)}
                  size="small"
                  placeholder={kpiValidation.inferredUnit || 'No unit detected'}
                  helperText={
                    kpiValidation.inferredUnit
                      ? `Detected: ${kpiValidation.inferredUnit} (leave blank to use)`
                      : 'Enter unit if applicable (e.g., nM, uM, %)'
                  }
                  margin="normal"
                />
              )}

              {/* Image column */}
              <FormControl fullWidth margin="normal">
                <InputLabel>Image File Column</InputLabel>
                <Select
                  value={imageColumn}
                  label="Image File Column"
                  onChange={(e) => setImageColumn(e.target.value)}
                >
                  <MenuItem value="">
                    <em>None</em>
                  </MenuItem>
                  {spreadsheetData.headers.map((header) => (
                    <MenuItem key={header} value={header}>
                      {header}
                    </MenuItem>
                  ))}
                </Select>
              </FormControl>

              <Divider sx={{ my: 2 }} />

              {/* Summary */}
              <Card variant="outlined" sx={{ mb: 2 }}>
                <CardContent sx={{ py: 1.5 }}>
                  <Typography variant="subtitle2" gutterBottom>
                    Import Summary
                  </Typography>
                  <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: 0.5 }}>
                    <Chip
                      size="small"
                      label={`${spreadsheetData.rows.length} rows`}
                      color="primary"
                    />
                    {compoundColumn && (
                      <Chip
                        size="small"
                        icon={<CheckCircle />}
                        label="Compound column set"
                        color="success"
                      />
                    )}
                    {kpiValidation.valid && kpiValidation.kpiValue && (
                      <Chip
                        size="small"
                        icon={<CheckCircle />}
                        label={`KPI: ${kpiValidation.kpiValue}`}
                        color="success"
                      />
                    )}
                    {(kpiUnitOverride || kpiValidation.inferredUnit) && (
                      <Chip
                        size="small"
                        label={`Unit: ${kpiUnitOverride || kpiValidation.inferredUnit}`}
                        color="info"
                      />
                    )}
                    {imageColumn && (
                      <Chip
                        size="small"
                        label="Images configured"
                        color="info"
                      />
                    )}
                  </Box>
                </CardContent>
              </Card>

              {/* Error display */}
              {importError && (
                <Alert severity="error" sx={{ mb: 2 }}>
                  {importError}
                </Alert>
              )}

              {/* Actions */}
              <Button
                variant="contained"
                fullWidth
                disabled={!canImport}
                onClick={handleImport}
                startIcon={isImporting ? <CircularProgress size={20} /> : <Science />}
              >
                {isImporting ? 'Importing...' : 'Create Assay & Import Data'}
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
              </Box>

              <TableContainer sx={{ maxHeight: 600 }}>
                <Table size="small" stickyHeader>
                  <TableHead>
                    <TableRow>
                      <TableCell sx={{ fontWeight: 600, bgcolor: 'grey.100', width: 50 }}>
                        #
                      </TableCell>
                      {spreadsheetData.headers.map((header) => {
                        let bgColor = 'grey.100';
                        let chipLabel = '';
                        let chipColor: 'primary' | 'success' | 'info' | 'default' = 'default';

                        if (header === compoundColumn) {
                          bgColor = 'primary.50';
                          chipLabel = 'Compound';
                          chipColor = 'primary';
                        } else if (header === kpiColumn) {
                          bgColor = 'warning.50';
                          chipLabel = 'KPI';
                          chipColor = 'default';
                        } else if (kpiValidation.kpiValue && header === kpiValidation.kpiValue) {
                          bgColor = 'success.50';
                          chipLabel = 'Value';
                          chipColor = 'success';
                        } else if (header === imageColumn) {
                          bgColor = 'info.50';
                          chipLabel = 'Image';
                          chipColor = 'info';
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
                      <TableRow key={idx} hover>
                        <TableCell sx={{ color: 'text.secondary' }}>{idx + 1}</TableCell>
                        {spreadsheetData.headers.map((header) => {
                          let bgColor: string | undefined;

                          if (header === compoundColumn) {
                            bgColor = 'primary.50';
                          } else if (header === kpiColumn) {
                            bgColor = 'warning.50';
                          } else if (kpiValidation.kpiValue && header === kpiValidation.kpiValue) {
                            bgColor = 'success.50';
                          } else if (header === imageColumn) {
                            bgColor = 'info.50';
                          }

                          const value = row[header];
                          const displayValue =
                            value !== null && value !== undefined ? String(value) : '-';

                          return (
                            <TableCell
                              key={header}
                              sx={{
                                bgcolor: bgColor,
                                fontFamily:
                                  header === kpiValidation.kpiValue ? 'monospace' : undefined,
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
