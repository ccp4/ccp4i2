'use client';

import { useState, useCallback, Suspense } from 'react';
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
  Grid2 as Grid,
  Divider,
  Card,
  CardContent,
  Skeleton,
  FormControlLabel,
  Checkbox,
  TextField,
  Tooltip,
  IconButton,
  Collapse,
  Autocomplete,
} from '@mui/material';
import {
  Upload,
  ArrowBack,
  Science,
  CheckCircle,
  Warning,
  Error as ErrorIcon,
  Cancel,
  ExpandMore,
  ExpandLess,
  CloudUpload,
  Link as LinkIcon,
  LinkOff,
} from '@mui/icons-material';
import { PageHeader } from '@/components/compounds/PageHeader';
import { routes } from '@/lib/compounds/routes';
import { apiUpload, useCompoundsApi } from '@/lib/compounds/api';

interface Target {
  id: string;
  name: string;
}

interface ParsedResult {
  compound_id: string;
  compound_matched: boolean;
  compound_reg_number: string | null;
  compound_db_id: string | null;
  species: string | null;
  assay_format: string | null;
  is_control: boolean;
  results: Record<string, any>;
  kpi_value: number | null;
  flags: string[];
}

interface ParseError {
  message: string;
  severity: string;
  row: number | null;
  column: string | null;
}

interface ParsePreviewResponse {
  status: string;
  filename: string;
  parser: {
    vendor: string;
    assay_type: string;
    protocol_slug: string;
    kpi_field: string;
  };
  metadata: Record<string, any>;
  results: ParsedResult[];
  errors: ParseError[];
  summary: {
    total: number;
    matched: number;
    unmatched: number;
    controls: number;
  };
}

interface ImportResponse {
  status: string;
  assay_id: string;
  created: number;
  skipped: number;
  errors_count: number;
  errors?: Array<{ compound_id: string; error: string }>;
}

// Loading skeleton
function ImportPageSkeleton() {
  return (
    <Container maxWidth="xl" sx={{ py: 4 }}>
      <Skeleton variant="text" width={300} height={40} />
      <Skeleton variant="text" width={200} sx={{ mb: 3 }} />
      <Skeleton variant="rectangular" height={200} />
    </Container>
  );
}

export default function ImportADMEPage() {
  return (
    <Suspense fallback={<ImportPageSkeleton />}>
      <ImportADMEContent />
    </Suspense>
  );
}

function ImportADMEContent() {
  const router = useRouter();
  const api = useCompoundsApi();

  // Fetch targets for the dropdown
  const { data: targetsData } = api.get<Target[]>('targets/');
  const targets = targetsData || [];

  // File state
  const [selectedFile, setSelectedFile] = useState<File | null>(null);
  const [isParsing, setIsParsing] = useState(false);
  const [parseError, setParseError] = useState<string | null>(null);

  // Preview state
  const [preview, setPreview] = useState<ParsePreviewResponse | null>(null);

  // Import options
  const [skipUnmatched, setSkipUnmatched] = useState(false);
  const [importComments, setImportComments] = useState('');
  const [selectedTarget, setSelectedTarget] = useState<string | null>(null);

  // Import state
  const [isImporting, setIsImporting] = useState(false);
  const [importError, setImportError] = useState<string | null>(null);

  // UI state
  const [showErrors, setShowErrors] = useState(false);
  const [showDetails, setShowDetails] = useState<string | null>(null);

  // Handle file selection
  const handleFileSelect = useCallback(async (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (!file) return;

    setSelectedFile(file);
    setParseError(null);
    setPreview(null);
    setIsParsing(true);

    try {
      const formData = new FormData();
      formData.append('file', file);

      const response = await apiUpload<ParsePreviewResponse>('assays/parse_adme_preview/', formData);

      if (response.status === 'error') {
        setParseError((response as any).error || 'Failed to parse file');
      } else {
        setPreview(response);
      }
    } catch (error: any) {
      console.error('Parse failed:', error);
      setParseError(error.message || 'Failed to parse file');
    } finally {
      setIsParsing(false);
    }
  }, []);

  // Handle file drop
  const handleDrop = useCallback((event: React.DragEvent<HTMLDivElement>) => {
    event.preventDefault();
    const file = event.dataTransfer.files?.[0];
    if (file) {
      // Trigger file select handler
      const dataTransfer = new DataTransfer();
      dataTransfer.items.add(file);
      const input = document.getElementById('adme-file-input') as HTMLInputElement;
      if (input) {
        input.files = dataTransfer.files;
        input.dispatchEvent(new Event('change', { bubbles: true }));
      }
    }
  }, []);

  // Handle clear
  const handleClear = useCallback(() => {
    setSelectedFile(null);
    setPreview(null);
    setParseError(null);
    setImportError(null);
    setSkipUnmatched(false);
    setImportComments('');
    setSelectedTarget(null);
  }, []);

  // Handle import
  const handleImport = useCallback(async () => {
    if (!preview) return;

    setIsImporting(true);
    setImportError(null);

    try {
      // Use FormData to include the original file for download later
      const formData = new FormData();
      formData.append('filename', preview.filename);
      formData.append('parser_slug', preview.parser.protocol_slug);
      formData.append('results', JSON.stringify(preview.results.filter(r => !r.is_control)));
      formData.append('skip_unmatched', String(skipUnmatched));
      if (importComments) {
        formData.append('comments', importComments);
      }
      if (selectedTarget) {
        formData.append('target', selectedTarget);
      }
      // Include the original file so it can be downloaded later
      if (selectedFile) {
        formData.append('file', selectedFile);
      }

      const response = await apiUpload<ImportResponse>('assays/import_adme/', formData);

      if (response.status === 'completed') {
        // Navigate to the new assay
        router.push(routes.assays.detail(response.assay_id));
      } else {
        setImportError(`Import completed with ${response.errors_count} errors`);
      }
    } catch (error: any) {
      console.error('Import failed:', error);
      setImportError(error.message || 'Import failed');
    } finally {
      setIsImporting(false);
    }
  }, [preview, skipUnmatched, importComments, selectedTarget, selectedFile, router]);

  // Calculate import count
  const importCount = preview
    ? preview.results.filter(r => !r.is_control && (r.compound_matched || !skipUnmatched)).length
    : 0;

  return (
    <Container maxWidth="xl" sx={{ py: 4 }}>
      <PageHeader
        breadcrumbs={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Assays', href: routes.assays.list() },
          { label: 'Import ADME Data' },
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
          Import ADME Data
        </Typography>
      </Box>

      <Alert severity="info" sx={{ mb: 3 }}>
        Import pre-analyzed ADME data from vendor Excel files. Supported formats:
        <Box component="ul" sx={{ m: 0, mt: 1, pl: 2 }}>
          <li><strong>NCU LM</strong> - Liver Microsome Stability (ADME-NCU-LM-YYYYMMDD.xlsx)</li>
          <li><strong>NCU BS</strong> - Blood/Serum Stability (ADME-NCU-BS-YYYYMMDD.xlsx)</li>
          <li><strong>NCU GSH</strong> - GSH Stability (ADME-NCU-GSH Stability-YYYYMMDD.xlsx)</li>
          <li><strong>NCU Caco-2</strong> - Permeability (ADME-NCU-Caco-2 Permeability-YYYYMMDD.xlsx)</li>
        </Box>
      </Alert>

      {!preview ? (
        // Step 1: File upload
        <Paper
          sx={{
            p: 4,
            textAlign: 'center',
            border: '2px dashed',
            borderColor: parseError ? 'error.main' : 'divider',
            bgcolor: 'background.default',
            cursor: 'pointer',
            transition: 'border-color 0.2s',
            '&:hover': {
              borderColor: 'primary.main',
            },
          }}
          onDrop={handleDrop}
          onDragOver={(e) => e.preventDefault()}
        >
          <input
            id="adme-file-input"
            type="file"
            accept=".xlsx,.xls"
            onChange={handleFileSelect}
            style={{ display: 'none' }}
          />

          {isParsing ? (
            <Box sx={{ py: 4 }}>
              <CircularProgress size={48} sx={{ mb: 2 }} />
              <Typography variant="h6">Parsing {selectedFile?.name}...</Typography>
              <Typography variant="body2" color="text.secondary">
                Detecting assay type and extracting data
              </Typography>
            </Box>
          ) : (
            <label htmlFor="adme-file-input" style={{ cursor: 'pointer', display: 'block' }}>
              <CloudUpload sx={{ fontSize: 64, color: 'action.active', mb: 2 }} />
              <Typography variant="h6" gutterBottom>
                Drop ADME Excel file here or click to browse
              </Typography>
              <Typography variant="body2" color="text.secondary">
                File type will be auto-detected from filename
              </Typography>
            </label>
          )}

          {parseError && (
            <Alert severity="error" sx={{ mt: 2, textAlign: 'left' }}>
              {parseError}
            </Alert>
          )}
        </Paper>
      ) : (
        // Step 2: Preview and import
        <Grid container spacing={3}>
          {/* Left panel: Summary and options */}
          <Grid size={{ xs: 12, md: 4 }}>
            <Paper sx={{ p: 3, position: 'sticky', top: 16 }}>
              {/* File info */}
              <Box sx={{ mb: 3, p: 2, bgcolor: 'success.50', borderRadius: 1 }}>
                <Typography variant="body2" fontWeight={600}>
                  {preview.filename}
                </Typography>
                <Chip
                  label={`${preview.parser.vendor} ${preview.parser.assay_type.replace(/_/g, ' ')}`}
                  size="small"
                  color="primary"
                  sx={{ mt: 1 }}
                />
                <Button size="small" onClick={handleClear} sx={{ mt: 1, ml: 1 }}>
                  Change File
                </Button>
              </Box>

              <Divider sx={{ my: 2 }} />

              {/* Summary */}
              <Typography variant="h6" gutterBottom>
                Summary
              </Typography>

              <Box sx={{ display: 'flex', flexDirection: 'column', gap: 1, mb: 2 }}>
                <Box sx={{ display: 'flex', justifyContent: 'space-between' }}>
                  <Typography variant="body2">Total compounds:</Typography>
                  <Typography variant="body2" fontWeight={600}>
                    {preview.summary.total - preview.summary.controls}
                  </Typography>
                </Box>
                <Box sx={{ display: 'flex', justifyContent: 'space-between' }}>
                  <Typography variant="body2" sx={{ display: 'flex', alignItems: 'center', gap: 0.5 }}>
                    <LinkIcon fontSize="small" color="success" /> Matched:
                  </Typography>
                  <Typography variant="body2" fontWeight={600} color="success.main">
                    {preview.summary.matched}
                  </Typography>
                </Box>
                <Box sx={{ display: 'flex', justifyContent: 'space-between' }}>
                  <Typography variant="body2" sx={{ display: 'flex', alignItems: 'center', gap: 0.5 }}>
                    <LinkOff fontSize="small" color="warning" /> Unmatched:
                  </Typography>
                  <Typography variant="body2" fontWeight={600} color="warning.main">
                    {preview.summary.unmatched}
                  </Typography>
                </Box>
                {preview.summary.controls > 0 && (
                  <Box sx={{ display: 'flex', justifyContent: 'space-between' }}>
                    <Typography variant="body2" color="text.secondary">
                      Controls (skipped):
                    </Typography>
                    <Typography variant="body2" color="text.secondary">
                      {preview.summary.controls}
                    </Typography>
                  </Box>
                )}
              </Box>

              {/* Errors */}
              {preview.errors.length > 0 && (
                <Box sx={{ mb: 2 }}>
                  <Button
                    size="small"
                    onClick={() => setShowErrors(!showErrors)}
                    startIcon={showErrors ? <ExpandLess /> : <ExpandMore />}
                    color="warning"
                  >
                    {preview.errors.length} warnings/errors
                  </Button>
                  <Collapse in={showErrors}>
                    <Box sx={{ mt: 1, maxHeight: 150, overflow: 'auto' }}>
                      {preview.errors.map((err, idx) => (
                        <Alert
                          key={idx}
                          severity={err.severity === 'error' ? 'error' : 'warning'}
                          sx={{ mb: 0.5, py: 0 }}
                        >
                          <Typography variant="caption">{err.message}</Typography>
                        </Alert>
                      ))}
                    </Box>
                  </Collapse>
                </Box>
              )}

              <Divider sx={{ my: 2 }} />

              {/* Options */}
              <Typography variant="h6" gutterBottom>
                Import Options
              </Typography>

              {/* Target selection */}
              <Autocomplete
                options={targets}
                getOptionLabel={(option) => option.name}
                value={targets.find((t) => t.id === selectedTarget) || null}
                onChange={(_, newValue) => setSelectedTarget(newValue?.id || null)}
                size="small"
                renderInput={(params) => (
                  <TextField
                    {...params}
                    label="Target"
                    helperText="Optional: target being tested"
                  />
                )}
                sx={{ mb: 2 }}
              />

              <FormControlLabel
                control={
                  <Checkbox
                    checked={skipUnmatched}
                    onChange={(e) => setSkipUnmatched(e.target.checked)}
                  />
                }
                label={
                  <Typography variant="body2">
                    Skip unmatched compounds ({preview.summary.unmatched})
                  </Typography>
                }
              />

              <TextField
                fullWidth
                multiline
                rows={2}
                label="Import notes (optional)"
                value={importComments}
                onChange={(e) => setImportComments(e.target.value)}
                size="small"
                sx={{ mt: 2 }}
              />

              <Divider sx={{ my: 2 }} />

              {/* Import button */}
              <Card variant="outlined" sx={{ mb: 2 }}>
                <CardContent sx={{ py: 1.5 }}>
                  <Typography variant="subtitle2">
                    Will import: <strong>{importCount}</strong> results
                  </Typography>
                  <Typography variant="caption" color="text.secondary">
                    KPI: {preview.parser.kpi_field.replace(/_/g, ' ')}
                  </Typography>
                </CardContent>
              </Card>

              {importError && (
                <Alert severity="error" sx={{ mb: 2 }}>
                  {importError}
                </Alert>
              )}

              <Button
                variant="contained"
                fullWidth
                disabled={isImporting || importCount === 0}
                onClick={handleImport}
                startIcon={isImporting ? <CircularProgress size={20} /> : <Science />}
              >
                {isImporting ? 'Importing...' : `Import ${importCount} Results`}
              </Button>
            </Paper>
          </Grid>

          {/* Right panel: Results preview */}
          <Grid size={{ xs: 12, md: 8 }}>
            <Paper sx={{ p: 3 }}>
              <Typography variant="h6" gutterBottom>
                Parsed Results
              </Typography>

              <TableContainer sx={{ maxHeight: 600 }}>
                <Table size="small" stickyHeader>
                  <TableHead>
                    <TableRow>
                      <TableCell sx={{ fontWeight: 600, bgcolor: 'grey.100' }}>
                        Compound
                      </TableCell>
                      <TableCell sx={{ fontWeight: 600, bgcolor: 'grey.100' }}>
                        Match
                      </TableCell>
                      {preview.results.some(r => r.species) && (
                        <TableCell sx={{ fontWeight: 600, bgcolor: 'grey.100' }}>
                          Species
                        </TableCell>
                      )}
                      <TableCell sx={{ fontWeight: 600, bgcolor: 'grey.100' }}>
                        {preview.parser.kpi_field.replace(/_/g, ' ')}
                      </TableCell>
                      <TableCell sx={{ fontWeight: 600, bgcolor: 'grey.100' }}>
                        Flags
                      </TableCell>
                      <TableCell sx={{ fontWeight: 600, bgcolor: 'grey.100', width: 50 }} />
                    </TableRow>
                  </TableHead>
                  <TableBody>
                    {preview.results
                      .filter(r => !r.is_control)
                      .map((result, idx) => (
                        <>
                          <TableRow
                            key={idx}
                            hover
                            sx={{
                              bgcolor: !result.compound_matched && skipUnmatched ? 'action.disabledBackground' : undefined,
                              opacity: !result.compound_matched && skipUnmatched ? 0.5 : 1,
                            }}
                          >
                            <TableCell>
                              <Typography variant="body2" fontWeight={500}>
                                {result.compound_id}
                              </Typography>
                              {result.compound_reg_number && (
                                <Typography variant="caption" color="text.secondary">
                                  {result.compound_reg_number}
                                </Typography>
                              )}
                            </TableCell>
                            <TableCell>
                              {result.compound_matched ? (
                                <Tooltip title="Matched to registry">
                                  <CheckCircle color="success" fontSize="small" />
                                </Tooltip>
                              ) : (
                                <Tooltip title="No match in registry">
                                  <Warning color="warning" fontSize="small" />
                                </Tooltip>
                              )}
                            </TableCell>
                            {preview.results.some(r => r.species) && (
                              <TableCell>
                                <Typography variant="body2">{result.species || '-'}</Typography>
                              </TableCell>
                            )}
                            <TableCell>
                              <Typography variant="body2" fontFamily="monospace">
                                {result.kpi_value != null ? result.kpi_value.toFixed(2) : '-'}
                              </Typography>
                            </TableCell>
                            <TableCell>
                              {result.flags.length > 0 ? (
                                <Box sx={{ display: 'flex', gap: 0.5, flexWrap: 'wrap' }}>
                                  {result.flags.slice(0, 2).map((flag, i) => (
                                    <Chip
                                      key={i}
                                      label={flag}
                                      size="small"
                                      color={flag.includes('unstable') || flag.includes('high') ? 'warning' : 'default'}
                                    />
                                  ))}
                                  {result.flags.length > 2 && (
                                    <Chip label={`+${result.flags.length - 2}`} size="small" />
                                  )}
                                </Box>
                              ) : (
                                <Typography variant="body2" color="text.secondary">-</Typography>
                              )}
                            </TableCell>
                            <TableCell>
                              <IconButton
                                size="small"
                                onClick={() => setShowDetails(showDetails === result.compound_id ? null : result.compound_id)}
                              >
                                {showDetails === result.compound_id ? <ExpandLess /> : <ExpandMore />}
                              </IconButton>
                            </TableCell>
                          </TableRow>
                          <TableRow>
                            <TableCell colSpan={6} sx={{ py: 0 }}>
                              <Collapse in={showDetails === result.compound_id}>
                                <Box sx={{ p: 2, bgcolor: 'grey.50' }}>
                                  <Typography variant="subtitle2" gutterBottom>
                                    Full Results
                                  </Typography>
                                  <Grid container spacing={1}>
                                    {Object.entries(result.results)
                                      .filter(([key]) => !['flags', 'source', 'time_course', 'monolayer_integrity'].includes(key))
                                      .map(([key, value]) => (
                                        <Grid size={{ xs: 6, sm: 4, md: 3 }} key={key}>
                                          <Typography variant="caption" color="text.secondary">
                                            {key.replace(/_/g, ' ')}
                                          </Typography>
                                          <Typography variant="body2" fontFamily="monospace">
                                            {typeof value === 'number' ? value.toFixed(3) : String(value ?? '-')}
                                          </Typography>
                                        </Grid>
                                      ))}
                                  </Grid>
                                </Box>
                              </Collapse>
                            </TableCell>
                          </TableRow>
                        </>
                      ))}
                  </TableBody>
                </Table>
              </TableContainer>
            </Paper>
          </Grid>
        </Grid>
      )}
    </Container>
  );
}
