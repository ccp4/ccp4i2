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
  Stepper,
  Step,
  StepLabel,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Chip,
  FormControl,
  FormControlLabel,
  InputLabel,
  Select,
  MenuItem,
  Autocomplete,
  TextField,
  IconButton,
  Tooltip,
  LinearProgress,
  Checkbox,
} from '@mui/material';
import {
  Upload,
  ArrowBack,
  ArrowForward,
  CheckCircle,
  Error,
  Warning,
  Download,
  Refresh,
  Inventory,
} from '@mui/icons-material';
import { PageHeader } from '@/components/compounds/PageHeader';
import { SpreadsheetUpload, SpreadsheetData, FieldMapping, SpreadsheetRow } from '@/components/compounds/SpreadsheetUpload';
import { MoleculeChip } from '@/components/compounds/MoleculeView';
import { useCompoundsApi, apiPost } from '@/lib/compounds/api';
import { routes } from '@/lib/compounds/routes';

interface Target {
  id: string;
  name: string;
}

interface Supplier {
  id: string;
  name: string;
  initials: string | null;
}

interface ValidationResult {
  row: number;
  valid: boolean;
  errors: string[];
  warnings: string[];
  data: Record<string, any>;
}

interface ImportResult {
  success: boolean;
  formatted_id?: string;
  id?: string;
  error?: string;
  batch_created?: boolean;
  batch_number?: number;
  batch_error?: string;
}

// Mapping of stereo_comment display labels and common supplier terms to internal keys.
// Keys are lowercase; lookup is case-insensitive.
const STEREO_COMMENT_MAP: Record<string, string> = {
  // Internal keys (pass through)
  'unset': 'unset',
  'achiral': 'achiral',
  'racemic': 'racemic',
  'single_unknown': 'single_unknown',
  'r_enantiomer': 'r_enantiomer',
  's_enantiomer': 's_enantiomer',
  'non_racemic_mixture': 'non_racemic_mixture',
  'four_diastereomers': 'four_diastereomers',
  'two_diastereomers': 'two_diastereomers',
  'single_diastereomer_unknown': 'single_diastereomer_unknown',
  'rr_diastereomer': 'rr_diastereomer',
  'rs_diastereomer': 'rs_diastereomer',
  'sr_diastereomer': 'sr_diastereomer',
  'ss_diastereomer': 'ss_diastereomer',
  'epimer_mixture': 'epimer_mixture',
  'ez_mixture': 'ez_mixture',
  'e_isomer': 'e_isomer',
  'z_isomer': 'z_isomer',
  // Display labels (from STEREO_CHOICES)
  'racemic mixture': 'racemic',
  'single enantiomer, configuration unknown': 'single_unknown',
  'r enantiomer': 'r_enantiomer',
  's enantiomer': 's_enantiomer',
  'non-racemic stereoisomer mixture': 'non_racemic_mixture',
  'mixture of 4 diastereoisomers': 'four_diastereomers',
  'mixture of 2 diastereoisomers': 'two_diastereomers',
  'single diastereoisomer, configuration unknown': 'single_diastereomer_unknown',
  'rr diastereoisomer': 'rr_diastereomer',
  'rs diastereoisomer': 'rs_diastereomer',
  'sr diastereoisomer': 'sr_diastereomer',
  'ss diastereoisomer': 'ss_diastereomer',
  'mixture of epimers': 'epimer_mixture',
  'mixture of e and z isomers': 'ez_mixture',
  'e isomer': 'e_isomer',
  'z isomer': 'z_isomer',
  // Common supplier terms (Enamine, etc.)
  'racemic or presumed racemic or meso': 'racemic',
  'racemic diastereomer with known relative stereochemistry': 'racemic',
  'meso': 'achiral',
  'racemate': 'racemic',
  'rac': 'racemic',
  'single enantiomer': 'single_unknown',
  'enantiopure': 'single_unknown',
  'optically active': 'single_unknown',
  'e/z mixture': 'ez_mixture',
  'cis/trans mixture': 'ez_mixture',
  'n/a': 'achiral',
  'not applicable': 'achiral',
  'none': 'unset',
};

/**
 * Normalize a stereo_comment value from a spreadsheet to an internal key.
 * Returns the mapped key and an optional warning if the value was transformed.
 */
function normalizeStereoComment(value: string): { key: string; warning?: string } {
  if (!value || !value.trim()) {
    return { key: 'unset' };
  }
  const normalized = value.trim().toLowerCase();
  const mapped = STEREO_COMMENT_MAP[normalized];
  if (mapped) {
    // If the input was already the internal key, no warning
    if (normalized === mapped) {
      return { key: mapped };
    }
    return { key: mapped, warning: `Stereo "${value}" mapped to "${mapped}"` };
  }
  // Fallback: try underscore normalization (e.g. "R Enantiomer" → "r_enantiomer")
  const underscored = normalized.replace(/[\s-]+/g, '_');
  if (STEREO_COMMENT_MAP[underscored]) {
    return { key: underscored, warning: `Stereo "${value}" mapped to "${underscored}"` };
  }
  return { key: 'unset', warning: `Unknown stereo "${value}" — defaulting to "unset"` };
}

const COMPOUND_FIELD_MAPPINGS: FieldMapping[] = [
  { field: 'smiles', label: 'SMILES', required: true },
  { field: 'target', label: 'Target' },
  { field: 'stereo_comment', label: 'Stereochemistry' },
  { field: 'supplier', label: 'Supplier' },
  { field: 'supplier_ref', label: 'Supplier Reference' },
  { field: 'comments', label: 'Comments' },
  { field: 'labbook_number', label: 'Lab Notebook #' },
  { field: 'page_number', label: 'Page #' },
  // Batch fields (optional - used when auto-creating batches)
  { field: 'batch_amount', label: 'Batch Amount (mg)' },
  { field: 'batch_salt_code', label: 'Batch Salt Code' },
];

const STEPS = ['Upload File', 'Map Columns', 'Validate', 'Import'];

export default function ImportCompoundsPage() {
  const router = useRouter();
  const api = useCompoundsApi();

  const [activeStep, setActiveStep] = useState(0);
  const [spreadsheetData, setSpreadsheetData] = useState<SpreadsheetData | null>(null);
  const [columnMapping, setColumnMapping] = useState<Record<string, string>>({});
  const [defaultTarget, setDefaultTarget] = useState<string | null>(null);
  const [defaultSupplier, setDefaultSupplier] = useState<string | null>(null);
  const [autoCreateBatches, setAutoCreateBatches] = useState(true);
  const [validationResults, setValidationResults] = useState<ValidationResult[]>([]);
  const [importing, setImporting] = useState(false);
  const [importProgress, setImportProgress] = useState(0);
  const [importResults, setImportResults] = useState<ImportResult[]>([]);
  const [error, setError] = useState<string | null>(null);

  // Fetch targets and suppliers
  const { data: targetsData } = api.get<Target[]>('targets/');
  const { data: suppliersData } = api.get<Supplier[]>('suppliers/');

  const targets = targetsData || [];
  const suppliers = suppliersData || [];

  // Handle spreadsheet data loaded
  const handleDataLoaded = useCallback((data: SpreadsheetData) => {
    setSpreadsheetData(data);
    setActiveStep(1);
    setError(null);
  }, []);

  // Handle column mapping complete
  const handleMappingComplete = useCallback(
    (mapping: Record<string, string>, data: SpreadsheetData) => {
      setColumnMapping(mapping);
      setSpreadsheetData(data);

      // Validate the data
      const results: ValidationResult[] = data.rows.map((row, idx) => {
        const errors: string[] = [];
        const warnings: string[] = [];
        const extractedData: Record<string, any> = {};

        // Extract SMILES (required)
        const smilesCol = mapping.smiles;
        if (smilesCol && row[smilesCol]) {
          extractedData.smiles = String(row[smilesCol]).trim();
        } else {
          errors.push('Missing SMILES structure');
        }

        // Extract other fields
        for (const [field, column] of Object.entries(mapping)) {
          if (field === 'smiles') continue;
          if (column && row[column]) {
            const rawValue = String(row[column]).trim();
            if (field === 'stereo_comment') {
              const { key, warning } = normalizeStereoComment(rawValue);
              extractedData[field] = key;
              if (warning) warnings.push(warning);
            } else {
              extractedData[field] = rawValue;
            }
          }
        }

        // Add default target/supplier if set
        if (defaultTarget && !extractedData.target) {
          extractedData.target = defaultTarget;
        }
        if (defaultSupplier && !extractedData.supplier) {
          extractedData.supplier = defaultSupplier;
        }

        // Validate target
        if (!extractedData.target) {
          warnings.push('No target specified - will use default');
        }

        // Basic SMILES validation (very simple check)
        if (extractedData.smiles && extractedData.smiles.length < 2) {
          errors.push('SMILES appears invalid (too short)');
        }

        return {
          row: idx + 1,
          valid: errors.length === 0,
          errors,
          warnings,
          data: extractedData,
        };
      });

      setValidationResults(results);
      setActiveStep(2);
    },
    [defaultTarget, defaultSupplier]
  );

  // Re-validate when defaults change
  const handleRevalidate = useCallback(() => {
    if (!spreadsheetData) return;
    handleMappingComplete(columnMapping, spreadsheetData);
  }, [spreadsheetData, columnMapping, handleMappingComplete]);

  // Start import
  const handleStartImport = useCallback(async () => {
    const validRows = validationResults.filter((r) => r.valid);
    if (validRows.length === 0) {
      setError('No valid rows to import');
      return;
    }

    setImporting(true);
    setImportProgress(0);
    setImportResults([]);
    setError(null);
    setActiveStep(3);

    const results: ImportResult[] = [];

    for (let i = 0; i < validRows.length; i++) {
      const validation = validRows[i];
      try {
        // Prepare compound data
        const compoundData: Record<string, any> = {
          smiles: validation.data.smiles,
        };

        // Map target name to ID
        if (validation.data.target) {
          const target = targets.find(
            (t) =>
              t.id === validation.data.target ||
              t.name.toLowerCase() === validation.data.target.toLowerCase()
          );
          if (target) {
            compoundData.target = target.id;
          }
        }
        if (!compoundData.target && defaultTarget) {
          compoundData.target = defaultTarget;
        }

        // Map supplier name to ID
        if (validation.data.supplier) {
          const supplier = suppliers.find(
            (s) =>
              s.id === validation.data.supplier ||
              s.name.toLowerCase() === validation.data.supplier.toLowerCase()
          );
          if (supplier) {
            compoundData.supplier = supplier.id;
          }
        }
        if (!compoundData.supplier && defaultSupplier) {
          compoundData.supplier = defaultSupplier;
        }

        // Copy other fields
        if (validation.data.supplier_ref) {
          compoundData.supplier_ref = validation.data.supplier_ref;
        }
        if (validation.data.comments) {
          compoundData.comments = validation.data.comments;
        }
        if (validation.data.stereo_comment) {
          compoundData.stereo_comment = validation.data.stereo_comment;
        }
        if (validation.data.labbook_number) {
          compoundData.labbook_number = parseInt(validation.data.labbook_number);
        }
        if (validation.data.page_number) {
          compoundData.page_number = parseInt(validation.data.page_number);
        }

        // Submit compound to API
        const result = await apiPost<{ id: string; formatted_id: string }>(
          'compounds/',
          compoundData
        );

        // Create batch if auto-create is enabled
        let batchCreated = false;
        let batchNumber: number | undefined;
        let batchError: string | undefined;

        if (autoCreateBatches) {
          try {
            const batchData: Record<string, any> = {
              compound: result.id,
            };

            // Use supplier from compound if available
            if (compoundData.supplier) {
              batchData.supplier = compoundData.supplier;
            }

            // Add batch-specific fields from spreadsheet if mapped
            if (validation.data.batch_amount) {
              batchData.amount = validation.data.batch_amount;
            }
            if (validation.data.batch_salt_code) {
              batchData.salt_code = validation.data.batch_salt_code;
            }

            const batchResult = await apiPost<{ id: string; batch_number: number }>(
              'batches/',
              batchData
            );
            batchCreated = true;
            batchNumber = batchResult.batch_number;
          } catch (batchErr) {
            batchError = String(batchErr) || 'Failed to create batch';
          }
        }

        results.push({
          success: true,
          formatted_id: result.formatted_id,
          id: result.id,
          batch_created: batchCreated,
          batch_number: batchNumber,
          batch_error: batchError,
        });
      } catch (e) {
        results.push({
          success: false,
          error: String(e) || 'Failed to import',
        });
      }

      setImportProgress(((i + 1) / validRows.length) * 100);
      setImportResults([...results]);
    }

    setImporting(false);
  }, [validationResults, targets, suppliers, defaultTarget, defaultSupplier, autoCreateBatches]);

  // Summary stats
  const validCount = useMemo(
    () => validationResults.filter((r) => r.valid).length,
    [validationResults]
  );
  const errorCount = useMemo(
    () => validationResults.filter((r) => !r.valid).length,
    [validationResults]
  );
  const successCount = useMemo(
    () => importResults.filter((r) => r.success).length,
    [importResults]
  );
  const failedCount = useMemo(
    () => importResults.filter((r) => !r.success).length,
    [importResults]
  );
  const batchesCreatedCount = useMemo(
    () => importResults.filter((r) => r.batch_created).length,
    [importResults]
  );

  const handleReset = useCallback(() => {
    setActiveStep(0);
    setSpreadsheetData(null);
    setColumnMapping({});
    setValidationResults([]);
    setImportResults([]);
    setImportProgress(0);
    setError(null);
  }, []);

  return (
    <Container maxWidth="lg" sx={{ py: 4 }}>
      <PageHeader
        breadcrumbs={[
          { label: 'Registry', href: routes.registry.targets() },
          { label: 'Bulk Import' },
        ]}
      />

      <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 3 }}>
        <Button
          component={Link}
          href={routes.registry.new()}
          startIcon={<ArrowBack />}
          size="small"
        >
          Single Compound
        </Button>
        <Typography variant="h4" sx={{ flex: 1 }}>
          Bulk Import Compounds
        </Typography>
      </Box>

      {/* Stepper */}
      <Paper sx={{ p: 2, mb: 3 }}>
        <Stepper activeStep={activeStep}>
          {STEPS.map((label) => (
            <Step key={label}>
              <StepLabel>{label}</StepLabel>
            </Step>
          ))}
        </Stepper>
      </Paper>

      {/* Step 0: Upload File */}
      {activeStep === 0 && (
        <SpreadsheetUpload
          title="Upload Compound Data"
          fieldMappings={COMPOUND_FIELD_MAPPINGS}
          onDataLoaded={handleDataLoaded}
          onMappingComplete={handleMappingComplete}
          showColumnMapping={true}
          previewRows={5}
        />
      )}

      {/* Step 1: Column Mapping (handled by SpreadsheetUpload) */}
      {activeStep === 1 && spreadsheetData && (
        <SpreadsheetUpload
          title="Map Columns"
          fieldMappings={COMPOUND_FIELD_MAPPINGS}
          onMappingComplete={handleMappingComplete}
          showColumnMapping={true}
          previewRows={10}
        />
      )}

      {/* Step 2: Validation */}
      {activeStep === 2 && (
        <Paper sx={{ p: 3 }}>
          <Typography variant="h6" gutterBottom>
            Validation Results
          </Typography>

          {/* Default values */}
          <Box sx={{ display: 'flex', gap: 2, mb: 2, flexWrap: 'wrap', alignItems: 'center' }}>
            <Autocomplete
              options={targets}
              getOptionLabel={(option) => option.name}
              value={targets.find((t) => t.id === defaultTarget) || null}
              onChange={(_, newValue) => setDefaultTarget(newValue?.id || null)}
              sx={{ minWidth: 250 }}
              renderInput={(params) => (
                <TextField {...params} label="Default Target" size="small" />
              )}
            />
            <Autocomplete
              options={suppliers}
              getOptionLabel={(option) => option.name}
              value={suppliers.find((s) => s.id === defaultSupplier) || null}
              onChange={(_, newValue) => setDefaultSupplier(newValue?.id || null)}
              sx={{ minWidth: 250 }}
              renderInput={(params) => (
                <TextField {...params} label="Default Supplier" size="small" />
              )}
            />
            <Button
              variant="outlined"
              onClick={handleRevalidate}
              startIcon={<Refresh />}
            >
              Re-validate
            </Button>
          </Box>

          {/* Auto-create batches option */}
          <Alert
            severity="info"
            icon={<Inventory />}
            sx={{ mb: 3 }}
            action={
              <FormControlLabel
                control={
                  <Checkbox
                    checked={autoCreateBatches}
                    onChange={(e) => setAutoCreateBatches(e.target.checked)}
                  />
                }
                label=""
              />
            }
          >
            <Typography variant="body2" component="span">
              <strong>Auto-create batches</strong> for imported compounds.
              {columnMapping.batch_amount || columnMapping.batch_salt_code ? (
                <> Batch data from spreadsheet columns will be used.</>
              ) : (
                <> Batches will inherit supplier from compound.</>
              )}
            </Typography>
          </Alert>

          {/* Summary */}
          <Box sx={{ display: 'flex', gap: 2, mb: 3 }}>
            <Chip
              icon={<CheckCircle />}
              label={`${validCount} valid`}
              color="success"
              variant="outlined"
            />
            <Chip
              icon={<Error />}
              label={`${errorCount} with errors`}
              color="error"
              variant="outlined"
            />
          </Box>

          {/* Validation table */}
          <TableContainer sx={{ maxHeight: 400 }}>
            <Table size="small" stickyHeader>
              <TableHead>
                <TableRow>
                  <TableCell sx={{ fontWeight: 600 }}>Row</TableCell>
                  <TableCell sx={{ fontWeight: 600 }}>Status</TableCell>
                  <TableCell sx={{ fontWeight: 600 }}>SMILES</TableCell>
                  <TableCell sx={{ fontWeight: 600 }}>Target</TableCell>
                  <TableCell sx={{ fontWeight: 600 }}>Issues</TableCell>
                </TableRow>
              </TableHead>
              <TableBody>
                {validationResults.slice(0, 50).map((result) => (
                  <TableRow
                    key={result.row}
                    sx={{
                      bgcolor: result.valid ? undefined : 'error.50',
                    }}
                  >
                    <TableCell>{result.row}</TableCell>
                    <TableCell>
                      {result.valid ? (
                        <CheckCircle color="success" fontSize="small" />
                      ) : (
                        <Error color="error" fontSize="small" />
                      )}
                    </TableCell>
                    <TableCell sx={{ maxWidth: 200, overflow: 'hidden', textOverflow: 'ellipsis' }}>
                      <Typography variant="caption" fontFamily="monospace">
                        {result.data.smiles || '-'}
                      </Typography>
                    </TableCell>
                    <TableCell>
                      {result.data.target || (
                        <Typography variant="caption" color="text.secondary">
                          (default)
                        </Typography>
                      )}
                    </TableCell>
                    <TableCell>
                      {result.errors.map((err, i) => (
                        <Chip
                          key={i}
                          label={err}
                          size="small"
                          color="error"
                          sx={{ mr: 0.5 }}
                        />
                      ))}
                      {result.warnings.map((warn, i) => (
                        <Chip
                          key={i}
                          label={warn}
                          size="small"
                          color="warning"
                          sx={{ mr: 0.5 }}
                        />
                      ))}
                    </TableCell>
                  </TableRow>
                ))}
              </TableBody>
            </Table>
          </TableContainer>

          {validationResults.length > 50 && (
            <Typography variant="caption" color="text.secondary" sx={{ mt: 1 }}>
              Showing first 50 of {validationResults.length} rows
            </Typography>
          )}

          {error && (
            <Alert severity="error" sx={{ mt: 2 }}>
              {error}
            </Alert>
          )}

          {/* Actions */}
          <Box sx={{ display: 'flex', justifyContent: 'flex-end', gap: 2, mt: 3 }}>
            <Button variant="outlined" onClick={handleReset}>
              Start Over
            </Button>
            <Button
              variant="contained"
              disabled={validCount === 0}
              onClick={handleStartImport}
              endIcon={<ArrowForward />}
            >
              Import {validCount} Compounds
            </Button>
          </Box>
        </Paper>
      )}

      {/* Step 3: Import Progress */}
      {activeStep === 3 && (
        <Paper sx={{ p: 3 }}>
          <Typography variant="h6" gutterBottom>
            {importing ? 'Importing...' : 'Import Complete'}
          </Typography>

          {importing && (
            <Box sx={{ mb: 3 }}>
              <LinearProgress variant="determinate" value={importProgress} />
              <Typography variant="body2" color="text.secondary" sx={{ mt: 1 }}>
                {Math.round(importProgress)}% complete ({importResults.length} of {validCount})
              </Typography>
            </Box>
          )}

          {!importing && (
            <>
              {/* Summary */}
              <Box sx={{ display: 'flex', gap: 2, mb: 3, flexWrap: 'wrap' }}>
                <Chip
                  icon={<CheckCircle />}
                  label={`${successCount} compounds imported`}
                  color="success"
                  variant="outlined"
                />
                {batchesCreatedCount > 0 && (
                  <Chip
                    icon={<Inventory />}
                    label={`${batchesCreatedCount} batches created`}
                    color="info"
                    variant="outlined"
                  />
                )}
                {failedCount > 0 && (
                  <Chip
                    icon={<Error />}
                    label={`${failedCount} failed`}
                    color="error"
                    variant="outlined"
                  />
                )}
              </Box>

              {/* Results table */}
              <TableContainer sx={{ maxHeight: 400, mb: 3 }}>
                <Table size="small" stickyHeader>
                  <TableHead>
                    <TableRow>
                      <TableCell sx={{ fontWeight: 600 }}>#</TableCell>
                      <TableCell sx={{ fontWeight: 600 }}>Status</TableCell>
                      <TableCell sx={{ fontWeight: 600 }}>Compound ID</TableCell>
                      <TableCell sx={{ fontWeight: 600 }}>Batch</TableCell>
                      <TableCell sx={{ fontWeight: 600 }}>Error</TableCell>
                    </TableRow>
                  </TableHead>
                  <TableBody>
                    {importResults.map((result, idx) => (
                      <TableRow
                        key={idx}
                        sx={{ bgcolor: result.success ? undefined : 'error.50' }}
                      >
                        <TableCell>{idx + 1}</TableCell>
                        <TableCell>
                          {result.success ? (
                            <CheckCircle color="success" fontSize="small" />
                          ) : (
                            <Error color="error" fontSize="small" />
                          )}
                        </TableCell>
                        <TableCell>
                          {result.formatted_id ? (
                            <Link href={routes.registry.compound(result.id!)}>
                              <Chip
                                label={result.formatted_id}
                                size="small"
                                color="success"
                                clickable
                              />
                            </Link>
                          ) : (
                            '-'
                          )}
                        </TableCell>
                        <TableCell>
                          {result.batch_created ? (
                            <Chip
                              icon={<Inventory />}
                              label={`#${result.batch_number}`}
                              size="small"
                              color="info"
                              variant="outlined"
                            />
                          ) : result.batch_error ? (
                            <Tooltip title={result.batch_error}>
                              <Chip
                                icon={<Warning />}
                                label="Failed"
                                size="small"
                                color="warning"
                                variant="outlined"
                              />
                            </Tooltip>
                          ) : (
                            '-'
                          )}
                        </TableCell>
                        <TableCell>
                          <Typography variant="caption" color="error">
                            {result.error || ''}
                          </Typography>
                        </TableCell>
                      </TableRow>
                    ))}
                  </TableBody>
                </Table>
              </TableContainer>

              {/* Actions */}
              <Box sx={{ display: 'flex', justifyContent: 'flex-end', gap: 2 }}>
                <Button variant="outlined" onClick={handleReset} startIcon={<Upload />}>
                  Import More
                </Button>
                <Button
                  variant="contained"
                  component={Link}
                  href={routes.registry.search()}
                >
                  Go to Registry
                </Button>
              </Box>
            </>
          )}
        </Paper>
      )}
    </Container>
  );
}
