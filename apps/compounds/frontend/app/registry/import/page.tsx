'use client';

import { useState, useCallback, useMemo, useEffect } from 'react';
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
  Description,
  Science,
} from '@mui/icons-material';
import { PageHeader } from '@/components/compounds/PageHeader';
import { SpreadsheetUpload, SpreadsheetData, FieldMapping, SpreadsheetRow, SpreadsheetHyperlinks } from '@/components/compounds/SpreadsheetUpload';
import { MoleculeChip } from '@/components/compounds/MoleculeView';
import { useCompoundsApi, apiPost, apiGet } from '@/lib/compounds/api';
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
  // Isomer series (for multi-stereocentre compounds)
  ...Object.fromEntries(
    Array.from({ length: 20 }, (_, i) => [`isomer_${i + 1}`, `isomer_${i + 1}`])
  ),
  ...Object.fromEntries(
    Array.from({ length: 20 }, (_, i) => [`isomer ${i + 1}`, `isomer_${i + 1}`])
  ),
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
  { field: 'supplier_ref', label: 'Supplier Reference / Peptide ID' },
  { field: 'comments', label: 'Comments' },
  { field: 'labbook_number', label: 'Lab Notebook #' },
  { field: 'page_number', label: 'Page #' },
  // ELN-linked fields (Phase 1)
  { field: 'eln_reference', label: 'ELN Reference (e.g., KF001)' },
  { field: 'helm_notation', label: 'HELM Notation' },
  { field: 'sequence_display', label: 'Sequence' },
  // Batch fields (optional - used when auto-creating batches)
  { field: 'batch_amount', label: 'Batch Amount (mg)' },
  { field: 'batch_salt_code', label: 'Batch Salt Code' },
];

/**
 * Categorize column as compound document or batch QC based on heuristics.
 * Returns 'compound' for structure files, 'batch' for QC/analysis files.
 */
function categorizeDocumentColumn(columnName: string): 'compound' | 'batch' {
  const lower = columnName.toLowerCase();
  // Compound documents: structure files (Chemdraw, etc.)
  if (lower.includes('chemdraw') || lower.includes('cdx') || lower.includes('structure') ||
      lower.includes('mol file') || lower.includes('.mol')) {
    return 'compound';
  }
  // Batch QC documents: analysis/QC files
  if (lower.includes('lcms') || lower.includes('lc-ms') || lower.includes('hrms') ||
      lower.includes('hr-ms') || lower.includes('mass') || lower.includes('purity') ||
      lower.includes('nmr') || lower.includes('spectrum') || lower.includes('spectra') ||
      lower.includes('chromatogram') || lower.includes('hplc') || lower.includes('qc')) {
    return 'batch';
  }
  // Default to compound (safer - user can reclassify)
  return 'compound';
}

/**
 * Infer compound document kind from column header name.
 * Returns 'chemdraw' | 'other' for compound documents.
 */
function inferCompoundDocumentKind(columnName: string): 'chemdraw' | 'other' {
  const lower = columnName.toLowerCase();
  if (lower.includes('chemdraw') || lower.includes('cdx') || lower.includes('structure')) {
    return 'chemdraw';
  }
  return 'other';
}

/**
 * Infer batch QC file kind from column header name.
 * Returns 'lcms' | 'hrms' | 'nmr' | 'hplc' | 'other'.
 */
function inferBatchQCKind(columnName: string): 'lcms' | 'hrms' | 'nmr' | 'hplc' | 'other' {
  const lower = columnName.toLowerCase();
  if (lower.includes('lcms') || lower.includes('lc-ms')) return 'lcms';
  if (lower.includes('hrms') || lower.includes('hr-ms')) return 'hrms';
  if (lower.includes('nmr')) return 'nmr';
  if (lower.includes('hplc')) return 'hplc';
  return 'other';
}

const STEPS = ['Upload File', 'Map Columns', 'Validate', 'Import'];

export default function ImportCompoundsPage() {
  const router = useRouter();
  const searchParams = useSearchParams();
  const api = useCompoundsApi();

  const [activeStep, setActiveStep] = useState(0);
  const [spreadsheetData, setSpreadsheetData] = useState<SpreadsheetData | null>(null);
  const [columnMapping, setColumnMapping] = useState<Record<string, string>>({});
  const [defaultTarget, setDefaultTarget] = useState<string | null>(null);
  const [defaultSupplier, setDefaultSupplier] = useState<string | null>(null);
  const [autoCreateBatches, setAutoCreateBatches] = useState(true);
  const [validationResults, setValidationResults] = useState<ValidationResult[]>([]);
  const [checkingDuplicates, setCheckingDuplicates] = useState(false);
  const [importing, setImporting] = useState(false);
  const [importProgress, setImportProgress] = useState(0);
  const [importResults, setImportResults] = useState<ImportResult[]>([]);
  const [error, setError] = useState<string | null>(null);

  // Document categorization: which hyperlinked columns go to compound vs batch
  const [compoundDocColumns, setCompoundDocColumns] = useState<string[]>([]);
  const [batchQCColumns, setBatchQCColumns] = useState<string[]>([]);
  // All hyperlinked columns detected in the spreadsheet
  const [allHyperlinkedColumns, setAllHyperlinkedColumns] = useState<string[]>([]);

  // Fetch targets and suppliers
  const { data: targetsData } = api.get<Target[]>('targets/');
  const { data: suppliersData } = api.get<Supplier[]>('suppliers/');

  const targets = targetsData || [];
  const suppliers = suppliersData || [];

  // If arrived with ?target=<id>, use it as the Default Target so rows that
  // don't name a target inherit it (e.g. coming from a target-filtered list).
  useEffect(() => {
    const targetParam = searchParams?.get('target');
    if (targetParam && targets.length && !defaultTarget) {
      const match = targets.find((t) => t.id === targetParam);
      if (match) setDefaultTarget(match.id);
    }
  }, [searchParams, targets, defaultTarget]);

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
        const rowHyperlinks = data.hyperlinks?.[idx] || {};

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
            // Extract hyperlink URL if present (for ELN reference)
            const hyperlink = rowHyperlinks[column];
            if (hyperlink) {
              // Store the hyperlink URL with a _url suffix
              extractedData[`${field}_url`] = hyperlink;
            }
          }
        }

        // Collect all hyperlinked columns for document creation
        // These are columns with URLs that aren't mapped to specific fields
        const hyperlinkedColumns: Array<{ column: string; label: string; url: string }> = [];
        for (const [column, url] of Object.entries(rowHyperlinks)) {
          if (url && String(url).startsWith('http')) {
            const displayText = row[column] ? String(row[column]).trim() : column;
            hyperlinkedColumns.push({ column, label: displayText, url });
          }
        }
        if (hyperlinkedColumns.length > 0) {
          extractedData._hyperlinkedColumns = hyperlinkedColumns;
        }

        // Add default target/supplier if set
        if (defaultTarget && !extractedData.target) {
          extractedData.target = defaultTarget;
        }
        if (defaultSupplier && !extractedData.supplier) {
          extractedData.supplier = defaultSupplier;
        }

        // Validate that target (from row or default) resolves to a known Target.
        // The import step silently drops unmapped values, so the server 400s
        // with "target: This field is required"; surface it here instead.
        if (!extractedData.target) {
          errors.push('No target specified — pick a Default Target above, or add a Target column');
        } else {
          const targetValue = String(extractedData.target);
          const resolved = targets.find(
            (t) => t.id === targetValue || t.name.toLowerCase() === targetValue.toLowerCase()
          );
          if (!resolved) {
            errors.push(
              `Target "${targetValue}" does not exist — create it first or pick a Default Target above`
            );
          }
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

      // Detect all unique hyperlinked columns across the spreadsheet
      // Exclude columns that are already mapped to known fields
      const knownMappedFields = new Set([
        'smiles', 'target', 'supplier', 'stereo_comment', 'comments',
        'labbook_number', 'page_number', 'helm_notation', 'sequence_display',
        'batch_amount', 'batch_salt_code', 'supplier_ref', 'eln_reference'
      ]);
      const mappedColumns = new Set(
        Object.entries(mapping)
          .filter(([field]) => knownMappedFields.has(field))
          .map(([, col]) => col)
      );

      const hyperlinkedColSet = new Set<string>();
      data.hyperlinks?.forEach(rowHyperlinks => {
        Object.keys(rowHyperlinks).forEach(col => {
          if (!mappedColumns.has(col)) {
            hyperlinkedColSet.add(col);
          }
        });
      });
      const detectedColumns = Array.from(hyperlinkedColSet);
      setAllHyperlinkedColumns(detectedColumns);

      // Auto-categorize using heuristics (only if not already set)
      if (detectedColumns.length > 0) {
        const newCompoundCols: string[] = [];
        const newBatchCols: string[] = [];
        detectedColumns.forEach(col => {
          const category = categorizeDocumentColumn(col);
          if (category === 'compound') {
            newCompoundCols.push(col);
          } else {
            newBatchCols.push(col);
          }
        });
        setCompoundDocColumns(newCompoundCols);
        setBatchQCColumns(newBatchCols);

        // Auto-enable batch creation if any batch QC columns detected
        if (newBatchCols.length > 0) {
          setAutoCreateBatches(true);
        }
      }

      setActiveStep(2);

      // Async duplicate check against existing compounds
      checkForDuplicates(results);
    },
    [defaultTarget, defaultSupplier, targets]
  );

  // Check SMILES against existing compounds in the database
  const checkForDuplicates = useCallback(async (results: ValidationResult[]) => {
    const rowsWithSmiles = results.filter((r) => r.data.smiles && r.valid);
    if (rowsWithSmiles.length === 0) return;

    setCheckingDuplicates(true);
    try {
      // Deduplicate by (smiles, stereo_comment) pairs to minimize API calls
      const pairSet = new Set<string>();
      const pairs: { smiles: string; stereo: string }[] = [];
      for (const r of rowsWithSmiles) {
        const smiles = r.data.smiles as string;
        const stereo = (r.data.stereo_comment as string) || '';
        const key = `${smiles}|||${stereo}`;
        if (!pairSet.has(key)) {
          pairSet.add(key);
          pairs.push({ smiles, stereo });
        }
      }

      // Maps "smiles|||stereo" -> existing formatted_id
      const duplicateMap: Record<string, string> = {};

      // Check each unique pair (in parallel batches of 10)
      for (let i = 0; i < pairs.length; i += 10) {
        const batch = pairs.slice(i, i + 10);
        const checks = batch.map(async ({ smiles, stereo }) => {
          try {
            let url = `compounds/resolve_by_smiles/?smiles=${encodeURIComponent(smiles)}`;
            if (stereo) {
              url += `&stereo_comment=${encodeURIComponent(stereo)}`;
            }
            const result = await apiGet<{
              found: boolean;
              compound: { formatted_id: string } | null;
            }>(url);
            if (result.found && result.compound) {
              duplicateMap[`${smiles}|||${stereo}`] = result.compound.formatted_id;
            }
          } catch {
            // Ignore - backend validation will catch it during import
          }
        });
        await Promise.all(checks);
      }

      // Update validation results with duplicate errors
      if (Object.keys(duplicateMap).length > 0) {
        setValidationResults((prev) =>
          prev.map((r) => {
            if (!r.data.smiles) return r;
            const key = `${r.data.smiles}|||${(r.data.stereo_comment as string) || ''}`;
            const existingId = duplicateMap[key];
            if (existingId) {
              return {
                ...r,
                valid: false,
                errors: [...r.errors, `Already registered as ${existingId}`],
              };
            }
            return r;
          })
        );
      }
    } finally {
      setCheckingDuplicates(false);
    }
  }, []);

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

        // ELN-linked fields
        if (validation.data.helm_notation) {
          compoundData.helm_notation = validation.data.helm_notation;
        }
        if (validation.data.sequence_display) {
          compoundData.sequence_display = validation.data.sequence_display;
        }

        // Handle ELN reference - lookup or create LabNotebookEntry
        // The display text (e.g., "KF001") is in validation.data.eln_reference
        // The hyperlink URL is in validation.data.eln_reference_url
        if (validation.data.eln_reference) {
          const elnRef = String(validation.data.eln_reference).trim();
          const elnUrl = validation.data.eln_reference_url || '';
          // Parse label like 'KF001' -> initials='KF', sequence=1
          const match = elnRef.match(/^([A-Za-z]{1,4})(\d+)$/);
          if (match && compoundData.supplier) {
            const initials = match[1].toUpperCase();
            const seqNum = parseInt(match[2], 10);
            try {
              // Check if notebook entry exists
              const existingEntries = await apiGet<{ results: { id: string }[] }>(
                `notebook-entries/?label=${encodeURIComponent(elnRef)}`
              );
              if (existingEntries.results && existingEntries.results.length > 0) {
                compoundData.notebook_entry = existingEntries.results[0].id;
              } else {
                // Create new notebook entry with the ELN URL from the hyperlink
                const newEntry = await apiPost<{ id: string }>(
                  'notebook-entries/',
                  {
                    supplier: compoundData.supplier,
                    sequence_number: seqNum,
                    url: elnUrl,
                  }
                );
                compoundData.notebook_entry = newEntry.id;
              }
            } catch {
              // Ignore - notebook entry creation is optional
            }
          }
        }

        // Submit compound to API
        const result = await apiPost<{ id: string; formatted_id: string }>(
          'compounds/',
          compoundData
        );

        // Create batch if auto-create is enabled (or if batch QC columns need it)
        let batchCreated = false;
        let batchNumber: number | undefined;
        let batchId: string | undefined;
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
            batchId = batchResult.id;
          } catch (batchErr) {
            batchError = String(batchErr) || 'Failed to create batch';
          }
        }

        // Create documents for hyperlinked columns based on user categorization
        if (validation.data._hyperlinkedColumns) {
          const hyperlinkedCols = validation.data._hyperlinkedColumns as Array<{
            column: string;
            label: string;
            url: string;
          }>;
          for (const doc of hyperlinkedCols) {
            // Skip ELN reference column - it's handled separately
            if (doc.column === columnMapping.eln_reference) continue;
            // Skip SMILES and other non-document columns
            const mappedField = Object.entries(columnMapping).find(([, col]) => col === doc.column);
            if (mappedField && ['smiles', 'target', 'supplier', 'stereo_comment', 'comments',
                               'labbook_number', 'page_number', 'helm_notation', 'sequence_display',
                               'batch_amount', 'batch_salt_code', 'supplier_ref', 'eln_reference'].includes(mappedField[0])) {
              continue;
            }

            // Check if this column is categorized as compound document
            if (compoundDocColumns.includes(doc.column)) {
              const kind = inferCompoundDocumentKind(doc.column);
              try {
                await apiPost('compound-documents/', {
                  compound: result.id,
                  kind,
                  label: doc.label,
                  url: doc.url,
                });
              } catch {
                // Ignore document creation errors - not critical
              }
            }
            // Check if this column is categorized as batch QC document
            else if (batchQCColumns.includes(doc.column) && batchId) {
              const kind = inferBatchQCKind(doc.column);
              try {
                await apiPost('batch-qc-files/', {
                  batch: batchId,
                  kind,
                  label: doc.label,
                  url: doc.url,
                });
              } catch {
                // Ignore QC file creation errors - not critical
              }
            }
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
  }, [validationResults, targets, suppliers, defaultTarget, defaultSupplier, autoCreateBatches, compoundDocColumns, batchQCColumns, columnMapping]);

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

      {/* Steps 0 + 1: single SpreadsheetUpload instance handles both upload
          and column mapping. Mounting a different instance per step loses the
          internal `data` state, making the user upload twice. */}
      {(activeStep === 0 || activeStep === 1) && (
        <SpreadsheetUpload
          title={activeStep === 0 ? 'Upload Compound Data' : 'Map Columns'}
          fieldMappings={COMPOUND_FIELD_MAPPINGS}
          onDataLoaded={handleDataLoaded}
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
            severity={batchQCColumns.length > 0 ? 'warning' : 'info'}
            icon={<Inventory />}
            sx={{ mb: 2 }}
            action={
              <FormControlLabel
                control={
                  <Checkbox
                    checked={autoCreateBatches}
                    onChange={(e) => setAutoCreateBatches(e.target.checked)}
                    disabled={batchQCColumns.length > 0}
                  />
                }
                label=""
              />
            }
          >
            <Typography variant="body2" component="span">
              <strong>Auto-create batches</strong> for imported compounds.
              {batchQCColumns.length > 0 ? (
                <> <strong>Required</strong> — QC documents need a batch to attach to.</>
              ) : columnMapping.batch_amount || columnMapping.batch_salt_code ? (
                <> Batch data from spreadsheet columns will be used.</>
              ) : (
                <> Batches will inherit supplier from compound.</>
              )}
            </Typography>
          </Alert>

          {/* Document Categorization - only show if hyperlinked columns detected */}
          {allHyperlinkedColumns.length > 0 && (
            <Paper variant="outlined" sx={{ p: 2, mb: 3 }}>
              <Typography variant="subtitle1" gutterBottom sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                <Description color="primary" />
                Document Columns
              </Typography>
              <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
                Hyperlinked columns are auto-categorized. Adjust as needed:
                <strong> Structure files</strong> attach to compounds, <strong>QC files</strong> attach to batches.
              </Typography>

              <Box sx={{ display: 'flex', gap: 2, flexWrap: 'wrap' }}>
                {/* Compound Documents */}
                <Autocomplete
                  multiple
                  options={allHyperlinkedColumns}
                  value={compoundDocColumns}
                  onChange={(_, newValue) => {
                    setCompoundDocColumns(newValue);
                    // Remove from batch QC if added here
                    setBatchQCColumns(prev => prev.filter(c => !newValue.includes(c)));
                  }}
                  renderTags={(value, getTagProps) =>
                    value.map((option, index) => (
                      <Chip
                        {...getTagProps({ index })}
                        key={option}
                        label={option}
                        size="small"
                        icon={<Science fontSize="small" />}
                        color="primary"
                        variant="outlined"
                      />
                    ))
                  }
                  renderInput={(params) => (
                    <TextField
                      {...params}
                      label="Compound Documents (structure files)"
                      size="small"
                      placeholder="Add column..."
                      helperText="Chemdraw, structure files"
                    />
                  )}
                  sx={{ minWidth: 300, flex: 1 }}
                />

                {/* Batch QC Documents */}
                <Autocomplete
                  multiple
                  options={allHyperlinkedColumns}
                  value={batchQCColumns}
                  onChange={(_, newValue) => {
                    setBatchQCColumns(newValue);
                    // Remove from compound docs if added here
                    setCompoundDocColumns(prev => prev.filter(c => !newValue.includes(c)));
                    // Auto-enable batch creation if QC docs selected
                    if (newValue.length > 0) {
                      setAutoCreateBatches(true);
                    }
                  }}
                  renderTags={(value, getTagProps) =>
                    value.map((option, index) => (
                      <Chip
                        {...getTagProps({ index })}
                        key={option}
                        label={option}
                        size="small"
                        icon={<Inventory fontSize="small" />}
                        color="secondary"
                        variant="outlined"
                      />
                    ))
                  }
                  renderInput={(params) => (
                    <TextField
                      {...params}
                      label="Batch QC Documents (analysis files)"
                      size="small"
                      placeholder="Add column..."
                      helperText="LCMS, HRMS, NMR, purity reports"
                    />
                  )}
                  sx={{ minWidth: 300, flex: 1 }}
                />
              </Box>
            </Paper>
          )}

          {/* Summary */}
          <Box sx={{ display: 'flex', gap: 2, mb: 3, alignItems: 'center' }}>
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
            {checkingDuplicates && (
              <Chip
                icon={<CircularProgress size={16} />}
                label="Checking for duplicates..."
                variant="outlined"
              />
            )}
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
              disabled={validCount === 0 || checkingDuplicates}
              onClick={handleStartImport}
              endIcon={<ArrowForward />}
            >
              {checkingDuplicates ? 'Checking duplicates...' : `Import ${validCount} Compounds`}
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
