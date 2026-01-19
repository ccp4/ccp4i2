'use client';

import React, { useEffect, useState } from 'react';
import {
  Container,
  Typography,
  Box,
  Paper,
  Button,
  Alert,
  CircularProgress,
  IconButton,
  Tooltip,
  Divider,
  FormControlLabel,
  Checkbox,
  Stack,
  TextField,
} from '@mui/material';
import {
  CloudUpload,
  CheckCircle,
  Error as ErrorIcon,
  Refresh,
  ArrowBack,
  DarkMode,
  LightMode,
} from '@mui/icons-material';
import Link from 'next/link';
import { authenticatedFetch, apiGet } from '../../../lib/users/api';
import { useTheme } from '../../../theme/theme-provider';

interface CurrentUser {
  id: number;
  username: string;
  email: string;
  first_name: string;
  last_name: string;
  display_name: string;
  is_admin: boolean;
}

interface ImportStatus {
  users: {
    total: number;
    with_legacy_username: number;
  };
  registry: {
    suppliers: number;
    targets: number;
    compounds: number;
    batches: number;
    batch_qc_files: number;
    compound_templates: number;
  };
  assays: {
    dilution_series: number;
    protocols: number;
    assays: number;
    data_series: number;
    analysis_results: number;
    hypotheses: number;
  };
}

interface CCP4i2ImportStatus {
  projects: {
    total: number;
    groups: number;
    memberships: number;
    tags: number;
  };
  jobs: {
    total: number;
    value_keys: number;
    float_values: number;
    char_values: number;
    xdata: number;
  };
  files: {
    total: number;
    types: number;
    uses: number;
    imports: number;
    exports: number;
  };
}

interface CCP4i2ImportResult {
  dry_run: boolean;
  loaded?: boolean;
  output?: string;
  stats?: Record<string, number>;
  errors: string[];
}

interface ImportResult {
  dry_run: boolean;
  loaded?: boolean;
  users?: {
    total_records: number;
    created: number;
    updated: number;
    skipped: number;
  };
  registry?: {
    total_records: number;
    by_model: Record<string, number>;
  };
  assays?: {
    total_records: number;
    by_model: Record<string, number>;
  };
  errors: string[];
}

export default function ImportPage() {
  const { mode, toggleTheme } = useTheme();
  const [currentUser, setCurrentUser] = useState<CurrentUser | null>(null);
  const [loading, setLoading] = useState(true);

  // Fixture import state (compounds)
  const [importStatus, setImportStatus] = useState<ImportStatus | null>(null);
  const [importStatusLoading, setImportStatusLoading] = useState(false);
  const [importing, setImporting] = useState(false);
  const [importResult, setImportResult] = useState<ImportResult | null>(null);
  const [importError, setImportError] = useState<string | null>(null);
  const [usersFile, setUsersFile] = useState<File | null>(null);
  const [registryFile, setRegistryFile] = useState<File | null>(null);
  const [assaysFile, setAssaysFile] = useState<File | null>(null);
  const [dryRun, setDryRun] = useState(true);

  // CCP4i2 fixture import state
  const [ccp4i2Status, setCCP4i2Status] = useState<CCP4i2ImportStatus | null>(null);
  const [ccp4i2StatusLoading, setCCP4i2StatusLoading] = useState(false);
  const [ccp4i2Importing, setCCP4i2Importing] = useState(false);
  const [ccp4i2Result, setCCP4i2Result] = useState<CCP4i2ImportResult | null>(null);
  const [ccp4i2Error, setCCP4i2Error] = useState<string | null>(null);
  const [ccp4i2File, setCCP4i2File] = useState<File | null>(null);
  const [ccp4i2DryRun, setCCP4i2DryRun] = useState(true);
  const [ccp4i2RemapFrom, setCCP4i2RemapFrom] = useState('');
  const [ccp4i2RemapTo, setCCP4i2RemapTo] = useState('');

  const fetchCurrentUser = async (): Promise<CurrentUser | null> => {
    try {
      return await apiGet<CurrentUser>('me');
    } catch (err) {
      console.error('Error fetching current user:', err);
      return null;
    }
  };

  // Fixture import functions
  const fetchImportStatus = async () => {
    setImportStatusLoading(true);
    setImportError(null);
    try {
      const response = await authenticatedFetch('/api/proxy/compounds/admin/import-status/');
      if (!response.ok) {
        throw new Error(`Failed to fetch status: ${response.statusText}`);
      }
      const data = await response.json();
      setImportStatus(data);
    } catch (err) {
      setImportError(err instanceof Error ? err.message : 'Failed to fetch status');
    } finally {
      setImportStatusLoading(false);
    }
  };

  const handleImport = async () => {
    if (!usersFile && !registryFile && !assaysFile) {
      setImportError('Please select at least one fixture file');
      return;
    }

    setImporting(true);
    setImportError(null);
    setImportResult(null);

    try {
      const formData = new FormData();
      if (usersFile) {
        formData.append('users_fixture', usersFile);
      }
      if (registryFile) {
        formData.append('registry_fixture', registryFile);
      }
      if (assaysFile) {
        formData.append('assays_fixture', assaysFile);
      }
      formData.append('dry_run', dryRun.toString());

      const response = await authenticatedFetch('/api/proxy/compounds/admin/import-legacy/', {
        method: 'POST',
        body: formData,
      });

      const data = await response.json();

      if (!response.ok) {
        setImportResult(data);
        if (data.errors?.length > 0) {
          setImportError(data.errors.join('; '));
        }
      } else {
        setImportResult(data);
        if (!dryRun) {
          await fetchImportStatus();
        }
      }
    } catch (err) {
      setImportError(err instanceof Error ? err.message : 'Import failed');
    } finally {
      setImporting(false);
    }
  };

  // CCP4i2 fixture import functions
  const fetchCCP4i2Status = async () => {
    setCCP4i2StatusLoading(true);
    setCCP4i2Error(null);
    try {
      const response = await authenticatedFetch('/api/proxy/ccp4i2/admin/import-status/');
      if (!response.ok) {
        throw new Error(`Failed to fetch status: ${response.statusText}`);
      }
      const data = await response.json();
      setCCP4i2Status(data);
    } catch (err) {
      setCCP4i2Error(err instanceof Error ? err.message : 'Failed to fetch status');
    } finally {
      setCCP4i2StatusLoading(false);
    }
  };

  const handleCCP4i2Import = async () => {
    if (!ccp4i2File) {
      setCCP4i2Error('Please select a CCP4i2 fixture file');
      return;
    }

    setCCP4i2Importing(true);
    setCCP4i2Error(null);
    setCCP4i2Result(null);

    try {
      const formData = new FormData();
      formData.append('ccp4i2_fixture', ccp4i2File);
      formData.append('dry_run', ccp4i2DryRun.toString());
      if (ccp4i2RemapFrom && ccp4i2RemapTo) {
        formData.append('remap_from', ccp4i2RemapFrom);
        formData.append('remap_to', ccp4i2RemapTo);
      }

      const response = await authenticatedFetch('/api/proxy/ccp4i2/admin/import-legacy/', {
        method: 'POST',
        body: formData,
      });

      const data = await response.json();

      if (!response.ok) {
        setCCP4i2Result(data);
        if (data.errors?.length > 0) {
          setCCP4i2Error(data.errors.join('; '));
        }
      } else {
        setCCP4i2Result(data);
        if (!ccp4i2DryRun) {
          await fetchCCP4i2Status();
        }
      }
    } catch (err) {
      setCCP4i2Error(err instanceof Error ? err.message : 'Import failed');
    } finally {
      setCCP4i2Importing(false);
    }
  };

  useEffect(() => {
    const loadData = async () => {
      setLoading(true);
      const user = await fetchCurrentUser();
      setCurrentUser(user);
      setLoading(false);
    };
    loadData();
    fetchImportStatus();
    fetchCCP4i2Status();
  }, []);

  if (loading) {
    return (
      <Container maxWidth="lg" sx={{ py: 4 }}>
        <Box sx={{ display: 'flex', justifyContent: 'center', py: 8 }}>
          <CircularProgress />
        </Box>
      </Container>
    );
  }

  if (!currentUser?.is_admin) {
    return (
      <Container maxWidth="lg" sx={{ py: 4 }}>
        <Alert severity="warning">
          You need platform admin access to import data. Contact an administrator if you need access.
        </Alert>
        <Box sx={{ mt: 2 }}>
          <Button component={Link} href="/admin" startIcon={<ArrowBack />}>
            Back to Admin
          </Button>
        </Box>
      </Container>
    );
  }

  return (
    <Container maxWidth="lg" sx={{ py: 4 }}>
      <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 4 }}>
        <CloudUpload sx={{ fontSize: 40, color: 'grey.600' }} />
        <Box>
          <Typography variant="h4" component="h1">
            Data Import
          </Typography>
          <Typography color="text.secondary">
            Import legacy data from fixtures
          </Typography>
        </Box>
        <Box sx={{ ml: 'auto', display: 'flex', alignItems: 'center', gap: 1 }}>
          <Tooltip title={mode === 'light' ? 'Switch to dark mode' : 'Switch to light mode'}>
            <IconButton onClick={toggleTheme} sx={{ color: 'text.secondary' }}>
              {mode === 'light' ? <DarkMode /> : <LightMode />}
            </IconButton>
          </Tooltip>
          <Button component={Link} href="/admin" variant="outlined" startIcon={<ArrowBack />}>
            Back to Admin
          </Button>
        </Box>
      </Box>

      {/* Compounds Legacy Fixture Import */}
      <Paper sx={{ p: 3, mb: 3 }}>
        <Box sx={{ display: 'flex', alignItems: 'center', mb: 2 }}>
          <Typography variant="h6">Import Compounds Legacy Fixtures</Typography>
          <Tooltip title="Refresh database counts">
            <IconButton onClick={fetchImportStatus} sx={{ ml: 1 }} disabled={importStatusLoading}>
              <Refresh />
            </IconButton>
          </Tooltip>
        </Box>

        <Typography variant="body2" color="text.secondary" paragraph>
          Upload JSON fixtures from legacy apps. Import users first to preserve FK relationships,
          then registry and assays fixtures.
        </Typography>

        <Divider sx={{ my: 2 }} />

        {/* Current Database Status */}
        {importStatus && (
          <Box sx={{ mb: 3 }}>
            <Typography variant="subtitle2" gutterBottom>
              Current Database Counts
            </Typography>
            <Box sx={{ display: 'flex', gap: 4, flexWrap: 'wrap' }}>
              <Box>
                <Typography variant="caption" color="text.secondary">Users</Typography>
                <Typography variant="body2">
                  {importStatus.users.total} total ({importStatus.users.with_legacy_username} imported)
                </Typography>
              </Box>
              <Box>
                <Typography variant="caption" color="text.secondary">Registry</Typography>
                <Typography variant="body2">
                  {importStatus.registry.compounds} compounds, {importStatus.registry.batches} batches, {importStatus.registry.targets} targets
                </Typography>
              </Box>
              <Box>
                <Typography variant="caption" color="text.secondary">Assays</Typography>
                <Typography variant="body2">
                  {importStatus.assays.assays} assays, {importStatus.assays.data_series} data series, {importStatus.assays.protocols} protocols
                </Typography>
              </Box>
            </Box>
          </Box>
        )}

        <Stack spacing={2}>
          {/* Users Fixture Upload */}
          <Box>
            <Typography variant="subtitle2" gutterBottom>
              Users Fixture (auth.json) - Import first for FK references
            </Typography>
            <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
              <Button
                variant="outlined"
                component="label"
                startIcon={<CloudUpload />}
                size="small"
              >
                Select File
                <input
                  type="file"
                  hidden
                  accept=".json"
                  onChange={(e) => setUsersFile(e.target.files?.[0] || null)}
                />
              </Button>
              {usersFile && (
                <Typography variant="body2" color="text.secondary">
                  {usersFile.name} ({(usersFile.size / 1024).toFixed(1)} KB)
                </Typography>
              )}
            </Box>
          </Box>

          {/* Registry Fixture Upload */}
          <Box>
            <Typography variant="subtitle2" gutterBottom>
              Registry Fixture (RegisterCompounds.json)
            </Typography>
            <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
              <Button
                variant="outlined"
                component="label"
                startIcon={<CloudUpload />}
                size="small"
              >
                Select File
                <input
                  type="file"
                  hidden
                  accept=".json"
                  onChange={(e) => setRegistryFile(e.target.files?.[0] || null)}
                />
              </Button>
              {registryFile && (
                <Typography variant="body2" color="text.secondary">
                  {registryFile.name} ({(registryFile.size / 1024).toFixed(1)} KB)
                </Typography>
              )}
            </Box>
          </Box>

          {/* Assays Fixture Upload */}
          <Box>
            <Typography variant="subtitle2" gutterBottom>
              Assays Fixture (AssayCompounds.json)
            </Typography>
            <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
              <Button
                variant="outlined"
                component="label"
                startIcon={<CloudUpload />}
                size="small"
              >
                Select File
                <input
                  type="file"
                  hidden
                  accept=".json"
                  onChange={(e) => setAssaysFile(e.target.files?.[0] || null)}
                />
              </Button>
              {assaysFile && (
                <Typography variant="body2" color="text.secondary">
                  {assaysFile.name} ({(assaysFile.size / 1024).toFixed(1)} KB)
                </Typography>
              )}
            </Box>
          </Box>

          {/* Options */}
          <FormControlLabel
            control={
              <Checkbox
                checked={dryRun}
                onChange={(e) => setDryRun(e.target.checked)}
                size="small"
              />
            }
            label="Dry run (validate without loading)"
          />

          {/* Import Button */}
          <Box>
            <Button
              variant="contained"
              onClick={handleImport}
              disabled={importing || (!usersFile && !registryFile && !assaysFile)}
              startIcon={importing ? <CircularProgress size={16} /> : <CloudUpload />}
              size="small"
            >
              {importing ? 'Processing...' : dryRun ? 'Validate' : 'Import'}
            </Button>
          </Box>

          {/* Error Display */}
          {importError && (
            <Alert severity="error" icon={<ErrorIcon />} onClose={() => setImportError(null)}>
              {importError}
            </Alert>
          )}

          {/* Import Result */}
          {importResult && !importError && (
            <Alert
              severity={importResult.errors?.length ? 'warning' : 'success'}
              icon={importResult.errors?.length ? <ErrorIcon /> : <CheckCircle />}
            >
              <Typography variant="subtitle2" gutterBottom>
                {importResult.dry_run ? 'Validation Result' : 'Import Result'}
              </Typography>
              {importResult.users && (
                <Typography variant="body2">
                  Users: {importResult.users.total_records} records ({importResult.users.created} new, {importResult.users.updated} updated)
                </Typography>
              )}
              {importResult.registry && (
                <Typography variant="body2">
                  Registry: {importResult.registry.total_records} records
                </Typography>
              )}
              {importResult.assays && (
                <Typography variant="body2">
                  Assays: {importResult.assays.total_records} records
                </Typography>
              )}
              {!importResult.dry_run && importResult.loaded && (
                <Typography variant="body2" sx={{ fontWeight: 'bold', mt: 1 }}>
                  Successfully loaded into database
                </Typography>
              )}
            </Alert>
          )}
        </Stack>
      </Paper>

      {/* CCP4i2 Legacy Fixture Import */}
      <Paper sx={{ p: 3 }}>
        <Box sx={{ display: 'flex', alignItems: 'center', mb: 2 }}>
          <Typography variant="h6">Import CCP4i2 Legacy Fixtures</Typography>
          <Tooltip title="Refresh database counts">
            <IconButton onClick={fetchCCP4i2Status} sx={{ ml: 1 }} disabled={ccp4i2StatusLoading}>
              <Refresh />
            </IconButton>
          </Tooltip>
        </Box>

        <Typography variant="body2" color="text.secondary" paragraph>
          Upload a CCP4i2 dumpdata JSON fixture file. This imports legacy projects, jobs, files,
          and related data with automatic UUID and timestamp conversion.
        </Typography>

        <Divider sx={{ my: 2 }} />

        {/* Current Database Status */}
        {ccp4i2Status && (
          <Box sx={{ mb: 3 }}>
            <Typography variant="subtitle2" gutterBottom>
              Current CCP4i2 Database Counts
            </Typography>
            <Box sx={{ display: 'flex', gap: 4, flexWrap: 'wrap' }}>
              <Box>
                <Typography variant="caption" color="text.secondary">Projects</Typography>
                <Typography variant="body2">
                  {ccp4i2Status.projects.total} projects, {ccp4i2Status.projects.groups} groups, {ccp4i2Status.projects.tags} tags
                </Typography>
              </Box>
              <Box>
                <Typography variant="caption" color="text.secondary">Jobs</Typography>
                <Typography variant="body2">
                  {ccp4i2Status.jobs.total} jobs, {ccp4i2Status.jobs.xdata} xdata
                </Typography>
              </Box>
              <Box>
                <Typography variant="caption" color="text.secondary">Files</Typography>
                <Typography variant="body2">
                  {ccp4i2Status.files.total} files, {ccp4i2Status.files.types} types, {ccp4i2Status.files.uses} uses
                </Typography>
              </Box>
            </Box>
          </Box>
        )}

        <Stack spacing={2}>
          {/* CCP4i2 Fixture Upload */}
          <Box>
            <Typography variant="subtitle2" gutterBottom>
              CCP4i2 Fixture (dumpdata JSON)
            </Typography>
            <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
              <Button
                variant="outlined"
                component="label"
                startIcon={<CloudUpload />}
                size="small"
              >
                Select File
                <input
                  type="file"
                  hidden
                  accept=".json"
                  onChange={(e) => setCCP4i2File(e.target.files?.[0] || null)}
                />
              </Button>
              {ccp4i2File && (
                <Typography variant="body2" color="text.secondary">
                  {ccp4i2File.name} ({(ccp4i2File.size / 1024).toFixed(1)} KB)
                </Typography>
              )}
            </Box>
          </Box>

          {/* Directory Remapping (optional) */}
          <Box>
            <Typography variant="subtitle2" gutterBottom>
              Directory Remapping (optional)
            </Typography>
            <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
              Remap file paths from legacy directory structure to new location.
            </Typography>
            <Box sx={{ display: 'flex', gap: 2, flexWrap: 'wrap' }}>
              <TextField
                size="small"
                label="Remap from"
                placeholder="/legacy/ccp4i2/projects"
                value={ccp4i2RemapFrom}
                onChange={(e) => setCCP4i2RemapFrom(e.target.value)}
                sx={{ minWidth: 250 }}
              />
              <TextField
                size="small"
                label="Remap to"
                placeholder="/new/ccp4i2/projects"
                value={ccp4i2RemapTo}
                onChange={(e) => setCCP4i2RemapTo(e.target.value)}
                sx={{ minWidth: 250 }}
              />
            </Box>
          </Box>

          {/* Options */}
          <FormControlLabel
            control={
              <Checkbox
                checked={ccp4i2DryRun}
                onChange={(e) => setCCP4i2DryRun(e.target.checked)}
                size="small"
              />
            }
            label="Dry run (validate without loading)"
          />

          {/* Import Button */}
          <Box>
            <Button
              variant="contained"
              onClick={handleCCP4i2Import}
              disabled={ccp4i2Importing || !ccp4i2File}
              startIcon={ccp4i2Importing ? <CircularProgress size={16} /> : <CloudUpload />}
              size="small"
            >
              {ccp4i2Importing ? 'Processing...' : ccp4i2DryRun ? 'Validate' : 'Import'}
            </Button>
          </Box>

          {/* Error Display */}
          {ccp4i2Error && (
            <Alert severity="error" icon={<ErrorIcon />} onClose={() => setCCP4i2Error(null)}>
              {ccp4i2Error}
            </Alert>
          )}

          {/* Import Result */}
          {ccp4i2Result && !ccp4i2Error && (
            <Alert
              severity={ccp4i2Result.errors?.length ? 'warning' : 'success'}
              icon={ccp4i2Result.errors?.length ? <ErrorIcon /> : <CheckCircle />}
            >
              <Typography variant="subtitle2" gutterBottom>
                {ccp4i2Result.dry_run ? 'Validation Result' : 'Import Result'}
              </Typography>
              {ccp4i2Result.stats && Object.keys(ccp4i2Result.stats).length > 0 && (
                <Typography variant="body2">
                  {Object.entries(ccp4i2Result.stats).map(([key, value]) =>
                    `${key}: ${value}`
                  ).join(', ')}
                </Typography>
              )}
              {!ccp4i2Result.dry_run && ccp4i2Result.loaded && (
                <Typography variant="body2" sx={{ fontWeight: 'bold', mt: 1 }}>
                  Successfully loaded into database
                </Typography>
              )}
              {ccp4i2Result.output && (
                <Box sx={{ mt: 1, maxHeight: 200, overflow: 'auto' }}>
                  <Typography variant="caption" component="pre" sx={{ whiteSpace: 'pre-wrap', fontFamily: 'monospace' }}>
                    {ccp4i2Result.output}
                  </Typography>
                </Box>
              )}
            </Alert>
          )}
        </Stack>
      </Paper>
    </Container>
  );
}
