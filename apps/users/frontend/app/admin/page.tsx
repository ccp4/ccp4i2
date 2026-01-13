'use client';

import { useEffect, useState } from 'react';
import {
  Container,
  Typography,
  Box,
  Paper,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Chip,
  Button,
  Alert,
  CircularProgress,
  IconButton,
  Tooltip,
  Divider,
  FormControlLabel,
  Checkbox,
  Stack,
} from '@mui/material';
import {
  AdminPanelSettings,
  PersonAdd,
  PersonRemove,
  Refresh,
  CloudUpload,
  CheckCircle,
  Error as ErrorIcon,
} from '@mui/icons-material';
import Link from 'next/link';
import { apiGet, apiPost } from '../../lib/users/api';

interface UserProfile {
  is_platform_admin: boolean;
  legacy_username: string;
  legacy_display_name: string;
  imported_at: string | null;
  first_login_at: string | null;
  last_seen_at: string | null;
}

interface User {
  id: number;
  username: string;
  email: string;
  display_name: string;
  is_admin: boolean;
  is_active: boolean;
  first_name?: string;
  last_name?: string;
  profile?: UserProfile;
}

interface CurrentUser {
  id: number;
  username: string;
  email: string;
  first_name: string;
  last_name: string;
  display_name: string;
  is_admin: boolean;
  profile: UserProfile;
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

export default function AdminPage() {
  const [currentUser, setCurrentUser] = useState<CurrentUser | null>(null);
  const [users, setUsers] = useState<User[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [actionLoading, setActionLoading] = useState<number | null>(null);

  // Fixture import state
  const [importStatus, setImportStatus] = useState<ImportStatus | null>(null);
  const [importStatusLoading, setImportStatusLoading] = useState(false);
  const [importing, setImporting] = useState(false);
  const [importResult, setImportResult] = useState<ImportResult | null>(null);
  const [importError, setImportError] = useState<string | null>(null);
  const [usersFile, setUsersFile] = useState<File | null>(null);
  const [registryFile, setRegistryFile] = useState<File | null>(null);
  const [assaysFile, setAssaysFile] = useState<File | null>(null);
  const [dryRun, setDryRun] = useState(true);

  const fetchCurrentUser = async (): Promise<CurrentUser | null> => {
    try {
      return await apiGet<CurrentUser>('me');
    } catch (err) {
      console.error('Error fetching current user:', err);
      return null;
    }
  };

  const fetchUsers = async (): Promise<User[]> => {
    try {
      return await apiGet<User[]>('users');
    } catch (err: any) {
      if (err.status === 403) {
        throw new Error('You do not have permission to view users. Admin access required.');
      }
      throw new Error('Failed to fetch users');
    }
  };

  const loadData = async () => {
    setLoading(true);
    setError(null);
    try {
      const user = await fetchCurrentUser();
      setCurrentUser(user);

      if (user?.is_admin) {
        const userList = await fetchUsers();
        setUsers(userList);
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : 'An error occurred');
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    loadData();
  }, []);

  const handleGrantAdmin = async (userId: number) => {
    setActionLoading(userId);
    try {
      await apiPost(`users/${userId}/grant_admin`);
      await loadData();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to grant admin');
    } finally {
      setActionLoading(null);
    }
  };

  const handleRevokeAdmin = async (userId: number) => {
    setActionLoading(userId);
    try {
      await apiPost(`users/${userId}/revoke_admin`);
      await loadData();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to revoke admin');
    } finally {
      setActionLoading(null);
    }
  };

  // Fixture import functions
  const fetchImportStatus = async () => {
    setImportStatusLoading(true);
    setImportError(null);
    try {
      const response = await fetch('/api/proxy/compounds/admin/import-status/');
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

      const response = await fetch('/api/proxy/compounds/admin/import-legacy/', {
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

  // Fetch import status on mount
  useEffect(() => {
    fetchImportStatus();
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

  return (
    <Container maxWidth="lg" sx={{ py: 4 }}>
      <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 4 }}>
        <AdminPanelSettings sx={{ fontSize: 40, color: 'grey.600' }} />
        <Box>
          <Typography variant="h4" component="h1">
            Platform Admin
          </Typography>
          <Typography color="text.secondary">
            User management and platform settings
          </Typography>
        </Box>
        <Box sx={{ ml: 'auto' }}>
          <Button component={Link} href="/" variant="outlined">
            Back to Apps
          </Button>
        </Box>
      </Box>

      {error && (
        <Alert severity="error" sx={{ mb: 3 }} onClose={() => setError(null)}>
          {error}
        </Alert>
      )}

      {/* Current User Info */}
      <Paper sx={{ p: 3, mb: 3 }}>
        <Typography variant="h6" gutterBottom>
          Current User
        </Typography>
        {currentUser ? (
          <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
            <Typography>
              <strong>{currentUser.display_name}</strong> ({currentUser.email || currentUser.username})
            </Typography>
            {currentUser.is_admin && (
              <Chip label="Platform Admin" color="primary" size="small" />
            )}
          </Box>
        ) : (
          <Typography color="text.secondary">Not authenticated</Typography>
        )}
      </Paper>

      {/* User Management */}
      {currentUser?.is_admin ? (
        <Paper sx={{ p: 3 }}>
          <Box sx={{ display: 'flex', alignItems: 'center', mb: 2 }}>
            <Typography variant="h6">User Management</Typography>
            <Tooltip title="Refresh">
              <IconButton onClick={loadData} sx={{ ml: 1 }}>
                <Refresh />
              </IconButton>
            </Tooltip>
          </Box>

          <TableContainer>
            <Table>
              <TableHead>
                <TableRow>
                  <TableCell>Name</TableCell>
                  <TableCell>Email</TableCell>
                  <TableCell>Status</TableCell>
                  <TableCell>Role</TableCell>
                  <TableCell>Last Seen</TableCell>
                  <TableCell align="right">Actions</TableCell>
                </TableRow>
              </TableHead>
              <TableBody>
                {users.map((user) => (
                  <TableRow key={user.id}>
                    <TableCell>{user.display_name || user.username}</TableCell>
                    <TableCell>{user.email}</TableCell>
                    <TableCell>
                      <Chip
                        label={user.is_active ? 'Active' : 'Inactive'}
                        color={user.is_active ? 'success' : 'default'}
                        size="small"
                      />
                    </TableCell>
                    <TableCell>
                      {user.is_admin ? (
                        <Chip label="Admin" color="primary" size="small" />
                      ) : (
                        <Chip label="User" variant="outlined" size="small" />
                      )}
                    </TableCell>
                    <TableCell>
                      {user.profile?.last_seen_at
                        ? new Date(user.profile.last_seen_at).toLocaleDateString()
                        : 'Never'}
                    </TableCell>
                    <TableCell align="right">
                      {actionLoading === user.id ? (
                        <CircularProgress size={24} />
                      ) : user.is_admin ? (
                        <Tooltip title="Revoke admin access">
                          <IconButton
                            onClick={() => handleRevokeAdmin(user.id)}
                            color="warning"
                            disabled={user.id === currentUser?.id}
                          >
                            <PersonRemove />
                          </IconButton>
                        </Tooltip>
                      ) : (
                        <Tooltip title="Grant admin access">
                          <IconButton
                            onClick={() => handleGrantAdmin(user.id)}
                            color="primary"
                          >
                            <PersonAdd />
                          </IconButton>
                        </Tooltip>
                      )}
                    </TableCell>
                  </TableRow>
                ))}
                {users.length === 0 && (
                  <TableRow>
                    <TableCell colSpan={6} align="center">
                      <Typography color="text.secondary">No users found</Typography>
                    </TableCell>
                  </TableRow>
                )}
              </TableBody>
            </Table>
          </TableContainer>
        </Paper>
      ) : (
        <Alert severity="warning">
          You need platform admin access to manage users. Contact an administrator if you need access.
        </Alert>
      )}

      {/* Legacy Fixture Import - Admin only */}
      {currentUser?.is_admin && (
        <Paper sx={{ p: 3, mt: 3 }}>
          <Box sx={{ display: 'flex', alignItems: 'center', mb: 2 }}>
            <Typography variant="h6">Import Legacy Fixtures</Typography>
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
      )}
    </Container>
  );
}
