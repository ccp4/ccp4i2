'use client';

import { useState, useCallback } from 'react';
import { useRouter, useParams } from 'next/navigation';
import Link from 'next/link';
import {
  Container,
  Paper,
  Typography,
  Box,
  TextField,
  Button,
  Alert,
  CircularProgress,
  Chip,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogContentText,
  DialogActions,
  Divider,
} from '@mui/material';
import {
  LocalShipping,
  ArrowBack,
  Save,
  Delete,
  Person,
  Warning,
} from '@mui/icons-material';
import { useSWRConfig } from 'swr';
import { PageHeader } from '@/components/compounds/PageHeader';
import { useCompoundsApi, apiPatch, apiDelete } from '@/lib/compounds/api';
import { routes } from '@/lib/compounds/routes';
import { Supplier } from '@/types/compounds/models';

export default function SupplierDetailPage() {
  const router = useRouter();
  const params = useParams();
  const supplierId = params.id as string;
  const { mutate } = useSWRConfig();
  const api = useCompoundsApi();

  const { data: supplier, isLoading, error: loadError } = api.get<Supplier>(`suppliers/${supplierId}/`);

  const [name, setName] = useState<string | null>(null);
  const [initials, setInitials] = useState<string | null>(null);
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [success, setSuccess] = useState<string | null>(null);
  const [deleteDialogOpen, setDeleteDialogOpen] = useState(false);
  const [deleting, setDeleting] = useState(false);

  // Use local state if set, otherwise use fetched data
  const displayName = name !== null ? name : (supplier?.name || '');
  const displayInitials = initials !== null ? initials : (supplier?.initials || '');

  const handleSave = useCallback(async () => {
    if (!displayName.trim()) {
      setError('Supplier name is required');
      return;
    }

    setSaving(true);
    setError(null);
    setSuccess(null);

    try {
      await apiPatch(`suppliers/${supplierId}/`, {
        name: displayName.trim(),
        initials: displayInitials.trim() || null,
      });

      // Invalidate caches
      mutate(`/api/proxy/compounds/suppliers/${supplierId}/`);
      mutate('/api/proxy/compounds/suppliers/');

      setSuccess('Supplier updated successfully');
      // Reset local state to use fetched data
      setName(null);
      setInitials(null);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to update supplier');
    } finally {
      setSaving(false);
    }
  }, [supplierId, displayName, displayInitials, mutate]);

  const handleDelete = useCallback(async () => {
    setDeleting(true);
    setError(null);

    try {
      await apiDelete(`suppliers/${supplierId}/`);

      // Invalidate cache and navigate back
      mutate('/api/proxy/compounds/suppliers/');
      router.push(routes.registry.suppliers());
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to delete supplier');
      setDeleteDialogOpen(false);
    } finally {
      setDeleting(false);
    }
  }, [supplierId, mutate, router]);

  if (isLoading) {
    return (
      <Container maxWidth="md" sx={{ py: 4 }}>
        <Box sx={{ display: 'flex', justifyContent: 'center', py: 8 }}>
          <CircularProgress />
        </Box>
      </Container>
    );
  }

  if (loadError || !supplier) {
    return (
      <Container maxWidth="md" sx={{ py: 4 }}>
        <Alert severity="error">
          {loadError?.message || 'Supplier not found'}
        </Alert>
        <Button
          component={Link}
          href={routes.registry.suppliers()}
          startIcon={<ArrowBack />}
          sx={{ mt: 2 }}
        >
          Back to Suppliers
        </Button>
      </Container>
    );
  }

  const hasChanges = (name !== null && name !== supplier.name) ||
                     (initials !== null && initials !== (supplier.initials || ''));

  const totalAssociations = (supplier.compound_count || 0) + (supplier.batch_count || 0);

  return (
    <Container maxWidth="md" sx={{ py: 4 }}>
      <PageHeader
        breadcrumbs={[
          { label: 'Registry', href: routes.registry.targets() },
          { label: 'Suppliers', href: routes.registry.suppliers() },
          { label: supplier.name },
        ]}
      />

      <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 3 }}>
        <Button
          component={Link}
          href={routes.registry.suppliers()}
          startIcon={<ArrowBack />}
          size="small"
        >
          Back to Suppliers
        </Button>
      </Box>

      <Paper sx={{ p: 3 }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 3 }}>
          <LocalShipping color="primary" sx={{ fontSize: 32 }} />
          <Box sx={{ flex: 1 }}>
            <Typography variant="h5">
              {supplier.name}
              {supplier.is_current_user && (
                <Chip
                  icon={<Person sx={{ fontSize: 14 }} />}
                  label="Your supplier"
                  size="small"
                  color="primary"
                  sx={{ ml: 2 }}
                />
              )}
            </Typography>
            <Typography color="text.secondary" variant="body2">
              {supplier.compound_count || 0} compounds, {supplier.batch_count || 0} batches
            </Typography>
          </Box>
        </Box>

        <Divider sx={{ mb: 3 }} />

        {error && (
          <Alert severity="error" sx={{ mb: 2 }} onClose={() => setError(null)}>
            {error}
          </Alert>
        )}

        {success && (
          <Alert severity="success" sx={{ mb: 2 }} onClose={() => setSuccess(null)}>
            {success}
          </Alert>
        )}

        <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2.5 }}>
          <TextField
            label="Supplier Name"
            value={displayName}
            onChange={(e) => setName(e.target.value)}
            fullWidth
            required
          />

          <TextField
            label="Initials"
            value={displayInitials}
            onChange={(e) => setInitials(e.target.value)}
            fullWidth
            placeholder="e.g., NCL, AST"
            helperText="Short code used in barcodes"
          />
        </Box>

        <Box sx={{ mt: 4, display: 'flex', justifyContent: 'space-between' }}>
          <Button
            variant="outlined"
            color="error"
            startIcon={<Delete />}
            onClick={() => setDeleteDialogOpen(true)}
          >
            Delete Supplier
          </Button>

          <Button
            variant="contained"
            startIcon={saving ? <CircularProgress size={16} /> : <Save />}
            onClick={handleSave}
            disabled={saving || !hasChanges}
          >
            {saving ? 'Saving...' : 'Save Changes'}
          </Button>
        </Box>
      </Paper>

      {/* Delete Confirmation Dialog */}
      <Dialog open={deleteDialogOpen} onClose={() => setDeleteDialogOpen(false)}>
        <DialogTitle sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Warning color="error" />
          Delete Supplier
        </DialogTitle>
        <DialogContent>
          <DialogContentText>
            Are you sure you want to delete <strong>{supplier.name}</strong>?
          </DialogContentText>

          {totalAssociations > 0 && (
            <Alert severity="warning" sx={{ mt: 2 }}>
              This supplier is associated with:
              <ul style={{ margin: '8px 0 0 0', paddingLeft: 20 }}>
                {(supplier.compound_count || 0) > 0 && (
                  <li>{supplier.compound_count} compound(s)</li>
                )}
                {(supplier.batch_count || 0) > 0 && (
                  <li>{supplier.batch_count} batch(es)</li>
                )}
              </ul>
              <Typography variant="body2" sx={{ mt: 1 }}>
                These will have their supplier field cleared (set to empty), but the compounds and batches themselves will <strong>not</strong> be deleted.
              </Typography>
            </Alert>
          )}
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setDeleteDialogOpen(false)} disabled={deleting}>
            Cancel
          </Button>
          <Button
            onClick={handleDelete}
            color="error"
            variant="contained"
            disabled={deleting}
            startIcon={deleting ? <CircularProgress size={16} /> : <Delete />}
          >
            {deleting ? 'Deleting...' : 'Delete'}
          </Button>
        </DialogActions>
      </Dialog>
    </Container>
  );
}
