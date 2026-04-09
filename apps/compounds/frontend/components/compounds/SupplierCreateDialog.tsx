'use client';

import { useState, useEffect } from 'react';
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  TextField,
  Box,
  Alert,
  IconButton,
  CircularProgress,
} from '@mui/material';
import { Close, LocalShipping, Add } from '@mui/icons-material';
import { apiPost } from '@/lib/compounds/api';
import type { Supplier } from '@/types/compounds/models';

interface SupplierCreateDialogProps {
  open: boolean;
  onClose: () => void;
  onCreated: (supplier: Supplier) => void;
}

export function SupplierCreateDialog({
  open,
  onClose,
  onCreated,
}: SupplierCreateDialogProps) {
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Form state
  const [name, setName] = useState('');
  const [initials, setInitials] = useState('');

  // Reset form when dialog opens
  useEffect(() => {
    if (open) {
      setName('');
      setInitials('');
      setError(null);
    }
  }, [open]);

  const handleSave = async () => {
    if (!name.trim()) {
      setError('Supplier name is required');
      return;
    }

    setSaving(true);
    setError(null);

    try {
      const newSupplier = await apiPost<Supplier>('suppliers/', {
        name: name.trim(),
        initials: initials.trim() || null,
      });

      onCreated(newSupplier);
      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to create supplier');
    } finally {
      setSaving(false);
    }
  };

  return (
    <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
      <DialogTitle sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <LocalShipping color="primary" />
          New Supplier
        </Box>
        <IconButton onClick={onClose} size="small">
          <Close />
        </IconButton>
      </DialogTitle>

      <DialogContent dividers>
        {error && (
          <Alert severity="error" sx={{ mb: 2 }}>
            {error}
          </Alert>
        )}

        <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2.5, pt: 1 }}>
          <TextField
            label="Supplier Name"
            value={name}
            onChange={(e) => setName(e.target.value)}
            fullWidth
            required
            autoFocus
            size="small"
            placeholder="e.g., Sigma-Aldrich, Internal Synthesis"
          />

          <TextField
            label="Initials (optional)"
            value={initials}
            onChange={(e) => setInitials(e.target.value)}
            fullWidth
            size="small"
            placeholder="e.g., NCL, AST"
            helperText="Short code used in barcodes"
          />
        </Box>
      </DialogContent>

      <DialogActions>
        <Button onClick={onClose} disabled={saving}>
          Cancel
        </Button>
        <Button
          variant="contained"
          onClick={handleSave}
          disabled={saving || !name.trim()}
          startIcon={saving ? <CircularProgress size={16} /> : <Add />}
        >
          {saving ? 'Creating...' : 'Create Supplier'}
        </Button>
      </DialogActions>
    </Dialog>
  );
}
