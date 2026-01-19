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
  Autocomplete,
} from '@mui/material';
import { Close, Inventory, Add } from '@mui/icons-material';
import { apiPost, useCompoundsApi } from '@/lib/compounds/api';
import type { Batch, Supplier } from '@/types/compounds/models';

interface BatchCreateDialogProps {
  open: boolean;
  onClose: () => void;
  onCreated: (batch: Batch) => void;
  compoundId: string;
  compoundFormattedId: string;
}

export function BatchCreateDialog({
  open,
  onClose,
  onCreated,
  compoundId,
  compoundFormattedId,
}: BatchCreateDialogProps) {
  const api = useCompoundsApi();
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Form state
  const [supplier, setSupplier] = useState<Supplier | null>(null);
  const [supplierRef, setSupplierRef] = useState('');
  const [amount, setAmount] = useState('');
  const [saltCode, setSaltCode] = useState('');
  const [labbookNumber, setLabbookNumber] = useState('');
  const [pageNumber, setPageNumber] = useState('');
  const [comments, setComments] = useState('');

  // Fetch suppliers for autocomplete
  const { data: suppliers } = api.get<Supplier[]>('suppliers/');

  // Reset form when dialog opens
  useEffect(() => {
    if (open) {
      setSupplier(null);
      setSupplierRef('');
      setAmount('');
      setSaltCode('');
      setLabbookNumber('');
      setPageNumber('');
      setComments('');
      setError(null);
    }
  }, [open]);

  const handleSave = async () => {
    setSaving(true);
    setError(null);

    try {
      const newBatch = await apiPost<Batch>('batches/', {
        compound: compoundId,
        supplier: supplier?.id || null,
        supplier_ref: supplierRef.trim() || null,
        amount: amount.trim() || null,
        salt_code: saltCode.trim() || null,
        labbook_number: labbookNumber ? parseInt(labbookNumber, 10) : null,
        page_number: pageNumber ? parseInt(pageNumber, 10) : null,
        comments: comments.trim() || null,
      });

      onCreated(newBatch);
      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to create batch');
    } finally {
      setSaving(false);
    }
  };

  return (
    <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
      <DialogTitle sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Inventory color="info" />
          New Batch for {compoundFormattedId}
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
          <Autocomplete
            options={suppliers || []}
            getOptionLabel={(option) => option.name}
            value={supplier}
            onChange={(_, newValue) => setSupplier(newValue)}
            renderInput={(params) => (
              <TextField
                {...params}
                label="Supplier"
                size="small"
                placeholder="Select a supplier"
              />
            )}
            isOptionEqualToValue={(option, value) => option.id === value.id}
          />

          <TextField
            label="Supplier Reference"
            value={supplierRef}
            onChange={(e) => setSupplierRef(e.target.value)}
            fullWidth
            size="small"
            placeholder="e.g., PO-12345, CAT-ABC123"
            helperText="Supplier's product code or purchase order"
          />

          <Box sx={{ display: 'flex', gap: 2 }}>
            <TextField
              label="Amount (mg)"
              value={amount}
              onChange={(e) => setAmount(e.target.value)}
              size="small"
              type="number"
              inputProps={{ step: '0.01', min: 0 }}
              sx={{ flex: 1 }}
            />
            <TextField
              label="Salt Code"
              value={saltCode}
              onChange={(e) => setSaltCode(e.target.value)}
              size="small"
              placeholder="e.g., HCl, TFA"
              sx={{ flex: 1 }}
            />
          </Box>

          <Box sx={{ display: 'flex', gap: 2 }}>
            <TextField
              label="Lab Book #"
              value={labbookNumber}
              onChange={(e) => setLabbookNumber(e.target.value)}
              size="small"
              type="number"
              inputProps={{ min: 1 }}
              sx={{ flex: 1 }}
            />
            <TextField
              label="Page #"
              value={pageNumber}
              onChange={(e) => setPageNumber(e.target.value)}
              size="small"
              type="number"
              inputProps={{ min: 1 }}
              sx={{ flex: 1 }}
            />
          </Box>

          <TextField
            label="Comments"
            value={comments}
            onChange={(e) => setComments(e.target.value)}
            fullWidth
            size="small"
            multiline
            rows={2}
            placeholder="Any additional notes about this batch"
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
          disabled={saving}
          startIcon={saving ? <CircularProgress size={16} /> : <Add />}
        >
          {saving ? 'Creating...' : 'Create Batch'}
        </Button>
      </DialogActions>
    </Dialog>
  );
}
