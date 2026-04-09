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
import { Close, Inventory, Save } from '@mui/icons-material';
import { apiPatch, useCompoundsApi } from '@/lib/compounds/api';
import type { Batch, Supplier } from '@/types/compounds/models';

interface BatchEditDialogProps {
  open: boolean;
  onClose: () => void;
  onSaved: (batch: Batch) => void;
  batch: Batch;
}

export function BatchEditDialog({
  open,
  onClose,
  onSaved,
  batch,
}: BatchEditDialogProps) {
  const api = useCompoundsApi();
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Form state
  const [supplier, setSupplier] = useState<Supplier | null>(null);
  const [supplierRef, setSupplierRef] = useState('');
  const [amount, setAmount] = useState('');
  const [saltCode, setSaltCode] = useState('');
  const [molecularWeight, setMolecularWeight] = useState('');
  const [labbookNumber, setLabbookNumber] = useState('');
  const [pageNumber, setPageNumber] = useState('');
  const [comments, setComments] = useState('');

  // Fetch suppliers for autocomplete
  const { data: suppliers } = api.get<Supplier[]>('suppliers/');

  // Initialize form from batch when dialog opens
  useEffect(() => {
    if (open && batch) {
      setSupplierRef(batch.supplier_ref || '');
      setAmount(batch.amount != null ? batch.amount : '');
      setSaltCode(batch.salt_code || '');
      setMolecularWeight(batch.molecular_weight != null ? String(batch.molecular_weight) : '');
      setLabbookNumber(batch.labbook_number != null ? String(batch.labbook_number) : '');
      setPageNumber(batch.page_number != null ? String(batch.page_number) : '');
      setComments(batch.comments || '');
      setError(null);
    }
  }, [open, batch]);

  // Set supplier once data is loaded
  useEffect(() => {
    if (suppliers && batch?.supplier) {
      const found = suppliers.find((s) => s.id === batch.supplier);
      setSupplier(found || null);
    } else if (suppliers && !batch?.supplier) {
      setSupplier(null);
    }
  }, [suppliers, batch?.supplier, open]);

  const handleSave = async () => {
    setSaving(true);
    setError(null);

    try {
      const payload: Record<string, any> = {
        supplier: supplier?.id || null,
        supplier_ref: supplierRef.trim() || null,
        amount: amount.toString().trim() || null,
        salt_code: saltCode.trim() || null,
        molecular_weight: molecularWeight.trim() ? parseFloat(molecularWeight) : null,
        labbook_number: labbookNumber ? parseInt(labbookNumber, 10) : null,
        page_number: pageNumber ? parseInt(pageNumber, 10) : null,
        comments: comments.trim() || null,
      };

      const updated = await apiPatch<Batch>(`batches/${batch.id}/`, payload);

      onSaved(updated);
      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to update batch');
    } finally {
      setSaving(false);
    }
  };

  return (
    <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
      <DialogTitle sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Inventory color="info" />
          Edit Batch #{batch.batch_number}
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

          <TextField
            label="Molecular Weight (salt)"
            value={molecularWeight}
            onChange={(e) => setMolecularWeight(e.target.value)}
            fullWidth
            size="small"
            type="number"
            inputProps={{ step: '0.01', min: 0 }}
            helperText="Molecular weight of the salt form"
          />

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
            rows={3}
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
          startIcon={saving ? <CircularProgress size={16} /> : <Save />}
        >
          {saving ? 'Saving...' : 'Save Changes'}
        </Button>
      </DialogActions>
    </Dialog>
  );
}
