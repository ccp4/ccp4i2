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
  MenuItem,
} from '@mui/material';
import { Close, Medication, Save, Science } from '@mui/icons-material';
import { apiPatch, useCompoundsApi } from '@/lib/compounds/api';
import { MoleculeChip } from '@/components/compounds/MoleculeView';
import type { Compound, Supplier, Target } from '@/types/compounds/models';

const STEREO_CHOICES = [
  { value: 'unset', label: 'Unset' },
  { value: 'achiral', label: 'Achiral' },
  { value: 'racemic', label: 'Racemic mixture' },
  { value: 'single_unknown', label: 'Single enantiomer, configuration unknown' },
  { value: 'r_enantiomer', label: 'R enantiomer' },
  { value: 's_enantiomer', label: 'S enantiomer' },
  { value: 'non_racemic_mixture', label: 'Non-racemic stereoisomer mixture' },
  { value: 'four_diastereomers', label: 'Mixture of 4 diastereoisomers' },
  { value: 'two_diastereomers', label: 'Mixture of 2 diastereoisomers' },
  { value: 'single_diastereomer_unknown', label: 'Single diastereoisomer, configuration unknown' },
  { value: 'rr_diastereomer', label: 'RR diastereoisomer' },
  { value: 'rs_diastereomer', label: 'RS diastereoisomer' },
  { value: 'sr_diastereomer', label: 'SR diastereoisomer' },
  { value: 'ss_diastereomer', label: 'SS diastereoisomer' },
  { value: 'epimer_mixture', label: 'Mixture of epimers' },
  { value: 'ez_mixture', label: 'Mixture of E and Z isomers' },
  { value: 'e_isomer', label: 'E isomer' },
  { value: 'z_isomer', label: 'Z isomer' },
  ...Array.from({ length: 20 }, (_, i) => ({
    value: `isomer_${i + 1}`,
    label: `Isomer ${i + 1}`,
  })),
];

interface CompoundEditDialogProps {
  open: boolean;
  onClose: () => void;
  onSaved: (compound: Compound) => void;
  compound: Compound;
}

export function CompoundEditDialog({
  open,
  onClose,
  onSaved,
  compound,
}: CompoundEditDialogProps) {
  const api = useCompoundsApi();
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Form state
  const [smiles, setSmiles] = useState('');
  const [target, setTarget] = useState<Target | null>(null);
  const [supplier, setSupplier] = useState<Supplier | null>(null);
  const [supplierRef, setSupplierRef] = useState('');
  const [stereoComment, setStereoComment] = useState('unset');
  const [labbookNumber, setLabbookNumber] = useState('');
  const [pageNumber, setPageNumber] = useState('');
  const [compoundNumber, setCompoundNumber] = useState('');
  const [comments, setComments] = useState('');

  // Fetch suppliers and targets for autocomplete
  const { data: suppliers } = api.get<Supplier[]>('suppliers/');
  const { data: targets } = api.get<Target[]>('targets/');

  // Initialize form from compound when dialog opens
  useEffect(() => {
    if (open && compound) {
      setSmiles(compound.smiles || '');
      setSupplierRef(compound.supplier_ref || '');
      setStereoComment(compound.stereo_comment || 'unset');
      setLabbookNumber(compound.labbook_number?.toString() || '');
      setPageNumber(compound.page_number?.toString() || '');
      setCompoundNumber(compound.compound_number?.toString() || '');
      setComments(compound.comments || '');
      setError(null);
    }
  }, [open, compound]);

  // Set supplier once data is loaded
  useEffect(() => {
    if (suppliers && compound?.supplier) {
      const found = suppliers.find((s) => s.id === compound.supplier);
      setSupplier(found || null);
    }
  }, [suppliers, compound?.supplier, open]);

  // Set target once data is loaded
  useEffect(() => {
    if (targets && compound?.target) {
      const found = targets.find((t) => t.id === compound.target);
      setTarget(found || null);
    }
  }, [targets, compound?.target, open]);

  const handleSave = async () => {
    setSaving(true);
    setError(null);

    try {
      const payload: Record<string, any> = {
        target: target?.id || null,
        supplier: supplier?.id || null,
        supplier_ref: supplierRef.trim() || null,
        stereo_comment: stereoComment,
        labbook_number: labbookNumber ? parseInt(labbookNumber, 10) : null,
        page_number: pageNumber ? parseInt(pageNumber, 10) : null,
        compound_number: compoundNumber ? parseInt(compoundNumber, 10) : null,
        comments: comments.trim() || null,
      };

      // Only send SMILES if it changed (triggers recomputation of canonical SMILES etc.)
      const trimmedSmiles = smiles.trim();
      if (trimmedSmiles && trimmedSmiles !== compound.smiles) {
        payload.smiles = trimmedSmiles;
      }

      const updated = await apiPatch<Compound>(`compounds/${compound.id}/`, payload);

      onSaved(updated);
      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to update compound');
    } finally {
      setSaving(false);
    }
  };

  return (
    <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
      <DialogTitle sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Medication color="secondary" />
          Edit {compound.formatted_id}
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
          {/* SMILES editing */}
          <Box>
            <Box sx={{ display: 'flex', gap: 2, alignItems: 'flex-start' }}>
              <TextField
                label="SMILES"
                value={smiles}
                onChange={(e) => setSmiles(e.target.value)}
                fullWidth
                size="small"
                helperText={
                  smiles.trim() !== compound.smiles
                    ? 'Structure will be recomputed on save'
                    : 'Update SMILES to refine stereochemistry'
                }
                InputProps={{
                  startAdornment: <Science sx={{ mr: 1, color: 'text.secondary' }} />,
                  sx: { fontFamily: 'monospace', fontSize: '0.85rem' },
                }}
              />
              {smiles.trim() && (
                <Box sx={{ flexShrink: 0 }}>
                  <MoleculeChip smiles={smiles.trim()} size={80} />
                </Box>
              )}
            </Box>
          </Box>

          <Autocomplete
            options={targets || []}
            getOptionLabel={(option) => option.name}
            value={target}
            onChange={(_, newValue) => setTarget(newValue)}
            renderInput={(params) => (
              <TextField
                {...params}
                label="Target"
                size="small"
                placeholder="Select a target"
              />
            )}
            isOptionEqualToValue={(option, value) => option.id === value.id}
          />

          <TextField
            label="Stereo Comment"
            select
            value={stereoComment}
            onChange={(e) => setStereoComment(e.target.value)}
            fullWidth
            size="small"
          >
            {STEREO_CHOICES.map((option) => (
              <MenuItem key={option.value} value={option.value}>
                {option.label}
              </MenuItem>
            ))}
          </TextField>

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
            <TextField
              label="Compound #"
              value={compoundNumber}
              onChange={(e) => setCompoundNumber(e.target.value)}
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
            placeholder="Any additional notes about this compound"
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
