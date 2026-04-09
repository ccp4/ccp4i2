'use client';

import { useState, useEffect } from 'react';
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  TextField,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  Box,
  Alert,
  IconButton,
  CircularProgress,
  Typography,
} from '@mui/material';
import { Close, Edit, Science } from '@mui/icons-material';
import { useCompoundsApi } from '@/lib/compounds/api';
import type { DilutionSeries } from '@/types/compounds/models';

interface DilutionSeriesEditDialogProps {
  open: boolean;
  onClose: () => void;
  dilutionSeries: DilutionSeries;
  onSave: () => void;
}

const UNIT_OPTIONS: { value: 'nM' | 'uM' | 'mM'; label: string }[] = [
  { value: 'nM', label: 'nM (nanomolar)' },
  { value: 'uM', label: 'µM (micromolar)' },
  { value: 'mM', label: 'mM (millimolar)' },
];

function parseConcentrations(input: string): number[] | null {
  const parts = input.split(/[,\s]+/).filter((p) => p.trim() !== '');
  const numbers: number[] = [];

  for (const part of parts) {
    const num = parseFloat(part.trim());
    if (isNaN(num) || num < 0) {
      return null;
    }
    numbers.push(num);
  }

  return numbers.length > 0 ? numbers : null;
}

export function DilutionSeriesEditDialog({
  open,
  onClose,
  dilutionSeries,
  onSave,
}: DilutionSeriesEditDialogProps) {
  const api = useCompoundsApi();
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Form state
  const [concentrationsInput, setConcentrationsInput] = useState('');
  const [unit, setUnit] = useState<'nM' | 'uM' | 'mM'>('nM');

  // Reset form when dilution series changes
  useEffect(() => {
    setConcentrationsInput(dilutionSeries.concentrations.join(', '));
    setUnit(dilutionSeries.unit);
    setError(null);
  }, [dilutionSeries]);

  const handleSave = async () => {
    const concentrations = parseConcentrations(concentrationsInput);

    if (!concentrations) {
      setError('Please enter valid concentration values (comma-separated numbers)');
      return;
    }

    setSaving(true);
    setError(null);

    try {
      await api.patch(`dilution-series/${dilutionSeries.id}/`, {
        concentrations,
        unit,
      });
      onSave();
      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to save');
    } finally {
      setSaving(false);
    }
  };

  const parsedConcentrations = parseConcentrations(concentrationsInput);
  const isValid = parsedConcentrations !== null && parsedConcentrations.length > 0;

  return (
    <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
      <DialogTitle sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Edit color="primary" />
          Edit Dilution Series
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
            label="Concentrations"
            value={concentrationsInput}
            onChange={(e) => setConcentrationsInput(e.target.value)}
            fullWidth
            required
            autoFocus
            size="small"
            placeholder="e.g., 100, 30, 10, 3, 1, 0.3, 0.1"
            helperText="Enter comma-separated concentration values (highest to lowest)"
          />

          {isValid && (
            <Typography
              variant="body2"
              sx={{ fontFamily: 'monospace', color: 'success.main', pl: 1 }}
            >
              {parsedConcentrations.length} points: {parsedConcentrations.join(', ')}
            </Typography>
          )}

          <FormControl fullWidth size="small">
            <InputLabel>Unit</InputLabel>
            <Select
              value={unit}
              onChange={(e) => setUnit(e.target.value as 'nM' | 'uM' | 'mM')}
              label="Unit"
            >
              {UNIT_OPTIONS.map((option) => (
                <MenuItem key={option.value} value={option.value}>
                  {option.label}
                </MenuItem>
              ))}
            </Select>
          </FormControl>
        </Box>
      </DialogContent>

      <DialogActions>
        <Button onClick={onClose} disabled={saving}>
          Cancel
        </Button>
        <Button
          variant="contained"
          onClick={handleSave}
          disabled={saving || !isValid}
          startIcon={saving ? <CircularProgress size={16} /> : null}
        >
          {saving ? 'Saving...' : 'Save Changes'}
        </Button>
      </DialogActions>
    </Dialog>
  );
}
