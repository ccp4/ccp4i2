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
  Typography,
  Alert,
  CircularProgress,
} from '@mui/material';
import { useCompoundsApi } from '@/lib/compounds/api';

export interface ProtocolThresholdInput {
  id: string;
  name: string;
  target_value?: number | null;
  poor_value?: number | null;
  threshold_scale?: 'log' | 'linear';
  kpi_unit?: string | null;
}

export interface ProtocolThresholdUpdate {
  target_value: number | null;
  poor_value: number | null;
  threshold_scale: 'log' | 'linear';
}

interface Props {
  open: boolean;
  protocol: ProtocolThresholdInput | null;
  onClose: () => void;
  /** Called with the saved values so callers can patch their local state. */
  onSaved: (id: string, update: ProtocolThresholdUpdate) => void;
}

/**
 * Minimal threshold-only edit dialog designed to be opened inline from
 * the aggregation views — touches only target_value, poor_value, and
 * threshold_scale. For full protocol editing use ProtocolEditDialog.
 */
export function ProtocolThresholdQuickEdit({
  open,
  protocol,
  onClose,
  onSaved,
}: Props) {
  const api = useCompoundsApi();
  const [targetStr, setTargetStr] = useState('');
  const [poorStr, setPoorStr] = useState('');
  const [scale, setScale] = useState<'log' | 'linear'>('log');
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (open && protocol) {
      setTargetStr(
        protocol.target_value != null ? String(protocol.target_value) : '',
      );
      setPoorStr(
        protocol.poor_value != null ? String(protocol.poor_value) : '',
      );
      setScale(protocol.threshold_scale || 'log');
      setError(null);
    }
  }, [open, protocol]);

  const handleSave = async () => {
    if (!protocol) return;

    const parse = (s: string): number | null => {
      if (s.trim() === '') return null;
      const n = Number(s);
      return Number.isFinite(n) ? n : null;
    };
    const targetValue = parse(targetStr);
    const poorValue = parse(poorStr);

    if (targetValue !== null && poorValue !== null && targetValue === poorValue) {
      setError(
        'Excellent and poor values must differ — direction of "better" is implied by their ordering.',
      );
      return;
    }

    setSaving(true);
    setError(null);
    try {
      await api.patch(`protocols/${protocol.id}/`, {
        target_value: targetValue,
        poor_value: poorValue,
        threshold_scale: scale,
      });
      onSaved(protocol.id, {
        target_value: targetValue,
        poor_value: poorValue,
        threshold_scale: scale,
      });
      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to save');
    } finally {
      setSaving(false);
    }
  };

  return (
    <Dialog open={open} onClose={onClose} maxWidth="xs" fullWidth>
      <DialogTitle sx={{ pb: 1 }}>
        <Typography variant="subtitle1" component="span" fontWeight={600}>
          Thresholds
        </Typography>
        {protocol && (
          <Typography
            variant="body2"
            color="text.secondary"
            sx={{ display: 'block', mt: 0.25 }}
          >
            {protocol.name}
          </Typography>
        )}
      </DialogTitle>
      <DialogContent dividers>
        {error && (
          <Alert severity="error" sx={{ mb: 2 }}>
            {error}
          </Alert>
        )}
        <Typography
          variant="caption"
          color="text.secondary"
          sx={{ display: 'block', mb: 1.5 }}
        >
          Absolute thresholds for colour-coding KPI values. Set excellent &lt;
          poor for potency metrics (IC50, EC50), or excellent &gt; poor for
          capacity metrics (solubility, permeability). Leave both blank to
          disable colouring.
          {protocol?.kpi_unit && (
            <>
              {' '}KPI unit: <strong>{protocol.kpi_unit}</strong>.
            </>
          )}
        </Typography>
        <Box sx={{ display: 'flex', gap: 1.5 }}>
          <TextField
            label="Excellent value"
            value={targetStr}
            onChange={(e) => setTargetStr(e.target.value)}
            type="number"
            size="small"
            fullWidth
            autoFocus
            inputProps={{ step: 'any' }}
          />
          <TextField
            label="Poor value"
            value={poorStr}
            onChange={(e) => setPoorStr(e.target.value)}
            type="number"
            size="small"
            fullWidth
            inputProps={{ step: 'any' }}
          />
          <FormControl size="small" sx={{ minWidth: 100 }}>
            <InputLabel>Scale</InputLabel>
            <Select
              value={scale}
              onChange={(e) => setScale(e.target.value as 'log' | 'linear')}
              label="Scale"
            >
              <MenuItem value="log">Log</MenuItem>
              <MenuItem value="linear">Linear</MenuItem>
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
          disabled={saving}
          startIcon={saving ? <CircularProgress size={16} /> : null}
        >
          {saving ? 'Saving…' : 'Save'}
        </Button>
      </DialogActions>
    </Dialog>
  );
}
