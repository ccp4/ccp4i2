/*
 * Copyright (C) 2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
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
} from '@mui/material';
import { Close, Science, Add } from '@mui/icons-material';
import { useCompoundsApi, apiPost } from '@/lib/api';
import type { Target } from '@/types/models';

interface TargetCreateDialogProps {
  open: boolean;
  onClose: () => void;
  onCreated: (target: Target) => void;
}

export function TargetCreateDialog({
  open,
  onClose,
  onCreated,
}: TargetCreateDialogProps) {
  const api = useCompoundsApi();
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Form state
  const [name, setName] = useState('');
  const [parentId, setParentId] = useState<string | null>(null);

  // Fetch existing targets for parent selection
  const { data: targets, isLoading: targetsLoading } = api.get<Target[]>('targets/');

  // Reset form when dialog opens
  useEffect(() => {
    if (open) {
      setName('');
      setParentId(null);
      setError(null);
    }
  }, [open]);

  const handleSave = async () => {
    if (!name.trim()) {
      setError('Target name is required');
      return;
    }

    setSaving(true);
    setError(null);

    try {
      const newTarget = await apiPost<Target>('targets/', {
        name: name.trim(),
        parent: parentId || null,
      });

      onCreated(newTarget);
      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to create target');
    } finally {
      setSaving(false);
    }
  };

  return (
    <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
      <DialogTitle sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Science color="primary" />
          New Target
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
            label="Target Name"
            value={name}
            onChange={(e) => setName(e.target.value)}
            fullWidth
            required
            autoFocus
            size="small"
            placeholder="e.g., EGFR, BCR-ABL, PD-L1"
          />

          <FormControl fullWidth size="small">
            <InputLabel>Parent Target (optional)</InputLabel>
            <Select
              value={parentId || ''}
              onChange={(e) => setParentId(e.target.value || null)}
              label="Parent Target (optional)"
            >
              <MenuItem value="">
                <em>None</em>
              </MenuItem>
              {targetsLoading ? (
                <MenuItem disabled>Loading...</MenuItem>
              ) : (
                targets?.map((target) => (
                  <MenuItem key={target.id} value={target.id}>
                    {target.name}
                  </MenuItem>
                ))
              )}
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
          disabled={saving || !name.trim()}
          startIcon={saving ? <CircularProgress size={16} /> : <Add />}
        >
          {saving ? 'Creating...' : 'Create Target'}
        </Button>
      </DialogActions>
    </Dialog>
  );
}
