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
  Typography,
  Alert,
  IconButton,
  CircularProgress,
  Autocomplete,
} from '@mui/material';
import { Close, Edit } from '@mui/icons-material';
import { useSWRConfig } from 'swr';
import { useCompoundsApi } from '@/lib/compounds/api';
import type { Assay, Target } from '@/types/compounds/models';

interface AssayEditDialogProps {
  open: boolean;
  onClose: () => void;
  assay: Assay;
  onSave: () => void;
}

export function AssayEditDialog({
  open,
  onClose,
  assay,
  onSave,
}: AssayEditDialogProps) {
  const api = useCompoundsApi();
  const { mutate } = useSWRConfig();
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Form state
  const [targetId, setTargetId] = useState<string | null>(assay.target || null);
  const [labbookNumber, setLabbookNumber] = useState<string>(
    assay.labbook_number?.toString() || ''
  );
  const [pageNumber, setPageNumber] = useState<string>(
    assay.page_number?.toString() || ''
  );
  const [comments, setComments] = useState(assay.comments || '');

  // Fetch available targets
  const { data: targets } = api.get<Target[]>('targets/');

  // Reset form when assay changes
  useEffect(() => {
    setTargetId(assay.target || null);
    setLabbookNumber(assay.labbook_number?.toString() || '');
    setPageNumber(assay.page_number?.toString() || '');
    setComments(assay.comments || '');
    setError(null);
  }, [assay]);

  const handleSave = async () => {
    setSaving(true);
    setError(null);

    try {
      await api.patch(`assays/${assay.id}/`, {
        target: targetId || null,
        labbook_number: labbookNumber ? parseInt(labbookNumber) : null,
        page_number: pageNumber ? parseInt(pageNumber) : null,
        comments: comments || null,
      });

      // Invalidate the assay cache
      mutate(`/api/proxy/compounds/assays/${assay.id}/`);

      onSave();
      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to save');
    } finally {
      setSaving(false);
    }
  };

  return (
    <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
      <DialogTitle sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Edit color="primary" />
          Edit Assay
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

        <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2.5 }}>
          {/* Assay file name (read-only info) */}
          <Box sx={{ p: 1.5, bgcolor: 'grey.100', borderRadius: 1 }}>
            <Typography variant="caption" color="text.secondary">
              Data File
            </Typography>
            <Typography variant="body2" fontWeight={500}>
              {assay.data_filename || 'No file'}
            </Typography>
          </Box>

          {/* Target */}
          <Autocomplete
            options={targets || []}
            getOptionLabel={(option) => option.name}
            value={targets?.find((t) => t.id === targetId) || null}
            onChange={(_, newValue) => setTargetId(newValue?.id || null)}
            size="small"
            renderInput={(params) => (
              <TextField
                {...params}
                label="Target"
                helperText="Target being tested in this assay"
              />
            )}
          />

          {/* Lab book number and page in a row */}
          <Box sx={{ display: 'flex', gap: 2 }}>
            <TextField
              label="Lab Book Number"
              value={labbookNumber}
              onChange={(e) => setLabbookNumber(e.target.value.replace(/\D/g, ''))}
              size="small"
              type="number"
              inputProps={{ min: 0 }}
              sx={{ flex: 1 }}
            />
            <TextField
              label="Page Number"
              value={pageNumber}
              onChange={(e) => setPageNumber(e.target.value.replace(/\D/g, ''))}
              size="small"
              type="number"
              inputProps={{ min: 0 }}
              sx={{ flex: 1 }}
            />
          </Box>

          {/* Comments */}
          <TextField
            label="Comments"
            value={comments}
            onChange={(e) => setComments(e.target.value)}
            fullWidth
            multiline
            rows={3}
            size="small"
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
          startIcon={saving ? <CircularProgress size={16} /> : null}
        >
          {saving ? 'Saving...' : 'Save Changes'}
        </Button>
      </DialogActions>
    </Dialog>
  );
}
