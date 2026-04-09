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
  Collapse,
} from '@mui/material';
import { Close, Add, GridOn, ExpandMore, ExpandLess } from '@mui/icons-material';
import { apiPost } from '@/lib/compounds/api';
import { PlateLayoutEditor } from './PlateLayoutEditor';
import type { PlateLayout, PlateLayoutRecord } from '@/types/compounds/models';

interface PlateLayoutCreateDialogProps {
  open: boolean;
  onClose: () => void;
  onCreated: (layout: PlateLayoutRecord) => void;
}

/**
 * Default plate layout for a 384-well plate
 */
const DEFAULT_LAYOUT: PlateLayout = {
  plate_format: 384,
  controls: {
    placement: 'edge_columns',
    max: { columns: [1, 2], rows: ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P'] },
    min: { columns: [23, 24], rows: ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P'] },
  },
  sample_region: {
    start_column: 3,
    end_column: 22,
    start_row: 'A',
    end_row: 'P',
  },
  dilution: {
    direction: 'horizontal',
    num_concentrations: 20,
  },
  replicate: {
    count: 2,
    pattern: 'adjacent_rows',
  },
  compound_source: {
    type: 'row_order',
  },
  spreadsheet_origin: {
    column: 'A',
    row: 1,
  },
};

export function PlateLayoutCreateDialog({
  open,
  onClose,
  onCreated,
}: PlateLayoutCreateDialogProps) {
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Form state
  const [name, setName] = useState('');
  const [description, setDescription] = useState('');
  const [config, setConfig] = useState<PlateLayout>(DEFAULT_LAYOUT);
  const [showEditor, setShowEditor] = useState(false);

  // Reset form when dialog opens
  useEffect(() => {
    if (open) {
      setName('');
      setDescription('');
      setConfig(DEFAULT_LAYOUT);
      setShowEditor(false);
      setError(null);
    }
  }, [open]);

  const handleSave = async () => {
    if (!name.trim()) {
      setError('Please enter a name for the plate layout');
      return;
    }

    setSaving(true);
    setError(null);

    try {
      const newLayout = await apiPost<PlateLayoutRecord>('plate-layouts/', {
        name: name.trim(),
        description: description.trim() || null,
        config,
      });

      onCreated(newLayout);
      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to create plate layout');
    } finally {
      setSaving(false);
    }
  };

  return (
    <Dialog
      open={open}
      onClose={onClose}
      maxWidth="lg"
      fullWidth
      PaperProps={{
        sx: { maxHeight: '90vh' }
      }}
    >
      <DialogTitle sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <GridOn color="primary" />
          New Plate Layout
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
            label="Name"
            value={name}
            onChange={(e) => setName(e.target.value)}
            fullWidth
            required
            autoFocus
            size="small"
            placeholder="e.g., Standard 384-well (20-point)"
            helperText="A unique name for this plate layout"
          />

          <TextField
            label="Description"
            value={description}
            onChange={(e) => setDescription(e.target.value)}
            fullWidth
            multiline
            rows={2}
            size="small"
            placeholder="Optional description of when to use this layout"
          />

          <Button
            variant="text"
            onClick={() => setShowEditor(!showEditor)}
            endIcon={showEditor ? <ExpandLess /> : <ExpandMore />}
            sx={{ alignSelf: 'flex-start' }}
          >
            {showEditor ? 'Hide' : 'Show'} Plate Configuration
          </Button>

          <Collapse in={showEditor}>
            <Box sx={{ border: 1, borderColor: 'divider', borderRadius: 1, p: 2 }}>
              <PlateLayoutEditor
                value={config}
                onChange={setConfig}
              />
            </Box>
          </Collapse>
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
          {saving ? 'Creating...' : 'Create Plate Layout'}
        </Button>
      </DialogActions>
    </Dialog>
  );
}
