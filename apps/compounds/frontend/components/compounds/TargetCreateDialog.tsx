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
  Autocomplete,
  Chip,
} from '@mui/material';
import { Close, Science, Add } from '@mui/icons-material';
import { useCompoundsApi, apiPost } from '@/lib/compounds/api';
import type { Target } from '@/types/compounds/models';

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
  const [geneSymbols, setGeneSymbols] = useState<string[]>([]);

  // Fetch existing targets for parent selection
  const { data: targets, isLoading: targetsLoading } = api.get<Target[]>('targets/');

  // Reset form when dialog opens
  useEffect(() => {
    if (open) {
      setName('');
      setParentId(null);
      setGeneSymbols([]);
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
      const payload: {
        name: string;
        parent: string | null;
        gene_symbols?: string[];
      } = {
        name: name.trim(),
        parent: parentId || null,
      };

      // Only include gene_symbols if the user typed at least one. Omitting
      // leaves the Target without linked Genes (still valid — annotator
      // may not know the mapping yet; the local backfill script covers that).
      const cleaned = geneSymbols
        .map((s) => s.trim().toUpperCase())
        .filter((s) => s.length > 0);
      if (cleaned.length > 0) {
        payload.gene_symbols = cleaned;
      }

      const newTarget = await apiPost<Target>('targets/', payload);

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

          {/* Gene symbols — free-form, chip-style. Each symbol creates or */}
          {/* links to a Gene row (HGNC-hydrated asynchronously). Unknown */}
          {/* symbols are accepted (useful for non-human targets / pre-HGNC */}
          {/* programmes); the hydrate_genes command later fills in aliases. */}
          <Autocomplete
            multiple
            freeSolo
            size="small"
            options={[]}
            value={geneSymbols}
            onChange={(_event, newValue) => {
              // newValue is a mix of strings (new entries) and whatever the
              // dropdown provides (none here). Normalise: split on commas too,
              // so a paste of "EGFR, ERBB2, MET" becomes three chips.
              const expanded: string[] = [];
              for (const item of newValue) {
                if (typeof item !== 'string') continue;
                for (const piece of item.split(/[\s,]+/)) {
                  const trimmed = piece.trim().toUpperCase();
                  if (trimmed && !expanded.includes(trimmed)) {
                    expanded.push(trimmed);
                  }
                }
              }
              setGeneSymbols(expanded);
            }}
            renderTags={(value, getTagProps) =>
              value.map((option, index) => {
                const tagProps = getTagProps({ index });
                return (
                  <Chip
                    variant="outlined"
                    label={option}
                    size="small"
                    {...tagProps}
                    key={option}
                  />
                );
              })
            }
            renderInput={(params) => (
              <TextField
                {...params}
                label="Gene Symbols (optional)"
                placeholder="e.g. EGFR, ERBB2 — press Enter or comma after each"
                helperText="HGNC-approved symbols. Case is normalised to uppercase. Aliases and names are filled in automatically after creation."
              />
            )}
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
          {saving ? 'Creating...' : 'Create Target'}
        </Button>
      </DialogActions>
    </Dialog>
  );
}
