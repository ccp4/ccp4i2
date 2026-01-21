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
  Divider,
  ListItemIcon,
  ListItemText,
} from '@mui/material';
import { Close, Description, Add } from '@mui/icons-material';
import { useSWRConfig } from 'swr';
import { DilutionSeriesCreateDialog } from './DilutionSeriesCreateDialog';
import { useCompoundsApi, apiPost } from '@/lib/compounds/api';
import type { Protocol, DilutionSeries, AnalysisMethod } from '@/types/compounds/models';

interface ProtocolCreateDialogProps {
  open: boolean;
  onClose: () => void;
  onCreated: (protocol: Protocol) => void;
}

const ANALYSIS_METHOD_OPTIONS: { value: AnalysisMethod; label: string }[] = [
  { value: 'hill_langmuir', label: 'Hill-Langmuir' },
  { value: 'hill_langmuir_fix_hill', label: 'Hill-Langmuir (fixed Hill)' },
  { value: 'hill_langmuir_fix_hill_minmax', label: 'Hill-Langmuir (fixed Hill/min/max)' },
  { value: 'hill_langmuir_fix_minmax', label: 'Hill-Langmuir (fixed min/max)' },
  { value: 'ms_intact', label: 'MS-Intact' },
  { value: 'table_of_values', label: 'Table of values' },
  { value: 'pharmaron_adme', label: 'Pharmaron ADME' },
];

export function ProtocolCreateDialog({
  open,
  onClose,
  onCreated,
}: ProtocolCreateDialogProps) {
  const api = useCompoundsApi();
  const { mutate } = useSWRConfig();
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [createDilutionsDialogOpen, setCreateDilutionsDialogOpen] = useState(false);

  // Form state
  const [name, setName] = useState('');
  const [analysisMethod, setAnalysisMethod] = useState<AnalysisMethod>('hill_langmuir');
  const [preferredDilutionsId, setPreferredDilutionsId] = useState<string | null>(null);
  const [comments, setComments] = useState('');

  // Fetch available dilution series
  const { data: dilutionSeries, isLoading: dilutionsLoading } = api.get<DilutionSeries[]>(
    'dilution-series/'
  );

  const handleDilutionSeriesCreated = (newSeries: DilutionSeries) => {
    mutate('/api/proxy/compounds/dilution-series/');
    setPreferredDilutionsId(newSeries.id);
  };

  // Reset form when dialog opens
  useEffect(() => {
    if (open) {
      setName('');
      setAnalysisMethod('hill_langmuir');
      setPreferredDilutionsId(null);
      setComments('');
      setError(null);
    }
  }, [open]);

  const handleSave = async () => {
    if (!name.trim()) {
      setError('Protocol name is required');
      return;
    }

    setSaving(true);
    setError(null);

    try {
      const newProtocol = await apiPost<Protocol>('protocols/', {
        name: name.trim(),
        analysis_method: analysisMethod,
        preferred_dilutions: preferredDilutionsId || null,
        comments: comments.trim() || null,
      });

      onCreated(newProtocol);
      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to create protocol');
    } finally {
      setSaving(false);
    }
  };

  return (
    <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
      <DialogTitle sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Description color="primary" />
          New Protocol
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
            label="Protocol Name"
            value={name}
            onChange={(e) => setName(e.target.value)}
            fullWidth
            required
            autoFocus
            size="small"
            placeholder="e.g., FP Competition Assay"
          />

          <FormControl fullWidth size="small">
            <InputLabel>Analysis Method</InputLabel>
            <Select
              value={analysisMethod}
              onChange={(e) => setAnalysisMethod(e.target.value as AnalysisMethod)}
              label="Analysis Method"
            >
              {ANALYSIS_METHOD_OPTIONS.map((option) => (
                <MenuItem key={option.value} value={option.value}>
                  {option.label}
                </MenuItem>
              ))}
            </Select>
          </FormControl>

          {/* Hide dilutions for ADME protocols - they use time points instead */}
          {analysisMethod !== 'pharmaron_adme' && (
            <FormControl fullWidth size="small">
              <InputLabel>Preferred Dilutions</InputLabel>
              <Select
                value={preferredDilutionsId || ''}
                onChange={(e) => {
                  const value = e.target.value;
                  if (value === '__create_new__') {
                    setCreateDilutionsDialogOpen(true);
                  } else {
                    setPreferredDilutionsId(value || null);
                  }
                }}
                label="Preferred Dilutions"
              >
                <MenuItem value="">
                  <em>None</em>
                </MenuItem>
                {dilutionsLoading ? (
                  <MenuItem disabled>Loading...</MenuItem>
                ) : (
                  dilutionSeries?.map((ds) => (
                    <MenuItem key={ds.id} value={ds.id}>
                      {ds.display_name || `${ds.concentrations.join(', ')} ${ds.unit}`}
                    </MenuItem>
                  ))
                )}
                <Divider />
                <MenuItem value="__create_new__">
                  <ListItemIcon>
                    <Add fontSize="small" />
                  </ListItemIcon>
                  <ListItemText>Create New...</ListItemText>
                </MenuItem>
              </Select>
            </FormControl>
          )}

          <TextField
            label="Comments"
            value={comments}
            onChange={(e) => setComments(e.target.value)}
            fullWidth
            multiline
            rows={3}
            size="small"
            placeholder="Optional notes about this protocol"
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
          {saving ? 'Creating...' : 'Create Protocol'}
        </Button>
      </DialogActions>

      <DilutionSeriesCreateDialog
        open={createDilutionsDialogOpen}
        onClose={() => setCreateDilutionsDialogOpen(false)}
        onCreated={handleDilutionSeriesCreated}
      />
    </Dialog>
  );
}
