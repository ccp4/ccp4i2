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
  Typography,
} from '@mui/material';
import { Close, Description, Add, GridOn } from '@mui/icons-material';
import { useSWRConfig } from 'swr';
import { DilutionSeriesCreateDialog } from './DilutionSeriesCreateDialog';
import { PlateLayoutCreateDialog } from './PlateLayoutCreateDialog';
import { FourPLConstraintsForm } from './FourPLConstraintsForm';
import { useCompoundsApi, apiPost } from '@/lib/compounds/api';
import type { Protocol, DilutionSeries, ImportType, FittingParameters, PlateLayoutRecord } from '@/types/compounds/models';

interface ProtocolCreateDialogProps {
  open: boolean;
  onClose: () => void;
  onCreated: (protocol: Protocol) => void;
}

const IMPORT_TYPE_OPTIONS: { value: ImportType; label: string; description?: string }[] = [
  { value: 'raw_data', label: 'Raw Data (Dose-Response)', description: 'Plate-based assay data for curve fitting' },
  { value: 'ms_intact', label: 'MS-Intact', description: 'Pre-analyzed MS-Intact data' },
  { value: 'table_of_values', label: 'Table of Values', description: 'Pre-analyzed data from external tools' },
  { value: 'pharmaron_adme', label: 'Pharmaron ADME', description: 'Pharmaron/NCU ADME Excel imports' },
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
  const [createLayoutDialogOpen, setCreateLayoutDialogOpen] = useState(false);

  // Form state
  const [name, setName] = useState('');
  const [importType, setImportType] = useState<ImportType>('raw_data');
  const [fittingParams, setFittingParams] = useState<FittingParameters>({});
  const [preferredDilutionsId, setPreferredDilutionsId] = useState<string | null>(null);
  const [plateLayoutId, setPlateLayoutId] = useState<string | null>(null);
  const [comments, setComments] = useState('');
  // Interpretation thresholds — recommended but not required at creation.
  // Direction is implied by ordering (target < poor → lower-better).
  const [targetValueStr, setTargetValueStr] = useState('');
  const [poorValueStr, setPoorValueStr] = useState('');
  const [thresholdScale, setThresholdScale] = useState<'log' | 'linear'>('log');

  // Fetch available dilution series
  const { data: dilutionSeries, isLoading: dilutionsLoading } = api.get<DilutionSeries[]>(
    'dilution-series/'
  );

  // Fetch available plate layouts
  const { data: plateLayouts, isLoading: layoutsLoading } = api.get<PlateLayoutRecord[]>(
    'plate-layouts/'
  );

  const handleDilutionSeriesCreated = (newSeries: DilutionSeries) => {
    mutate('/api/proxy/compounds/dilution-series/');
    setPreferredDilutionsId(newSeries.id);
  };

  const handlePlateLayoutCreated = (newLayout: PlateLayoutRecord) => {
    mutate('/api/proxy/compounds/plate-layouts/');
    setPlateLayoutId(newLayout.id);
  };

  // Reset form when dialog opens
  useEffect(() => {
    if (open) {
      setName('');
      setImportType('raw_data');
      setFittingParams({});
      setPreferredDilutionsId(null);
      setPlateLayoutId(null);
      setComments('');
      setTargetValueStr('');
      setPoorValueStr('');
      setThresholdScale('log');
      setError(null);
    }
  }, [open]);

  const handleSave = async () => {
    if (!name.trim()) {
      setError('Protocol name is required');
      return;
    }

    // Parse thresholds (empty → null)
    const parseThreshold = (s: string): number | null => {
      if (s.trim() === '') return null;
      const n = Number(s);
      return Number.isFinite(n) ? n : null;
    };
    const targetValue = parseThreshold(targetValueStr);
    const poorValue = parseThreshold(poorValueStr);

    if (targetValue !== null && poorValue !== null && targetValue === poorValue) {
      setError('Excellent and poor values must differ — the direction of "better" is implied by their ordering.');
      return;
    }

    setSaving(true);
    setError(null);

    try {
      // Build fitting_parameters - only include if there are constraints
      const hasFittingConstraints = fittingParams.fix_hill != null ||
        fittingParams.fix_top != null ||
        fittingParams.fix_bottom != null;

      const newProtocol = await apiPost<Protocol>('protocols/', {
        name: name.trim(),
        import_type: importType,
        preferred_dilutions: preferredDilutionsId || null,
        plate_layout: plateLayoutId || null,
        fitting_parameters: hasFittingConstraints ? fittingParams : null,
        comments: comments.trim() || null,
        target_value: targetValue,
        poor_value: poorValue,
        threshold_scale: thresholdScale,
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
            <InputLabel>Import Type</InputLabel>
            <Select
              value={importType}
              onChange={(e) => {
                setImportType(e.target.value as ImportType);
                // Clear fitting params when switching away from raw_data
                if (e.target.value !== 'raw_data') {
                  setFittingParams({});
                }
              }}
              label="Import Type"
            >
              {IMPORT_TYPE_OPTIONS.map((option) => (
                <MenuItem key={option.value} value={option.value}>
                  {option.label}
                </MenuItem>
              ))}
            </Select>
          </FormControl>

          {/* Show curve fitting constraints for raw data protocols */}
          {importType === 'raw_data' && (
            <FourPLConstraintsForm
              value={fittingParams}
              onChange={setFittingParams}
            />
          )}

          {/* Hide dilutions for ADME protocols - they use time points instead */}
          {importType !== 'pharmaron_adme' && (
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

          <FormControl fullWidth size="small">
            <InputLabel>Plate Layout</InputLabel>
            <Select
              value={plateLayoutId || ''}
              onChange={(e) => {
                const value = e.target.value;
                if (value === '__create_new__') {
                  setCreateLayoutDialogOpen(true);
                } else {
                  setPlateLayoutId(value || null);
                }
              }}
              label="Plate Layout"
            >
              <MenuItem value="">
                <em>None</em>
              </MenuItem>
              {layoutsLoading ? (
                <MenuItem disabled>Loading...</MenuItem>
              ) : (
                plateLayouts?.map((layout) => (
                  <MenuItem key={layout.id} value={layout.id}>
                    <ListItemIcon>
                      <GridOn fontSize="small" />
                    </ListItemIcon>
                    <ListItemText>
                      {layout.name}
                      {layout.plate_format && ` (${layout.plate_format}-well)`}
                    </ListItemText>
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

          <Divider />

          {/* Interpretation Thresholds — recommended but not required */}
          <Box>
            <Typography variant="subtitle2" color="text.secondary" gutterBottom>
              Value Interpretation (recommended)
            </Typography>
            <Typography variant="caption" color="text.secondary" sx={{ display: 'block', mb: 1.5 }}>
              Absolute thresholds for colour-coding KPI values in aggregation views.
              Direction of &quot;better&quot; is implied by ordering: set excellent &lt; poor
              for potency-style metrics (IC50, EC50), or excellent &gt; poor for
              capacity-style metrics (solubility, permeability).
              Leave blank if you&apos;ll decide later — the protocol will render
              uncoloured in aggregation views until set.
            </Typography>
            <Box sx={{ display: 'flex', gap: 1.5 }}>
              <TextField
                label="Excellent value"
                value={targetValueStr}
                onChange={(e) => setTargetValueStr(e.target.value)}
                type="number"
                size="small"
                fullWidth
                inputProps={{ step: 'any' }}
              />
              <TextField
                label="Poor value"
                value={poorValueStr}
                onChange={(e) => setPoorValueStr(e.target.value)}
                type="number"
                size="small"
                fullWidth
                inputProps={{ step: 'any' }}
              />
              <FormControl size="small" sx={{ minWidth: 110 }}>
                <InputLabel>Scale</InputLabel>
                <Select
                  value={thresholdScale}
                  onChange={(e) => setThresholdScale(e.target.value as 'log' | 'linear')}
                  label="Scale"
                >
                  <MenuItem value="log">Log</MenuItem>
                  <MenuItem value="linear">Linear</MenuItem>
                </Select>
              </FormControl>
            </Box>
          </Box>

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

      <PlateLayoutCreateDialog
        open={createLayoutDialogOpen}
        onClose={() => setCreateLayoutDialogOpen(false)}
        onCreated={handlePlateLayoutCreated}
      />
    </Dialog>
  );
}
