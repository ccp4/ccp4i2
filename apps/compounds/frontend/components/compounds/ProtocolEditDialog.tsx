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
  IconButton,
  Divider,
  CircularProgress,
} from '@mui/material';
import { Close, Edit, Code } from '@mui/icons-material';
import { TightBindingParametersForm } from './TightBindingParametersForm';
import { FittingMethodEditDialog } from './FittingMethodEditDialog';
import { useCompoundsApi } from '@/lib/compounds/api';
import type { Protocol, FittingMethod, FittingParameters, DilutionSeries } from '@/types/compounds/models';

interface ProtocolEditDialogProps {
  open: boolean;
  onClose: () => void;
  protocol: Protocol;
  onSave: () => void;
}

export function ProtocolEditDialog({
  open,
  onClose,
  protocol,
  onSave,
}: ProtocolEditDialogProps) {
  const api = useCompoundsApi();
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [scriptEditorOpen, setScriptEditorOpen] = useState(false);

  // Form state
  const [name, setName] = useState(protocol.name);
  const [fittingMethodId, setFittingMethodId] = useState<string | null>(
    protocol.fitting_method || null
  );
  const [fittingParams, setFittingParams] = useState<FittingParameters>(
    protocol.fitting_parameters || {}
  );
  const [preferredDilutionsId, setPreferredDilutionsId] = useState<string | null>(
    protocol.preferred_dilutions || null
  );
  const [pherastarTable, setPherastarTable] = useState(protocol.pherastar_table || '');
  const [comments, setComments] = useState(protocol.comments || '');

  // Fetch available fitting methods
  const { data: fittingMethods, isLoading: methodsLoading } = api.get<FittingMethod[]>(
    'fitting-methods/'
  );

  // Fetch available dilution series
  const { data: dilutionSeries, isLoading: dilutionsLoading } = api.get<DilutionSeries[]>(
    'dilution-series/'
  );

  // Reset form when protocol changes
  useEffect(() => {
    setName(protocol.name);
    setFittingMethodId(protocol.fitting_method || null);
    setFittingParams(protocol.fitting_parameters || {});
    setPreferredDilutionsId(protocol.preferred_dilutions || null);
    setPherastarTable(protocol.pherastar_table || '');
    setComments(protocol.comments || '');
    setError(null);
  }, [protocol]);

  // Check if selected fitting method is tight-binding
  const selectedMethod = fittingMethods?.find((m) => m.id === fittingMethodId);
  const isTightBinding =
    selectedMethod?.slug?.includes('tight-binding') ||
    selectedMethod?.name?.toLowerCase().includes('wang');

  // Get the unit from selected dilution series for the parameters form
  const selectedDilutions = dilutionSeries?.find((ds) => ds.id === preferredDilutionsId);
  const concentrationUnit = selectedDilutions?.unit || 'nM';

  const handleSave = async () => {
    setSaving(true);
    setError(null);

    try {
      const response = await fetch(`/api/proxy/compounds/protocols/${protocol.id}/`, {
        method: 'PATCH',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          name,
          fitting_method: fittingMethodId || null,
          fitting_parameters: fittingParams || {},
          preferred_dilutions: preferredDilutionsId || null,
          pherastar_table: pherastarTable || null,
          comments: comments || null,
        }),
      });

      if (!response.ok) {
        const data = await response.json().catch(() => ({}));
        // DRF returns validation errors as {field_name: ["error"]}
        const errorMessages: string[] = [];
        for (const [field, errors] of Object.entries(data)) {
          if (Array.isArray(errors)) {
            errorMessages.push(`${field}: ${errors.join(', ')}`);
          } else if (typeof errors === 'string') {
            errorMessages.push(`${field}: ${errors}`);
          }
        }
        throw new Error(errorMessages.length > 0 ? errorMessages.join('; ') : 'Failed to save protocol');
      }

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
          Edit Protocol
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
          {/* Protocol Name */}
          <TextField
            label="Protocol Name"
            value={name}
            onChange={(e) => setName(e.target.value)}
            fullWidth
            required
            size="small"
          />

          <Divider />

          {/* Fitting Method */}
          <Typography variant="subtitle2" color="text.secondary">
            Analysis Settings
          </Typography>

          <Box sx={{ display: 'flex', gap: 1, alignItems: 'flex-start' }}>
            <FormControl fullWidth size="small">
              <InputLabel>Fitting Method</InputLabel>
              <Select
                value={fittingMethodId || ''}
                onChange={(e) => {
                  setFittingMethodId(e.target.value || null);
                  // Clear fitting params when method changes
                  setFittingParams({});
                }}
                label="Fitting Method"
              >
                <MenuItem value="">
                  <em>Default (4PL)</em>
                </MenuItem>
                {methodsLoading ? (
                  <MenuItem disabled>Loading...</MenuItem>
                ) : (
                  fittingMethods?.map((method) => (
                    <MenuItem key={method.id} value={method.id}>
                      {method.name}
                    </MenuItem>
                  ))
                )}
              </Select>
            </FormControl>
            {fittingMethodId && (
              <Button
                variant="outlined"
                size="small"
                startIcon={<Code />}
                onClick={() => setScriptEditorOpen(true)}
                sx={{ minWidth: 'auto', whiteSpace: 'nowrap', mt: 0.5 }}
              >
                Edit Script
              </Button>
            )}
          </Box>

          {/* Preferred Dilutions */}
          <FormControl fullWidth size="small">
            <InputLabel>Preferred Dilutions</InputLabel>
            <Select
              value={preferredDilutionsId || ''}
              onChange={(e) => setPreferredDilutionsId(e.target.value || null)}
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
                    {ds.display_name || `${ds.concentrations.length} points (${ds.unit})`}
                  </MenuItem>
                ))
              )}
            </Select>
          </FormControl>

          {/* Tight-binding parameters (conditional) */}
          {isTightBinding && (
            <Box sx={{ pl: 1, borderLeft: 3, borderColor: 'secondary.main' }}>
              <TightBindingParametersForm
                value={fittingParams}
                onChange={setFittingParams}
                unit={concentrationUnit}
              />
            </Box>
          )}

          <Divider />

          {/* Other Settings */}
          <Typography variant="subtitle2" color="text.secondary">
            Other Settings
          </Typography>

          <TextField
            label="PHERAstar Table"
            value={pherastarTable}
            onChange={(e) => setPherastarTable(e.target.value)}
            fullWidth
            size="small"
            placeholder="e.g., FI 485 520"
          />

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
          disabled={saving || !name.trim()}
          startIcon={saving ? <CircularProgress size={16} /> : null}
        >
          {saving ? 'Saving...' : 'Save Changes'}
        </Button>
      </DialogActions>

      {/* Fitting Method Script Editor */}
      {fittingMethodId && (
        <FittingMethodEditDialog
          open={scriptEditorOpen}
          onClose={() => setScriptEditorOpen(false)}
          fittingMethodId={fittingMethodId}
        />
      )}
    </Dialog>
  );
}
