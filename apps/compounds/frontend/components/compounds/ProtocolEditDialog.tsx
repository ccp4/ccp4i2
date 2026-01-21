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
  ListItemIcon,
  ListItemText,
} from '@mui/material';
import { Close, Edit, Code, Add } from '@mui/icons-material';
import { useSWRConfig } from 'swr';
import { TightBindingParametersForm } from './TightBindingParametersForm';
import { FittingMethodEditDialog } from './FittingMethodEditDialog';
import { DilutionSeriesCreateDialog } from './DilutionSeriesCreateDialog';
import { useCompoundsApi } from '@/lib/compounds/api';
import type { Protocol, FittingMethod, FittingParameters, DilutionSeries, AnalysisMethod } from '@/types/compounds/models';

const ANALYSIS_METHOD_OPTIONS: { value: AnalysisMethod; label: string }[] = [
  { value: 'hill_langmuir', label: 'Hill-Langmuir' },
  { value: 'hill_langmuir_fix_hill', label: 'Hill-Langmuir (fixed Hill)' },
  { value: 'hill_langmuir_fix_hill_minmax', label: 'Hill-Langmuir (fixed Hill/min/max)' },
  { value: 'hill_langmuir_fix_minmax', label: 'Hill-Langmuir (fixed min/max)' },
  { value: 'ms_intact', label: 'MS-Intact' },
  { value: 'table_of_values', label: 'Table of values' },
  { value: 'pharmaron_adme', label: 'Pharmaron ADME' },
];

// Helper to check if analysis method is Pharmaron ADME (handles both underscore and hyphen variants)
function isAdmeProtocol(analysisMethod: string | undefined): boolean {
  if (!analysisMethod) return false;
  const normalized = analysisMethod.toLowerCase().replace(/-/g, '_');
  return normalized === 'pharmaron_adme' || normalized === 'adme';
}

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
  const { mutate } = useSWRConfig();
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [scriptEditorOpen, setScriptEditorOpen] = useState(false);
  const [createDilutionsDialogOpen, setCreateDilutionsDialogOpen] = useState(false);

  // Form state
  const [name, setName] = useState(protocol.name);
  const [analysisMethod, setAnalysisMethod] = useState<AnalysisMethod>(
    protocol.analysis_method as AnalysisMethod || 'hill_langmuir'
  );
  const [fittingMethodId, setFittingMethodId] = useState<string | null>(
    protocol.fitting_method || null
  );
  const [fittingParams, setFittingParams] = useState<FittingParameters>(
    protocol.fitting_parameters || {}
  );
  const [preferredDilutionsId, setPreferredDilutionsId] = useState<string | null>(
    protocol.preferred_dilutions || null
  );
  const [comments, setComments] = useState(protocol.comments || '');

  // Fetch available fitting methods
  const { data: fittingMethods, isLoading: methodsLoading } = api.get<FittingMethod[]>(
    'fitting-methods/'
  );

  // Fetch available dilution series
  const { data: dilutionSeries, isLoading: dilutionsLoading } = api.get<DilutionSeries[]>(
    'dilution-series/'
  );

  const handleDilutionSeriesCreated = (newSeries: DilutionSeries) => {
    mutate('/api/proxy/compounds/dilution-series/');
    setPreferredDilutionsId(newSeries.id);
  };

  // Reset form when protocol changes
  useEffect(() => {
    setName(protocol.name);
    setAnalysisMethod(protocol.analysis_method as AnalysisMethod || 'hill_langmuir');
    setFittingMethodId(protocol.fitting_method || null);
    setFittingParams(protocol.fitting_parameters || {});
    setPreferredDilutionsId(protocol.preferred_dilutions || null);
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
      await api.patch(`protocols/${protocol.id}/`, {
        name,
        analysis_method: analysisMethod,
        fitting_method: fittingMethodId || null,
        fitting_parameters: fittingParams || {},
        preferred_dilutions: preferredDilutionsId || null,
        comments: comments || null,
      });
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

          {/* Analysis Method */}
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

          {/* Analysis Settings - hidden for ADME protocols */}
          {!isAdmeProtocol(analysisMethod) && (
            <>
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
                        {ds.display_name || `${ds.concentrations.length} points (${ds.unit})`}
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
            </>
          )}

          {/* ADME-specific info */}
          {isAdmeProtocol(analysisMethod) && (
            <>
              <Divider />
              <Alert severity="info" sx={{ mt: 1 }}>
                This is a Pharmaron ADME protocol. Data is imported from Excel files
                (ADME-NCU-*.xlsx) rather than plate-based assays.
              </Alert>
            </>
          )}

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

      <DilutionSeriesCreateDialog
        open={createDilutionsDialogOpen}
        onClose={() => setCreateDilutionsDialogOpen(false)}
        onCreated={handleDilutionSeriesCreated}
      />
    </Dialog>
  );
}
