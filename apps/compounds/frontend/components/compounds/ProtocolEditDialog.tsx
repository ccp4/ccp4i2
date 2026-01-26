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
  Autocomplete,
} from '@mui/material';
import { Close, Edit, Code, Add, Update } from '@mui/icons-material';
import { useSWRConfig } from 'swr';
import { TightBindingParametersForm } from './TightBindingParametersForm';
import { FittingMethodEditDialog } from './FittingMethodEditDialog';
import { DilutionSeriesCreateDialog } from './DilutionSeriesCreateDialog';
import { useCompoundsApi } from '@/lib/compounds/api';
import { useAuth } from '@/lib/compounds/auth-context';
import type { Protocol, FittingMethod, FittingParameters, DilutionSeries, AnalysisMethod, Target } from '@/types/compounds/models';

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
  const { canAdminister } = useAuth();
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [scriptEditorOpen, setScriptEditorOpen] = useState(false);
  const [createDilutionsDialogOpen, setCreateDilutionsDialogOpen] = useState(false);

  // Target propagation confirmation dialog state
  const [showPropagateDialog, setShowPropagateDialog] = useState(false);
  const [propagating, setPropagating] = useState(false);
  const [propagateResult, setPropagateResult] = useState<{ updated: number } | null>(null);

  // Dilution series propagation dialog state (admin only)
  const [showDilutionPropagateDialog, setShowDilutionPropagateDialog] = useState(false);
  const [dilutionPropagating, setDilutionPropagating] = useState(false);
  const [dilutionPropagateResult, setDilutionPropagateResult] = useState<{
    updated_series: number;
    deleted_analyses: number;
  } | null>(null);

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
  const [targetId, setTargetId] = useState<string | null>(protocol.target || null);
  const [originalTargetId] = useState<string | null>(protocol.target || null);
  const [originalDilutionsId] = useState<string | null>(protocol.preferred_dilutions || null);
  const [comments, setComments] = useState(protocol.comments || '');

  // Fetch available fitting methods
  const { data: fittingMethods, isLoading: methodsLoading } = api.get<FittingMethod[]>(
    'fitting-methods/'
  );

  // Fetch available dilution series
  const { data: dilutionSeries, isLoading: dilutionsLoading } = api.get<DilutionSeries[]>(
    'dilution-series/'
  );

  // Fetch available targets
  const { data: targets } = api.get<Target[]>('targets/');

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
    setTargetId(protocol.target || null);
    setComments(protocol.comments || '');
    setError(null);
    setPropagateResult(null);
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
        target: targetId || null,
        comments: comments || null,
      });

      // If target changed and new target is set, offer to propagate to assays
      const targetChanged = targetId !== originalTargetId;
      const hasAssays = (protocol.assays_count || 0) > 0;

      if (targetChanged && targetId && hasAssays) {
        setShowPropagateDialog(true);
        return;
      }

      // If dilution series changed and user is admin, offer to propagate to data series
      const dilutionsChanged = preferredDilutionsId !== originalDilutionsId;

      if (dilutionsChanged && preferredDilutionsId && hasAssays && canAdminister) {
        setShowDilutionPropagateDialog(true);
        return;
      }

      onSave();
      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to save');
    } finally {
      setSaving(false);
    }
  };

  const handlePropagate = async () => {
    setPropagating(true);
    setError(null);

    try {
      const result = await api.post<{ updated: number }>(`protocols/${protocol.id}/propagate_target/`, {});
      setPropagateResult(result);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to update assays');
    } finally {
      setPropagating(false);
    }
  };

  const handleClosePropagateDialog = () => {
    setShowPropagateDialog(false);
    setPropagateResult(null);

    // Check if dilution series also changed (for admin users)
    const dilutionsChanged = preferredDilutionsId !== originalDilutionsId;
    const hasAssays = (protocol.assays_count || 0) > 0;

    if (dilutionsChanged && preferredDilutionsId && hasAssays && canAdminister) {
      setShowDilutionPropagateDialog(true);
      return;
    }

    onSave();
    onClose();
  };

  const handleDilutionPropagate = async () => {
    setDilutionPropagating(true);
    setError(null);

    try {
      const result = await api.post<{
        updated_series: number;
        deleted_analyses: number;
      }>(`protocols/${protocol.id}/propagate_dilutions/`, {});
      setDilutionPropagateResult(result);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to update data series');
    } finally {
      setDilutionPropagating(false);
    }
  };

  const handleCloseDilutionPropagateDialog = () => {
    setShowDilutionPropagateDialog(false);
    setDilutionPropagateResult(null);
    onSave();
    onClose();
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

          {/* Default Target */}
          <Autocomplete
            options={targets || []}
            getOptionLabel={(option) => option.name}
            value={targets?.find((t) => t.id === targetId) || null}
            onChange={(_, newValue) => setTargetId(newValue?.id || null)}
            size="small"
            renderInput={(params) => (
              <TextField
                {...params}
                label="Default Target"
                helperText="Default target for assays using this protocol"
              />
            )}
          />

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

      {/* Propagate Target Confirmation Dialog */}
      <Dialog open={showPropagateDialog} onClose={handleClosePropagateDialog} maxWidth="xs" fullWidth>
        <DialogTitle sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Update color="primary" />
          Update Existing Assays?
        </DialogTitle>
        <DialogContent>
          {propagateResult ? (
            <Alert severity="success" sx={{ mt: 1 }}>
              Updated {propagateResult.updated} assay{propagateResult.updated !== 1 ? 's' : ''} to use the new target.
            </Alert>
          ) : (
            <>
              <Typography variant="body2" gutterBottom>
                You changed the default target for this protocol. Would you like to update
                all existing assays ({protocol.assays_count}) to use this target as well?
              </Typography>
              <Typography variant="body2" color="text.secondary" sx={{ mt: 1 }}>
                New target: <strong>{targets?.find((t) => t.id === targetId)?.name}</strong>
              </Typography>
            </>
          )}
          {error && (
            <Alert severity="error" sx={{ mt: 2 }}>
              {error}
            </Alert>
          )}
        </DialogContent>
        <DialogActions>
          {propagateResult ? (
            <Button variant="contained" onClick={handleClosePropagateDialog}>
              Done
            </Button>
          ) : (
            <>
              <Button onClick={handleClosePropagateDialog} disabled={propagating}>
                No, Keep Existing
              </Button>
              <Button
                variant="contained"
                onClick={handlePropagate}
                disabled={propagating}
                startIcon={propagating ? <CircularProgress size={16} /> : <Update />}
              >
                {propagating ? 'Updating...' : 'Yes, Update All'}
              </Button>
            </>
          )}
        </DialogActions>
      </Dialog>

      {/* Propagate Dilution Series Confirmation Dialog (Admin Only) */}
      <Dialog open={showDilutionPropagateDialog} onClose={handleCloseDilutionPropagateDialog} maxWidth="sm" fullWidth>
        <DialogTitle sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Update color="warning" />
          Update Data Series Dilutions?
        </DialogTitle>
        <DialogContent>
          {dilutionPropagateResult ? (
            <Alert severity="success" sx={{ mt: 1 }}>
              Updated {dilutionPropagateResult.updated_series} data series.
              {dilutionPropagateResult.deleted_analyses > 0 && (
                <> Deleted {dilutionPropagateResult.deleted_analyses} analysis result{dilutionPropagateResult.deleted_analyses !== 1 ? 's' : ''}.</>
              )}
            </Alert>
          ) : (
            <>
              <Typography variant="body2" gutterBottom>
                You changed the preferred dilution series for this protocol. Would you like to update
                all data series in assays using this protocol?
              </Typography>
              <Alert severity="warning" sx={{ mt: 2 }}>
                <strong>Warning:</strong> This will delete all existing analysis results for affected
                data series. They will need to be re-analyzed with the new dilution series.
              </Alert>
              <Typography variant="body2" color="text.secondary" sx={{ mt: 2 }}>
                New dilution series: <strong>{dilutionSeries?.find((ds) => ds.id === preferredDilutionsId)?.display_name || 'Selected series'}</strong>
              </Typography>
              <Typography variant="body2" color="text.secondary" sx={{ mt: 1 }}>
                Affected assays: <strong>{protocol.assays_count}</strong>
              </Typography>
            </>
          )}
          {error && (
            <Alert severity="error" sx={{ mt: 2 }}>
              {error}
            </Alert>
          )}
        </DialogContent>
        <DialogActions>
          {dilutionPropagateResult ? (
            <Button variant="contained" onClick={handleCloseDilutionPropagateDialog}>
              Done
            </Button>
          ) : (
            <>
              <Button onClick={handleCloseDilutionPropagateDialog} disabled={dilutionPropagating}>
                No, Keep Existing
              </Button>
              <Button
                variant="contained"
                color="warning"
                onClick={handleDilutionPropagate}
                disabled={dilutionPropagating}
                startIcon={dilutionPropagating ? <CircularProgress size={16} /> : <Update />}
              >
                {dilutionPropagating ? 'Updating...' : 'Yes, Update All'}
              </Button>
            </>
          )}
        </DialogActions>
      </Dialog>
    </Dialog>
  );
}
