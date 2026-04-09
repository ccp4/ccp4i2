'use client';

import { useState, useEffect } from 'react';
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  Box,
  Alert,
  IconButton,
  CircularProgress,
  Typography,
  List,
  ListItem,
  ListItemButton,
  ListItemIcon,
  ListItemText,
  Radio,
  Chip,
  Divider,
} from '@mui/material';
import { Close, Science, Sync, Check } from '@mui/icons-material';
import { useCompoundsApi } from '@/lib/compounds/api';
import type { DilutionSeries } from '@/types/compounds/models';

interface DilutionSeriesSelectDialogProps {
  open: boolean;
  onClose: () => void;
  dataSeriesId: string;
  currentDilutionSeriesId?: string;
  protocolDilutionSeriesId?: string;
  onSave: () => void;
}

export function DilutionSeriesSelectDialog({
  open,
  onClose,
  dataSeriesId,
  currentDilutionSeriesId,
  protocolDilutionSeriesId,
  onSave,
}: DilutionSeriesSelectDialogProps) {
  const api = useCompoundsApi();
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [selectedId, setSelectedId] = useState<string | null>(null);

  // Fetch all available dilution series
  const { data: dilutionSeriesList, isLoading } = api.get<DilutionSeries[]>('dilution-series/');

  // Reset selection when dialog opens
  useEffect(() => {
    if (open) {
      setSelectedId(currentDilutionSeriesId || null);
      setError(null);
    }
  }, [open, currentDilutionSeriesId]);

  const handleSave = async () => {
    if (!selectedId) {
      setError('Please select a dilution series');
      return;
    }

    setSaving(true);
    setError(null);

    try {
      await api.patch(`data-series/${dataSeriesId}/`, {
        dilution_series: selectedId,
      });
      onSave();
      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to update dilution series');
    } finally {
      setSaving(false);
    }
  };

  const handleSyncFromProtocol = async () => {
    if (!protocolDilutionSeriesId) {
      setError('Protocol has no preferred dilution series configured');
      return;
    }

    setSaving(true);
    setError(null);

    try {
      await api.patch(`data-series/${dataSeriesId}/`, {
        dilution_series: protocolDilutionSeriesId,
      });
      onSave();
      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to sync dilution series');
    } finally {
      setSaving(false);
    }
  };

  const formatConcentrations = (series: DilutionSeries) => {
    if (!series.concentrations || series.concentrations.length === 0) {
      return 'No concentrations';
    }
    const concs = series.concentrations;
    if (concs.length <= 5) {
      return concs.join(', ') + ` ${series.unit}`;
    }
    return `${concs[0]} - ${concs[concs.length - 1]} ${series.unit} (${concs.length} points)`;
  };

  const hasChanged = selectedId !== currentDilutionSeriesId;

  return (
    <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
      <DialogTitle sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Science color="primary" />
          Select Dilution Series
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

        {/* Sync from Protocol option */}
        {protocolDilutionSeriesId && protocolDilutionSeriesId !== currentDilutionSeriesId && (
          <>
            <Box sx={{ mb: 2 }}>
              <Typography variant="subtitle2" color="text.secondary" gutterBottom>
                Quick Action
              </Typography>
              <Button
                variant="outlined"
                startIcon={saving ? <CircularProgress size={16} /> : <Sync />}
                onClick={handleSyncFromProtocol}
                disabled={saving}
                fullWidth
                sx={{ justifyContent: 'flex-start', textTransform: 'none' }}
              >
                Sync from Protocol's Preferred Dilutions
              </Button>
            </Box>
            <Divider sx={{ my: 2 }} />
          </>
        )}

        <Typography variant="subtitle2" color="text.secondary" gutterBottom>
          Available Dilution Series
        </Typography>

        {isLoading ? (
          <Box sx={{ display: 'flex', justifyContent: 'center', py: 4 }}>
            <CircularProgress />
          </Box>
        ) : dilutionSeriesList && dilutionSeriesList.length > 0 ? (
          <List sx={{ maxHeight: 300, overflow: 'auto' }}>
            {dilutionSeriesList.map((series) => {
              const isSelected = selectedId === series.id;
              const isCurrent = currentDilutionSeriesId === series.id;
              const isProtocol = protocolDilutionSeriesId === series.id;

              return (
                <ListItem key={series.id} disablePadding>
                  <ListItemButton
                    onClick={() => setSelectedId(series.id)}
                    selected={isSelected}
                    sx={{
                      borderRadius: 1,
                      mb: 0.5,
                      border: isSelected ? '2px solid' : '1px solid',
                      borderColor: isSelected ? 'primary.main' : 'divider',
                    }}
                  >
                    <ListItemIcon sx={{ minWidth: 40 }}>
                      <Radio checked={isSelected} size="small" />
                    </ListItemIcon>
                    <ListItemText
                      primary={
                        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                          {series.display_name || `Series ${series.id.slice(0, 8)}`}
                          {isCurrent && (
                            <Chip label="Current" size="small" color="info" variant="outlined" />
                          )}
                          {isProtocol && !isCurrent && (
                            <Chip label="Protocol" size="small" color="success" variant="outlined" />
                          )}
                        </Box>
                      }
                      secondary={formatConcentrations(series)}
                      secondaryTypographyProps={{ sx: { fontFamily: 'monospace', fontSize: '0.75rem' } }}
                    />
                  </ListItemButton>
                </ListItem>
              );
            })}
          </List>
        ) : (
          <Typography color="text.secondary" sx={{ py: 2, textAlign: 'center' }}>
            No dilution series available
          </Typography>
        )}
      </DialogContent>

      <DialogActions>
        <Button onClick={onClose} disabled={saving}>
          Cancel
        </Button>
        <Button
          variant="contained"
          onClick={handleSave}
          disabled={saving || !selectedId || !hasChanged}
          startIcon={saving ? <CircularProgress size={16} /> : <Check />}
        >
          {saving ? 'Saving...' : 'Apply Selection'}
        </Button>
      </DialogActions>
    </Dialog>
  );
}
