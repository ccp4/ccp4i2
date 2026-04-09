'use client';

import { useState } from 'react';
import {
  Box,
  TextField,
  Typography,
  Alert,
  InputAdornment,
  Grid2 as Grid,
  Tooltip,
  IconButton,
  Divider,
  FormControlLabel,
  Checkbox,
} from '@mui/material';
import { HelpOutline } from '@mui/icons-material';
import type { FittingParameters } from '@/types/compounds/models';

interface TightBindingParametersFormProps {
  value: FittingParameters;
  onChange: (params: FittingParameters) => void;
  disabled?: boolean;
  unit?: string;  // Concentration unit from DilutionSeries (e.g., "uM", "nM")
}

export function TightBindingParametersForm({
  value,
  onChange,
  disabled = false,
  unit = 'nM',  // Default to nM if not specified
}: TightBindingParametersFormProps) {
  const [errors, setErrors] = useState<Record<string, string>>({});

  const validateAndUpdate = (field: keyof FittingParameters, strValue: string) => {
    const numValue = strValue === '' ? undefined : parseFloat(strValue);

    const newErrors = { ...errors };
    if (numValue !== undefined && (isNaN(numValue) || numValue <= 0)) {
      newErrors[field] = 'Must be a positive number';
    } else {
      delete newErrors[field];
    }
    setErrors(newErrors);

    onChange({
      ...value,
      [field]: numValue,
    });
  };

  const isComplete =
    value.protein_conc !== undefined &&
    value.ligand_conc !== undefined &&
    value.ligand_kd !== undefined;

  return (
    <Box>
      <Box sx={{ display: 'flex', alignItems: 'center', mb: 2 }}>
        <Typography variant="subtitle2" color="text.secondary">
          Tight-Binding Competition Parameters
        </Typography>
        <Tooltip title="Required for Wang equation analysis. These values are typically determined from your assay conditions.">
          <IconButton size="small" sx={{ ml: 0.5 }}>
            <HelpOutline fontSize="small" />
          </IconButton>
        </Tooltip>
      </Box>

      {!isComplete && (
        <Alert severity="info" sx={{ mb: 2 }}>
          All three parameters are required for tight-binding analysis.
        </Alert>
      )}

      <Grid container spacing={2}>
        <Grid size={{ xs: 12, md: 4 }}>
          <TextField
            label="Protein Concentration [P]t"
            value={value.protein_conc ?? ''}
            onChange={(e) => validateAndUpdate('protein_conc', e.target.value)}
            error={!!errors.protein_conc}
            helperText={errors.protein_conc || 'Total protein in assay'}
            slotProps={{
              input: {
                endAdornment: <InputAdornment position="end">{unit}</InputAdornment>,
              },
              htmlInput: { min: 0, step: 'any' },
            }}
            type="number"
            fullWidth
            disabled={disabled}
            size="small"
          />
        </Grid>
        <Grid size={{ xs: 12, md: 4 }}>
          <TextField
            label="Ligand Concentration [L]t"
            value={value.ligand_conc ?? ''}
            onChange={(e) => validateAndUpdate('ligand_conc', e.target.value)}
            error={!!errors.ligand_conc}
            helperText={errors.ligand_conc || 'Total labeled competitor'}
            slotProps={{
              input: {
                endAdornment: <InputAdornment position="end">{unit}</InputAdornment>,
              },
              htmlInput: { min: 0, step: 'any' },
            }}
            type="number"
            fullWidth
            disabled={disabled}
            size="small"
          />
        </Grid>
        <Grid size={{ xs: 12, md: 4 }}>
          <TextField
            label="Ligand Kd"
            value={value.ligand_kd ?? ''}
            onChange={(e) => validateAndUpdate('ligand_kd', e.target.value)}
            error={!!errors.ligand_kd}
            helperText={errors.ligand_kd || 'Kd of labeled competitor'}
            slotProps={{
              input: {
                endAdornment: <InputAdornment position="end">{unit}</InputAdornment>,
              },
              htmlInput: { min: 0, step: 'any' },
            }}
            type="number"
            fullWidth
            disabled={disabled}
            size="small"
          />
        </Grid>
      </Grid>

      {/* Asymptote Constraints Section */}
      <Divider sx={{ my: 2 }} />

      <Box sx={{ display: 'flex', alignItems: 'center', mb: 1 }}>
        <Typography variant="caption" color="text.secondary">
          Asymptote Constraints
        </Typography>
        <Tooltip title="Optional constraints for curve fitting. Hard constraints force exact values; soft restraints guide the fit while allowing deviation.">
          <IconButton size="small" sx={{ ml: 0.5 }}>
            <HelpOutline fontSize="small" />
          </IconButton>
        </Tooltip>
      </Box>

      <Box sx={{ display: 'flex', flexDirection: 'column', gap: 0.5 }}>
        {/* Hard constraints */}
        <FormControlLabel
          control={
            <Checkbox
              checked={value.fix_top === true}
              onChange={(e) =>
                onChange({
                  ...value,
                  fix_top: e.target.checked ? true : null,
                })
              }
              disabled={disabled}
              size="small"
            />
          }
          label={
            <Box sx={{ display: 'flex', alignItems: 'center' }}>
              <Typography variant="body2">
                Fix top to max control
              </Typography>
              <Tooltip title="Hard constraint: top asymptote will be exactly equal to the max control value">
                <HelpOutline sx={{ fontSize: 14, ml: 0.5, color: 'text.secondary' }} />
              </Tooltip>
            </Box>
          }
        />

        <FormControlLabel
          control={
            <Checkbox
              checked={value.fix_bottom === true}
              onChange={(e) =>
                onChange({
                  ...value,
                  fix_bottom: e.target.checked ? true : null,
                })
              }
              disabled={disabled}
              size="small"
            />
          }
          label={
            <Box sx={{ display: 'flex', alignItems: 'center' }}>
              <Typography variant="body2">
                Fix bottom to min control
              </Typography>
              <Tooltip title="Hard constraint: bottom asymptote will be exactly equal to the min control value">
                <HelpOutline sx={{ fontSize: 14, ml: 0.5, color: 'text.secondary' }} />
              </Tooltip>
            </Box>
          }
        />

        {/* Soft constraints */}
        <FormControlLabel
          control={
            <Checkbox
              checked={value.restrain_to_controls === true}
              onChange={(e) =>
                onChange({
                  ...value,
                  restrain_to_controls: e.target.checked ? true : null,
                })
              }
              disabled={disabled || (value.fix_top === true && value.fix_bottom === true)}
              size="small"
            />
          }
          label={
            <Box sx={{ display: 'flex', alignItems: 'center' }}>
              <Typography
                variant="body2"
                color={(value.fix_top === true && value.fix_bottom === true) ? 'text.disabled' : 'text.primary'}
              >
                Restrain to controls (soft)
              </Typography>
              <Tooltip title="Soft constraint: guides asymptotes toward control values using pseudo data points, while allowing deviation if data suggests different values. Only applies to asymptotes not already hard-fixed above.">
                <HelpOutline sx={{ fontSize: 14, ml: 0.5, color: 'text.secondary' }} />
              </Tooltip>
            </Box>
          }
        />
        {(value.fix_top === true && value.fix_bottom === true) && value.restrain_to_controls && (
          <Typography variant="caption" color="warning.main" sx={{ ml: 4 }}>
            Soft restraints have no effect when both asymptotes are hard-fixed
          </Typography>
        )}
      </Box>
    </Box>
  );
}
