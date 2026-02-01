'use client';

import { useState } from 'react';
import {
  Box,
  TextField,
  Typography,
  Alert,
  InputAdornment,
  Grid,
  Tooltip,
  IconButton,
} from '@mui/material';
import { HelpOutline } from '@mui/icons-material';
import type { FittingParameters } from '@/types/models';

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
    </Box>
  );
}
