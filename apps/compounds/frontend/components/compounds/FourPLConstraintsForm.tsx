'use client';

import {
  Box,
  FormControlLabel,
  Checkbox,
  Typography,
  Tooltip,
  IconButton,
} from '@mui/material';
import { HelpOutline } from '@mui/icons-material';
import type { FittingParameters } from '@/types/compounds/models';

interface FourPLConstraintsFormProps {
  value: FittingParameters;
  onChange: (params: FittingParameters) => void;
  disabled?: boolean;
}

/**
 * Form for configuring 4PL curve fitting constraints.
 *
 * Options:
 * - Fix Hill coefficient to 1.0 (standard dose-response assumption)
 * - Fix top asymptote to control max value
 * - Fix bottom asymptote to control min value
 */
export function FourPLConstraintsForm({
  value,
  onChange,
  disabled = false,
}: FourPLConstraintsFormProps) {
  return (
    <Box>
      <Box sx={{ display: 'flex', alignItems: 'center', mb: 1 }}>
        <Typography variant="subtitle2" color="text.secondary">
          Curve Fitting Constraints
        </Typography>
        <Tooltip title="Optional constraints for 4PL curve fitting. These can improve fits for certain assay types.">
          <IconButton size="small" sx={{ ml: 0.5 }}>
            <HelpOutline fontSize="small" />
          </IconButton>
        </Tooltip>
      </Box>

      <Box sx={{ display: 'flex', flexDirection: 'column', gap: 0.5 }}>
        <FormControlLabel
          control={
            <Checkbox
              checked={value.fix_hill === 1.0}
              onChange={(e) =>
                onChange({
                  ...value,
                  fix_hill: e.target.checked ? 1.0 : null,
                })
              }
              disabled={disabled}
              size="small"
            />
          }
          label={
            <Typography variant="body2">
              Fix Hill coefficient to 1.0
            </Typography>
          }
        />

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
            <Typography variant="body2">
              Fix top asymptote to max control value
            </Typography>
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
            <Typography variant="body2">
              Fix bottom asymptote to min control value
            </Typography>
          }
        />
      </Box>
    </Box>
  );
}
