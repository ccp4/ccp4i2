'use client';

import {
  Box,
  FormControlLabel,
  Checkbox,
  Typography,
  Tooltip,
  IconButton,
  Divider,
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
 * Hard constraints:
 * - Fix Hill coefficient to 1.0 (standard dose-response assumption)
 * - Fix top asymptote to control max value
 * - Fix bottom asymptote to control min value
 *
 * Soft constraints:
 * - Restrain asymptotes to control values via pseudo data points
 */
export function FourPLConstraintsForm({
  value,
  onChange,
  disabled = false,
}: FourPLConstraintsFormProps) {
  // Determine if soft restraints would have any effect
  // (only useful if at least one asymptote is not hard-fixed)
  const bothAsymptotesFixed = value.fix_top === true && value.fix_bottom === true;

  return (
    <Box>
      <Box sx={{ display: 'flex', alignItems: 'center', mb: 1 }}>
        <Typography variant="subtitle2" color="text.secondary">
          Curve Fitting Constraints
        </Typography>
        <Tooltip title="Optional constraints for 4PL curve fitting. Hard constraints force exact values; soft restraints guide the fit while allowing deviation.">
          <IconButton size="small" sx={{ ml: 0.5 }}>
            <HelpOutline fontSize="small" />
          </IconButton>
        </Tooltip>
      </Box>

      {/* Hill coefficient */}
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

        <Divider sx={{ my: 1 }} />

        {/* Asymptote constraints section */}
        <Typography variant="caption" color="text.secondary" sx={{ mb: 0.5 }}>
          Asymptote Constraints
        </Typography>

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
              disabled={disabled || bothAsymptotesFixed}
              size="small"
            />
          }
          label={
            <Box sx={{ display: 'flex', alignItems: 'center' }}>
              <Typography
                variant="body2"
                color={bothAsymptotesFixed ? 'text.disabled' : 'text.primary'}
              >
                Restrain to controls (soft)
              </Typography>
              <Tooltip title="Soft constraint: guides asymptotes toward control values using pseudo data points, while allowing deviation if data suggests different values. Only applies to asymptotes not already hard-fixed above.">
                <HelpOutline sx={{ fontSize: 14, ml: 0.5, color: 'text.secondary' }} />
              </Tooltip>
            </Box>
          }
        />
        {bothAsymptotesFixed && value.restrain_to_controls && (
          <Typography variant="caption" color="warning.main" sx={{ ml: 4 }}>
            Soft restraints have no effect when both asymptotes are hard-fixed
          </Typography>
        )}
      </Box>
    </Box>
  );
}
