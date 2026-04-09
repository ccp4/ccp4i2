'use client';

import {
  Box,
  Typography,
  FormControlLabel,
  Checkbox,
  Button,
  Tooltip,
  IconButton,
  Alert,
  CircularProgress,
} from '@mui/material';
import { HelpOutline, RestartAlt } from '@mui/icons-material';
import { useCompoundsApi } from '@/lib/compounds/api';

interface FlagInfo {
  id: string;
  label: string;
  description: string;
}

interface ValidationFlagsResponse {
  flags: FlagInfo[];
  default_invalidating: string[];
}

interface ValidationRulesFormProps {
  value: string[];  // Currently selected invalidating flags
  onChange: (flags: string[]) => void;
  disabled?: boolean;
}

export function ValidationRulesForm({
  value,
  onChange,
  disabled = false,
}: ValidationRulesFormProps) {
  const api = useCompoundsApi();

  // Fetch available flags from backend
  const { data, isLoading, error } = api.get<ValidationFlagsResponse>('validation-flags/');

  const handleToggle = (flagId: string) => {
    if (value.includes(flagId)) {
      onChange(value.filter((f) => f !== flagId));
    } else {
      onChange([...value, flagId]);
    }
  };

  const handleResetToDefaults = () => {
    if (data?.default_invalidating) {
      onChange(data.default_invalidating);
    }
  };

  const isUsingDefaults = data?.default_invalidating &&
    value.length === data.default_invalidating.length &&
    value.every((f) => data.default_invalidating.includes(f));

  if (isLoading) {
    return (
      <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, py: 2 }}>
        <CircularProgress size={20} />
        <Typography variant="body2" color="text.secondary">
          Loading validation options...
        </Typography>
      </Box>
    );
  }

  if (error || !data) {
    return (
      <Alert severity="error" sx={{ mt: 1 }}>
        Failed to load validation options
      </Alert>
    );
  }

  return (
    <Box>
      <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', mb: 1 }}>
        <Box sx={{ display: 'flex', alignItems: 'center' }}>
          <Typography variant="subtitle2" color="text.secondary">
            Validation Rules
          </Typography>
          <Tooltip title="Select which fitting quality flags should cause a measurement to be marked as invalid. Unchecked flags will appear as warnings but allow the measurement to remain valid.">
            <IconButton size="small" sx={{ ml: 0.5 }}>
              <HelpOutline fontSize="small" />
            </IconButton>
          </Tooltip>
        </Box>
        {!isUsingDefaults && (
          <Button
            size="small"
            startIcon={<RestartAlt />}
            onClick={handleResetToDefaults}
            disabled={disabled}
          >
            Reset to Defaults
          </Button>
        )}
      </Box>

      <Typography variant="caption" color="text.secondary" sx={{ display: 'block', mb: 1.5 }}>
        Checked flags will cause automatic invalidation:
      </Typography>

      <Box sx={{ display: 'flex', flexDirection: 'column', gap: 0.5 }}>
        {data.flags.map((flag) => (
          <Tooltip key={flag.id} title={flag.description} placement="right">
            <FormControlLabel
              control={
                <Checkbox
                  size="small"
                  checked={value.includes(flag.id)}
                  onChange={() => handleToggle(flag.id)}
                  disabled={disabled}
                />
              }
              label={
                <Typography variant="body2">
                  {flag.label}
                  <Typography
                    component="span"
                    variant="caption"
                    color="text.secondary"
                    sx={{ ml: 1 }}
                  >
                    ({flag.description})
                  </Typography>
                </Typography>
              }
              sx={{ ml: 0 }}
            />
          </Tooltip>
        ))}
      </Box>
    </Box>
  );
}
