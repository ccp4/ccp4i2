'use client';

import {
  Box,
  IconButton,
  InputAdornment,
  TextField,
} from '@mui/material';
import { Search, Clear } from '@mui/icons-material';
import {
  AggregationResponse,
  ConcentrationDisplayMode,
  isCompactResponse,
} from '@/types/compounds/aggregation';
import { ConcentrationDisplaySelector, ThresholdLegend } from './shared';

interface Props {
  effectiveData: AggregationResponse | null | undefined;
  searchTerm: string;
  onSearchTermChange: (next: string) => void;
  displayMode: ConcentrationDisplayMode;
  onDisplayChange: (next: ConcentrationDisplayMode) => void;
}

/** Top-of-table controls: threshold legend (compact only), search,
 *  concentration display selector. */
export function AggregationToolbar({
  effectiveData,
  searchTerm,
  onSearchTermChange,
  displayMode,
  onDisplayChange,
}: Props) {
  return (
    <Box sx={{ display: 'flex', justifyContent: 'flex-end', alignItems: 'center', mb: 1, flexShrink: 0, gap: 1 }}>
      {effectiveData && isCompactResponse(effectiveData) && (
        <ThresholdLegend protocols={effectiveData.protocols} />
      )}
      <TextField
        size="small"
        placeholder="Search compounds..."
        value={searchTerm}
        onChange={(e) => onSearchTermChange(e.target.value)}
        slotProps={{
          input: {
            startAdornment: (
              <InputAdornment position="start">
                <Search fontSize="small" color="action" />
              </InputAdornment>
            ),
            endAdornment: searchTerm ? (
              <InputAdornment position="end">
                <IconButton size="small" onClick={() => onSearchTermChange('')} edge="end">
                  <Clear fontSize="small" />
                </IconButton>
              </InputAdornment>
            ) : null,
          },
        }}
        sx={{ width: 220 }}
      />
      <ConcentrationDisplaySelector value={displayMode} onChange={onDisplayChange} />
    </Box>
  );
}
