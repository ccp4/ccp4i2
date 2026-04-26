'use client';

import { Paper, LinearProgress, Typography } from '@mui/material';

type StatusKind = 'loading' | 'no-selection' | 'no-data';

const MESSAGE: Record<Exclude<StatusKind, 'loading'>, string> = {
  'no-selection': 'Select a target, protocol, or compound above to see results.',
  'no-data': 'No data found matching your criteria.',
};

/**
 * Shared placeholder Paper for the aggregation table's three pre-data
 * states: loading spinner, "select something", and "no rows matched".
 */
export function AggregationStatus({
  state,
  fillHeight = false,
}: {
  state: StatusKind;
  fillHeight?: boolean;
}) {
  const sx = { p: 3, ...(fillHeight && { height: '100%' }) };
  if (state === 'loading') {
    return (
      <Paper sx={sx}>
        <LinearProgress />
        <Typography sx={{ mt: 2 }} color="text.secondary" align="center">
          Running aggregation query...
        </Typography>
      </Paper>
    );
  }
  return (
    <Paper sx={sx}>
      <Typography color="text.secondary" align="center">
        {MESSAGE[state]}
      </Typography>
    </Paper>
  );
}
