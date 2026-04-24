'use client';

import { Fragment, useMemo } from 'react';
import { Box, Paper, Typography, Tooltip, Alert } from '@mui/material';
import { useRouter } from 'next/navigation';
import type {
  CompactAggregationResponse,
  CompactRow,
} from '@/types/compounds/aggregation';
import type { ScorecardConfig } from '@/types/compounds/models';
import { evaluateScorecard, type AxisEvaluation } from '@/lib/compounds/scorecard';
import { MoleculeChip } from '../MoleculeView';
import { routes } from '@/lib/compounds/routes';

interface Props {
  data: CompactAggregationResponse;
  scorecardConfig: ScorecardConfig | null | undefined;
  fillHeight?: boolean;
  searchTerm?: string;
}

/**
 * Bullet-chart small-multiples view: one row per compound, one cell per
 * scorecard axis. Each cell is a horizontal bar showing the compound's
 * normalised score against that axis, on a poor→excellent gradient.
 *
 * Designed for SAR-glance reading across many compounds at once — the
 * Cards-view spider gets tangled past ~6 compounds, but a bullets table
 * scales linearly. Same scoring logic, different layout.
 */
export function BulletsView({ data, scorecardConfig, fillHeight, searchTerm = '' }: Props) {
  const router = useRouter();

  const rows = useMemo(() => {
    const all = data.data as CompactRow[];
    if (!searchTerm) return all;
    const term = searchTerm.toLowerCase();
    return all.filter((r) => r.formatted_id?.toLowerCase().includes(term));
  }, [data.data, searchTerm]);

  if (!scorecardConfig?.axes?.length) {
    return (
      <Paper sx={{ p: 3 }}>
        <Alert severity="info">
          The bullets view needs a scorecard. Configure one on the target&apos;s
          page (Scorecard button) and pick this format again.
        </Alert>
      </Paper>
    );
  }

  const axes = scorecardConfig.axes;

  return (
    <Paper
      sx={{
        p: 2,
        overflow: 'auto',
        ...(fillHeight && { height: '100%' }),
      }}
    >
      <Box
        sx={{
          display: 'grid',
          gridTemplateColumns: `64px minmax(140px, max-content) repeat(${axes.length}, minmax(120px, 1fr))`,
          rowGap: 0.5,
          columnGap: 1.5,
          alignItems: 'center',
        }}
      >
        {/* Header row */}
        <Box />
        <Box />
        {axes.map((axis, i) => (
          <Tooltip key={i} title={axis.label} arrow>
            <Typography
              variant="caption"
              sx={{
                fontWeight: 600,
                color: 'text.secondary',
                whiteSpace: 'nowrap',
                overflow: 'hidden',
                textOverflow: 'ellipsis',
                pb: 0.5,
                borderBottom: 1,
                borderColor: 'divider',
              }}
            >
              {axis.label}
            </Typography>
          </Tooltip>
        ))}

        {/* Body rows */}
        {rows.map((row) => {
          const evals = evaluateScorecard(scorecardConfig, row);
          return (
            <Fragment key={row.compound_id}>
              <Box
                sx={{ cursor: 'pointer' }}
                onClick={() => router.push(routes.registry.compound(row.compound_id))}
              >
                {row.smiles ? (
                  <MoleculeChip smiles={row.smiles} size={56} />
                ) : (
                  <Box
                    sx={{
                      width: 56,
                      height: 56,
                      bgcolor: 'grey.100',
                      borderRadius: 0.5,
                      display: 'flex',
                      alignItems: 'center',
                      justifyContent: 'center',
                    }}
                  >
                    <Typography variant="caption" color="text.secondary">
                      —
                    </Typography>
                  </Box>
                )}
              </Box>
              <Typography
                variant="body2"
                sx={{
                  fontFamily: 'monospace',
                  fontWeight: 500,
                  cursor: 'pointer',
                  '&:hover': { textDecoration: 'underline' },
                }}
                onClick={() => router.push(routes.registry.compound(row.compound_id))}
              >
                {row.formatted_id}
              </Typography>
              {evals.map((ev, i) => (
                <BulletCell key={i} evaluation={ev} />
              ))}
            </Fragment>
          );
        })}
      </Box>
    </Paper>
  );
}

// ---------------------------------------------------------------------------

const GREEN_HUE = 140;
const RED_HUE = 5;

function hueAt(t: number): number {
  return RED_HUE + (GREEN_HUE - RED_HUE) * t;
}

/**
 * One bar in the small-multiples grid. Faint background gradient (poor→excellent)
 * with a coloured fill from 0 to t. Tooltip carries the raw value + tier so the
 * cell stays compact while remaining inspectable.
 */
function BulletCell({ evaluation }: { evaluation: AxisEvaluation }) {
  const { value, t } = evaluation;

  if (t == null) {
    return (
      <Tooltip title="No data" arrow>
        <Box
          sx={{
            height: 22,
            border: '1px dashed',
            borderColor: 'divider',
            borderRadius: 0.5,
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
          }}
        >
          <Typography variant="caption" color="text.disabled" sx={{ fontSize: '0.65rem' }}>
            —
          </Typography>
        </Box>
      </Tooltip>
    );
  }

  const fillColour = `hsl(${hueAt(t).toFixed(1)}, 65%, 62%)`;
  const bgGradient = `linear-gradient(to right, hsl(${RED_HUE}, 50%, 92%), hsl(${(RED_HUE + GREEN_HUE) / 2}, 50%, 92%), hsl(${GREEN_HUE}, 50%, 92%))`;

  const tooltip = `${formatValue(value)}  ${tierLabel(t)}`;

  return (
    <Tooltip title={tooltip} arrow>
      <Box
        sx={{
          position: 'relative',
          height: 22,
          background: bgGradient,
          border: 1,
          borderColor: 'divider',
          borderRadius: 0.5,
          overflow: 'hidden',
        }}
      >
        <Box
          sx={{
            position: 'absolute',
            left: 0,
            top: 0,
            bottom: 0,
            width: `${(t * 100).toFixed(1)}%`,
            backgroundColor: fillColour,
          }}
        />
        <Typography
          variant="caption"
          sx={{
            position: 'relative',
            display: 'block',
            textAlign: 'center',
            lineHeight: '22px',
            fontFamily: 'monospace',
            fontSize: '0.7rem',
            fontWeight: 500,
            color: t > 0.5 ? 'rgba(0,0,0,0.85)' : 'text.secondary',
          }}
        >
          {value != null ? formatValue(value) : ''}
        </Typography>
      </Box>
    </Tooltip>
  );
}

function formatValue(v: number | null): string {
  if (v == null || !Number.isFinite(v)) return '—';
  if (Math.abs(v) >= 100 || Math.abs(v) < 0.1) return v.toPrecision(3);
  return v.toFixed(2);
}

function tierLabel(t: number): string {
  if (t >= 1) return '(excellent)';
  if (t >= 2 / 3) return '(good)';
  if (t >= 1 / 3) return '(mid)';
  if (t > 0) return '(poor)';
  return '(failing)';
}
