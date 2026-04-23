'use client';

/**
 * Dev-only fixture page for debugging the scatter plot / aggregation views
 * without needing a backend or a deployed container.
 *
 * Usage:
 *   cd apps/compounds/frontend && npm run dev
 *   open http://localhost:3001/__dev/scatter
 *
 * Swap the output format via the route's `?format=` param:
 *   ?format=cards    (default)
 *   ?format=pivot
 *   ?format=compact  (omit the param)
 *
 * To reproduce the current scatter crash, load the page, click the
 * bubble-chart icon on any protocol column, then pick axes.
 */

import { useSearchParams } from 'next/navigation';
import { useMemo } from 'react';
import { Box } from '@mui/material';
import { AggregationTable } from '@/components/compounds/AggregationTable';
import type {
  CompactAggregationResponse,
  AggregationType,
  OutputFormat,
  ProtocolInfo,
  CompactRow,
} from '@/types/compounds/aggregation';
import type { ScorecardConfig } from '@/types/compounds/models';

// --- Fixture ----------------------------------------------------------------

const PROTOCOLS: ProtocolInfo[] = [
  {
    id: 'proto-tm',
    name: 'mEGFR TR-FRET TM',
    kpi_unit: 'nM',
    target_value: 10,
    poor_value: 10000,
    threshold_scale: 'log',
  },
  {
    id: 'proto-wt',
    name: 'mEGFR TR-FRET WT',
    kpi_unit: 'nM',
    target_value: 10,
    poor_value: 10000,
    threshold_scale: 'log',
  },
  {
    id: 'proto-pc9',
    name: 'mEGFR Cellular HTRF - PC9',
    kpi_unit: 'nM',
    target_value: 500,
    poor_value: 30000,
    threshold_scale: 'log',
  },
  {
    id: 'proto-hlm',
    name: 'HLM (Astex) 1 uM',
    kpi_unit: 'µL/min/mg',
    target_value: 10,
    poor_value: 50,
    threshold_scale: 'log',
  },
  // Unconfigured protocol — renders uncoloured
  {
    id: 'proto-logd',
    name: 'HPLC LogD',
    kpi_unit: null,
  },
];

const COMPOUNDS: Array<{
  id: string;
  formatted_id: string;
  smiles: string;
  values: Record<string, number | null>;
}> = [
  { id: 'c1', formatted_id: 'NCL-00031070', smiles: 'CC(=O)N',        values: { 'proto-tm': 11,   'proto-wt': 320,   'proto-pc9': 100,  'proto-hlm': 12, 'proto-logd': 2.4 } },
  { id: 'c2', formatted_id: 'NCL-00031071', smiles: 'CCOC',           values: { 'proto-tm': 38,   'proto-wt': 970,   'proto-pc9': 300,  'proto-hlm': 18, 'proto-logd': 1.8 } },
  { id: 'c3', formatted_id: 'NCL-00031072', smiles: 'CCN',            values: { 'proto-tm': 7.1,  'proto-wt': 46,    'proto-pc9': 100,  'proto-hlm': 8,  'proto-logd': 2.0 } },
  { id: 'c4', formatted_id: 'NCL-00031136', smiles: 'c1ccccc1',       values: { 'proto-tm': 277,  'proto-wt': 6711,  'proto-pc9': 370,  'proto-hlm': 22, 'proto-logd': 3.1 } },
  { id: 'c5', formatted_id: 'NCL-00031032', smiles: 'Cc1ccccc1',      values: { 'proto-tm': 1.5,  'proto-wt': 1.8,   'proto-pc9': 9.6,  'proto-hlm': 6,  'proto-logd': 2.5 } },
  { id: 'c6', formatted_id: 'NCL-00031200', smiles: 'CCC(=O)O',       values: { 'proto-tm': 55,   'proto-wt': 420,   'proto-pc9': 1500, 'proto-hlm': 35, 'proto-logd': 2.2 } },
  { id: 'c7', formatted_id: 'NCL-00031201', smiles: 'CCCN',           values: { 'proto-tm': 120,  'proto-wt': 2800,  'proto-pc9': 4500, 'proto-hlm': 45, 'proto-logd': 1.5 } },
  { id: 'c8', formatted_id: 'NCL-00031202', smiles: 'CCCCO',          values: { 'proto-tm': 3.2,  'proto-wt': 68,    'proto-pc9': 45,   'proto-hlm': 10, 'proto-logd': 2.8 } },
  { id: 'c9', formatted_id: 'NCL-00031203', smiles: 'CNC',            values: { 'proto-tm': 800,  'proto-wt': null,  'proto-pc9': null, 'proto-hlm': null, 'proto-logd': null } },
];

function buildCompactResponse(): CompactAggregationResponse {
  // Molecular properties so the Lipinski axis evaluates against the fixture.
  // Mixed to exercise all four tiers of the Lipinski count.
  const lipinskiSeed = [
    { molecular_weight: 420, clogp: 3.8, hbd: 2, hba: 7 },   // passes all 4
    { molecular_weight: 480, clogp: 4.2, hbd: 3, hba: 8 },   // passes all 4
    { molecular_weight: 360, clogp: 2.9, hbd: 1, hba: 5 },   // passes all 4
    { molecular_weight: 560, clogp: 5.3, hbd: 4, hba: 11 },  // fails 4 (high MW, LogP, HBA)
    { molecular_weight: 320, clogp: 1.5, hbd: 1, hba: 4 },   // passes all 4
    { molecular_weight: 610, clogp: 6.1, hbd: 3, hba: 9 },   // fails 2
    { molecular_weight: 700, clogp: 7.0, hbd: 6, hba: 13 },  // fails all 4
    { molecular_weight: 380, clogp: 3.0, hbd: 2, hba: 6 },   // passes all 4
    { molecular_weight: 450, clogp: 4.5, hbd: 2, hba: 7 },   // passes all 4
  ];

  const data: CompactRow[] = COMPOUNDS.map((c, i) => ({
    compound_id: c.id,
    formatted_id: c.formatted_id,
    smiles: c.smiles,
    target_name: 'EGFR',
    properties: lipinskiSeed[i % lipinskiSeed.length],
    protocols: Object.fromEntries(
      PROTOCOLS.map((p) => {
        const v = c.values[p.id];
        return [
          p.id,
          v == null
            ? { geomean: null, count: 0, tested: 0, no_analysis: 0, invalid: 0, unassigned: 0 }
            : { geomean: v, count: 1, stdev: null, stdev_log: null, tested: 1, no_analysis: 0, invalid: 0, unassigned: 0 },
        ];
      }),
    ),
  }));

  const totalMeasurements = data.reduce(
    (acc, row) => acc + Object.values(row.protocols).reduce((s, p) => s + (p.count ?? 0), 0),
    0,
  );

  return {
    meta: {
      compound_count: data.length,
      row_count: data.length,
      protocol_count: PROTOCOLS.length,
      total_measurements: totalMeasurements,
    },
    protocols: PROTOCOLS,
    data,
  };
}

// --- Page -------------------------------------------------------------------

export default function ScatterDebugPage() {
  const searchParams = useSearchParams();
  const formatParam = searchParams?.get('format');
  const outputFormat: OutputFormat =
    formatParam === 'pivot' ? 'pivot' :
    formatParam === 'cards' ? 'cards' :
    formatParam === 'compact' ? 'compact' :
    'cards';
  const aggregations: AggregationType[] = useMemo(() => ['geomean', 'count'], []);
  const data = useMemo(buildCompactResponse, []);

  // Fixture scorecard exercising all four axis kinds — used to reproduce
  // Cards-view crashes (React error #310 etc.) locally with unminified traces.
  const scorecardConfig: ScorecardConfig = useMemo(
    () => ({
      axes: [
        {
          kind: 'protocol',
          label: 'TM potency',
          protocol_id: 'proto-tm',
          target_value: 10,
          poor_value: 10000,
          threshold_scale: 'log',
        },
        {
          kind: 'ratio',
          label: 'WT/TM selectivity',
          numerator_id: 'proto-wt',
          denominator_id: 'proto-tm',
          target_value: 100,
          poor_value: 1,
          threshold_scale: 'log',
        },
        {
          kind: 'worst_of',
          label: 'cellular worst',
          protocol_ids: ['proto-pc9'],
          target_value: 500,
          poor_value: 30000,
          threshold_scale: 'log',
        },
        {
          kind: 'lipinski',
          label: 'Lipinski',
          target_value: 4,
          poor_value: 1,
          threshold_scale: 'linear',
        },
      ],
    }),
    [],
  );

  return (
    <Box sx={{ p: 2, height: '100vh', display: 'flex', flexDirection: 'column' }}>
      <Box sx={{ mb: 1, fontSize: 14, color: 'text.secondary', fontFamily: 'monospace' }}>
        Dev fixture — output=<b>{outputFormat}</b> — swap with ?format=cards|pivot|compact
      </Box>
      <Box sx={{ flex: 1, minHeight: 0 }}>
        <AggregationTable
          data={data}
          aggregations={aggregations}
          outputFormat={outputFormat}
          scorecardConfig={scorecardConfig}
          fillHeight
        />
      </Box>
    </Box>
  );
}
