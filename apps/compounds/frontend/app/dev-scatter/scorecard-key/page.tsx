'use client';

/**
 * Dev fixture for ScorecardKeyTable + the editor's ScorecardKeyPanel.
 *
 *   cd apps/compounds/frontend && npm run dev
 *   open http://localhost:3001/dev-scatter/scorecard-key
 *
 * Sample config exercises all four axis kinds (protocol, ratio, worst_of,
 * lipinski) with realistic-length names so we can spot truncation /
 * width / overlap issues without the Azure deploy round-trip.
 */

import { useState } from 'react';
import {
  Box,
  Container,
  IconButton,
  Paper,
  Slider,
  Tooltip,
  Typography,
} from '@mui/material';
import { Check, ContentCopy } from '@mui/icons-material';
import html2canvas from 'html2canvas';
import { useRef } from 'react';
import { ScorecardKeyTable } from '@/components/compounds/ScorecardKeyTable';
import type { ScorecardConfig } from '@/types/compounds/models';
import type { ProtocolInfo } from '@/types/compounds/aggregation';

// --- Fixture ----------------------------------------------------------------

const PROTOCOLS: ProtocolInfo[] = [
  { id: 'p-tm',     name: 'mEGFR TR-FRET TM',                        kpi_unit: 'nM' },
  { id: 'p-wt',     name: 'mEGFR TR-FRET WT',                        kpi_unit: 'nM' },
  { id: 'p-pc9',    name: 'mEGFR Cellular HTRF - PC9',               kpi_unit: 'nM' },
  { id: 'p-h1975',  name: 'mEGFR Cellular HTRF - H1975',             kpi_unit: 'nM' },
  { id: 'p-h1975g', name: 'mEGFR Cellular Cyquant - H1975',          kpi_unit: 'nM' },
  { id: 'p-a431',   name: 'mEGFR Cellular Cyquant - A431',           kpi_unit: 'nM' },
  { id: 'p-baf-d', name: 'mEGFR BaF3 CellTitreGlo del19_C797S GI50', kpi_unit: 'nM' },
  { id: 'p-baf-l', name: 'mEGFR BaF3 CellTitreGlo L858R_C797S GI50', kpi_unit: 'nM' },
];

const SAMPLE_CONFIG: ScorecardConfig = {
  axes: [
    {
      kind: 'protocol',
      label: 'FRET potency',
      sector: 'potency',
      protocol_id: 'p-tm',
      target_value: 10,
      poor_value: 10000,
      threshold_scale: 'log',
    },
    {
      kind: 'protocol',
      label: 'H1975 engage',
      sector: 'potency',
      protocol_id: 'p-h1975',
      target_value: 10,
      poor_value: 1000,
      threshold_scale: 'log',
    },
    {
      kind: 'protocol',
      label: 'H1975 growth',
      sector: 'potency',
      protocol_id: 'p-h1975g',
      target_value: 0.01,
      poor_value: 100000,
      threshold_scale: 'log',
    },
    {
      kind: 'ratio',
      label: 'FRET select',
      sector: 'selectivity',
      numerator_id: 'p-wt',
      denominator_id: 'p-tm',
      target_value: 100,
      poor_value: 1,
      threshold_scale: 'log',
    },
    {
      kind: 'ratio',
      label: 'H1975 grow select',
      sector: 'selectivity',
      numerator_id: 'p-a431',
      denominator_id: 'p-h1975g',
      target_value: 100,
      poor_value: 1,
      threshold_scale: 'log',
    },
    {
      kind: 'ratio',
      label: 'H1975 engage select',
      sector: 'selectivity',
      numerator_id: 'p-a431',
      denominator_id: 'p-h1975',
      target_value: 100,
      poor_value: 1,
      threshold_scale: 'log',
    },
    {
      kind: 'worst_of',
      label: 'BAF Selectivity',
      sector: 'selectivity',
      protocol_ids: ['p-baf-d', 'p-baf-l'],
      target_value: 100,
      poor_value: -2,
      threshold_scale: 'log',
    },
    {
      kind: 'lipinski',
      label: 'Lipinski compliance',
      sector: 'phys-props',
      target_value: 4,
      poor_value: 1,
      threshold_scale: 'linear',
    },
  ],
};

// --- Page -------------------------------------------------------------------

export default function ScorecardKeyDevPage() {
  const [width, setWidth] = useState(640);
  const keyRef = useRef<HTMLDivElement>(null);
  const [copied, setCopied] = useState(false);

  const handleCopy = () => {
    const node = keyRef.current;
    if (!node) return;
    const blobPromise = (async () => {
      if (typeof document !== 'undefined' && document.fonts?.ready) {
        await document.fonts.ready;
      }
      const canvas = await html2canvas(node, { backgroundColor: '#ffffff', scale: 2 });
      return new Promise<Blob>((resolve, reject) => {
        canvas.toBlob((b) => (b ? resolve(b) : reject(new Error('toBlob'))), 'image/png');
      });
    })();
    navigator.clipboard
      .write([new ClipboardItem({ 'image/png': blobPromise })])
      .then(() => {
        setCopied(true);
        setTimeout(() => setCopied(false), 1800);
      })
      .catch((err) => console.error('Copy failed:', err));
  };

  return (
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <Typography variant="h4" gutterBottom>
        Dev fixture: ScorecardKeyTable
      </Typography>
      <Typography variant="body2" color="text.secondary" sx={{ mb: 3 }}>
        Drag the slider to resize the key panel. Sample config exercises all
        four axis kinds with realistic protocol names so truncation /
        wrapping behaviour can be inspected without an Azure deploy.
      </Typography>

      <Paper sx={{ p: 2, mb: 3 }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
          <Typography variant="body2" sx={{ minWidth: 80 }}>
            Width: {width}px
          </Typography>
          <Slider
            min={320}
            max={1100}
            step={20}
            value={width}
            onChange={(_, v) => setWidth(typeof v === 'number' ? v : v[0])}
            sx={{ flex: 1, maxWidth: 400 }}
          />
        </Box>
      </Paper>

      <Paper sx={{ p: 2, width, transition: 'width 0.15s' }}>
        <Box sx={{ display: 'flex', alignItems: 'center', mb: 1 }}>
          <Typography variant="subtitle1" fontWeight={600} sx={{ flex: 1 }}>
            Scorecard key
          </Typography>
          <Tooltip title={copied ? 'Copied!' : 'Copy as image'} arrow>
            <IconButton onClick={handleCopy} size="small">
              {copied ? <Check fontSize="small" color="success" /> : <ContentCopy fontSize="small" />}
            </IconButton>
          </Tooltip>
        </Box>
        <Box
          ref={keyRef}
          sx={{
            bgcolor: '#fff',
            border: '1px solid',
            borderColor: 'divider',
            borderRadius: 1,
            p: 1.5,
          }}
        >
          <ScorecardKeyTable config={SAMPLE_CONFIG} protocols={PROTOCOLS} />
        </Box>
      </Paper>
    </Container>
  );
}
