'use client';

/**
 * TargetGenesSection — displays (and, when the user has edit rights,
 * inline-edits) the genes linked to a Target.
 *
 * Read mode: a row of chips, one per linked Gene. Tooltips expose full
 * HGNC metadata (name, aliases, UniProt/HGNC cross-ref links, hydration
 * status). Chip colour reflects hydration_status so partly-hydrated or
 * failed entries are visible at a glance.
 *
 * Edit mode: an Autocomplete with multiple+freeSolo on the same shape
 * the TargetCreateDialog uses. On save it PATCHes {gene_symbols: [...]}
 * to the target. Unknown symbols create pending Gene rows server-side;
 * the hydrate_genes management command fills them in asynchronously.
 */

import { useEffect, useState } from 'react';
import {
  Autocomplete,
  Box,
  Button,
  Chip,
  CircularProgress,
  IconButton,
  Link as MuiLink,
  Stack,
  TextField,
  Tooltip,
  Typography,
} from '@mui/material';
import { Cancel, Edit, Save, Biotech } from '@mui/icons-material';
import type { Gene, HydrationStatus, Target } from '@/types/compounds/models';
import { apiPatch } from '@/lib/compounds/api';

interface TargetGenesSectionProps {
  target: Target;
  canEdit: boolean;
  onUpdated: (updated: Target) => void;
}

const STATUS_COLOR: Record<HydrationStatus, 'default' | 'success' | 'warning' | 'error' | 'info'> = {
  ok: 'success',
  manual: 'info',
  partial: 'warning',
  pending: 'default',
  failed: 'error',
};

const STATUS_LABEL: Record<HydrationStatus, string> = {
  ok: 'HGNC-hydrated',
  manual: 'Manually curated',
  partial: 'Partial hydration',
  pending: 'Awaiting HGNC hydration',
  failed: 'HGNC lookup failed (non-human or unknown symbol)',
};

function GeneChipTooltip({ gene }: { gene: Gene }) {
  return (
    <Box sx={{ p: 0.5, maxWidth: 320 }}>
      <Typography variant="subtitle2" sx={{ fontWeight: 600 }}>
        {gene.symbol}
        {gene.name && (
          <Typography component="span" variant="body2" sx={{ ml: 1, color: 'grey.300' }}>
            — {gene.name}
          </Typography>
        )}
      </Typography>

      {gene.aliases.length > 0 && (
        <Typography variant="caption" sx={{ display: 'block', mt: 0.5 }}>
          <strong>Aliases:</strong> {gene.aliases.join(', ')}
        </Typography>
      )}

      <Typography variant="caption" sx={{ display: 'block', mt: 0.5 }}>
        <strong>Hydration:</strong> {STATUS_LABEL[gene.hydration_status]}
      </Typography>

      {(gene.hgnc_id || gene.uniprot_ids.length > 0 || gene.ensembl_gene_id) && (
        <Stack direction="row" spacing={1} sx={{ mt: 1 }} flexWrap="wrap">
          {gene.hgnc_id && (
            <MuiLink
              href={`https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/${gene.hgnc_id}`}
              target="_blank"
              rel="noopener"
              sx={{ fontSize: '0.75rem' }}
            >
              HGNC
            </MuiLink>
          )}
          {gene.uniprot_ids.slice(0, 2).map((u) => (
            <MuiLink
              key={u}
              href={`https://www.uniprot.org/uniprotkb/${u}`}
              target="_blank"
              rel="noopener"
              sx={{ fontSize: '0.75rem' }}
            >
              UniProt {u}
            </MuiLink>
          ))}
          {gene.ensembl_gene_id && (
            <MuiLink
              href={`https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=${gene.ensembl_gene_id}`}
              target="_blank"
              rel="noopener"
              sx={{ fontSize: '0.75rem' }}
            >
              Ensembl
            </MuiLink>
          )}
        </Stack>
      )}
    </Box>
  );
}

export function TargetGenesSection({ target, canEdit, onUpdated }: TargetGenesSectionProps) {
  const [editing, setEditing] = useState(false);
  const [symbols, setSymbols] = useState<string[]>([]);
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Reset local edit state whenever the target changes or edit mode toggles.
  useEffect(() => {
    setSymbols((target.genes ?? []).map((g) => g.symbol));
    setError(null);
  }, [target.id, editing, target.genes]);

  const genes = target.genes ?? [];

  const handleStartEdit = () => setEditing(true);
  const handleCancel = () => {
    setEditing(false);
    setError(null);
  };

  const handleSave = async () => {
    setSaving(true);
    setError(null);
    try {
      const cleaned = symbols
        .map((s) => s.trim().toUpperCase())
        .filter((s) => s.length > 0);
      const updated = await apiPatch<Target>(`targets/${target.id}/`, {
        gene_symbols: cleaned,
      });
      onUpdated(updated);
      setEditing(false);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to update genes');
    } finally {
      setSaving(false);
    }
  };

  return (
    <Box sx={{ mt: 2 }}>
      <Stack direction="row" alignItems="center" spacing={1} sx={{ mb: 1 }}>
        <Biotech fontSize="small" color="action" />
        <Typography variant="subtitle2" sx={{ fontWeight: 600 }}>
          Genes
        </Typography>
        {!editing && canEdit && (
          <Tooltip title="Edit gene symbols">
            <IconButton size="small" onClick={handleStartEdit}>
              <Edit fontSize="inherit" />
            </IconButton>
          </Tooltip>
        )}
      </Stack>

      {editing ? (
        <Stack spacing={1}>
          <Autocomplete
            multiple
            freeSolo
            size="small"
            options={[]}
            value={symbols}
            onChange={(_event, newValue) => {
              const expanded: string[] = [];
              for (const item of newValue) {
                if (typeof item !== 'string') continue;
                for (const piece of item.split(/[\s,]+/)) {
                  const trimmed = piece.trim().toUpperCase();
                  if (trimmed && !expanded.includes(trimmed)) expanded.push(trimmed);
                }
              }
              setSymbols(expanded);
            }}
            renderTags={(value, getTagProps) =>
              value.map((option, index) => {
                const tagProps = getTagProps({ index });
                return (
                  <Chip
                    variant="outlined"
                    label={option}
                    size="small"
                    {...tagProps}
                    key={option}
                  />
                );
              })
            }
            renderInput={(params) => (
              <TextField
                {...params}
                placeholder="e.g. EGFR, ERBB2 — press Enter or comma after each"
                helperText="Unknown symbols create pending Gene rows. Submit [] to clear all genes."
              />
            )}
          />
          {error && (
            <Typography variant="caption" color="error">
              {error}
            </Typography>
          )}
          <Stack direction="row" spacing={1}>
            <Button
              variant="contained"
              size="small"
              startIcon={saving ? <CircularProgress size={14} /> : <Save />}
              onClick={handleSave}
              disabled={saving}
            >
              {saving ? 'Saving...' : 'Save'}
            </Button>
            <Button
              variant="outlined"
              size="small"
              startIcon={<Cancel />}
              onClick={handleCancel}
              disabled={saving}
            >
              Cancel
            </Button>
          </Stack>
        </Stack>
      ) : genes.length === 0 ? (
        <Typography variant="body2" color="text.secondary" sx={{ fontStyle: 'italic' }}>
          No genes linked.
          {canEdit && ' Use the edit button above to add HGNC symbols (e.g. EGFR, SKP2).'}
        </Typography>
      ) : (
        <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: 0.75 }}>
          {genes.map((gene) => (
            <Tooltip
              key={gene.id}
              title={<GeneChipTooltip gene={gene} />}
              arrow
              placement="top"
            >
              <Chip
                label={gene.symbol}
                size="small"
                color={STATUS_COLOR[gene.hydration_status]}
                variant={gene.hydration_status === 'ok' ? 'filled' : 'outlined'}
              />
            </Tooltip>
          ))}
        </Box>
      )}
    </Box>
  );
}
