'use client';

import { useState, useEffect, useMemo } from 'react';
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  TextField,
  Box,
  Alert,
  IconButton,
  CircularProgress,
  ToggleButtonGroup,
  ToggleButton,
  Autocomplete,
  Typography,
  Divider,
  FormControlLabel,
  Checkbox,
} from '@mui/material';
import { Close, Delete, Warning, SwapHoriz } from '@mui/icons-material';
import { useCompoundsApi, apiGet, apiPost } from '@/lib/compounds/api';
import type { Target } from '@/types/compounds/models';

interface ReassignableTarget {
  id: string;
  name: string;
}

interface DeletionPreview {
  target: { id: string; name: string };
  compound_count: number;
  template_count: number;
  hypothesis_count: number;
  assay_count: number;
  protocol_count: number;
  child_target_count: number;
  cascade_impact: {
    batches: number;
    batch_qc_files: number;
    compound_documents: number;
    assay_data_series: number;
  };
  reassignable_targets: ReassignableTarget[];
  affected_sample: {
    compounds: string[];
  };
}

interface ReassignResult {
  compounds_reassigned: number;
  target_deleted: boolean;
  replacement_target: { id: string; name: string };
}

interface CascadeResult {
  compounds_deleted: number;
  batches_deleted: number;
  files_removed: number;
  target_deleted: boolean;
}

type Mode = 'reassign' | 'cascade';

interface TargetDeletionDialogProps {
  open: boolean;
  target: Target | null;
  onClose: () => void;
  onDeleted: (summary: string) => void;
}

export function TargetDeletionDialog({
  open,
  target,
  onClose,
  onDeleted,
}: TargetDeletionDialogProps) {
  const api = useCompoundsApi();
  const [preview, setPreview] = useState<DeletionPreview | null>(null);
  const [previewLoading, setPreviewLoading] = useState(false);
  const [previewError, setPreviewError] = useState<string | null>(null);

  const [mode, setMode] = useState<Mode>('reassign');
  const [replacement, setReplacement] = useState<ReassignableTarget | null>(null);
  const [typedName, setTypedName] = useState('');
  const [typedCount, setTypedCount] = useState('');
  const [verifiedCheckbox, setVerifiedCheckbox] = useState(false);
  const [countMismatch, setCountMismatch] = useState(false);

  const [submitting, setSubmitting] = useState(false);
  const [submitError, setSubmitError] = useState<string | null>(null);

  const loadPreview = async (targetId: string) => {
    setPreviewLoading(true);
    setPreviewError(null);
    try {
      const data = await apiGet<DeletionPreview>(
        `targets/${targetId}/deletion_preview/`
      );
      setPreview(data);
    } catch (err) {
      setPreviewError(
        err instanceof Error ? err.message : 'Failed to load deletion preview'
      );
    } finally {
      setPreviewLoading(false);
    }
  };

  useEffect(() => {
    if (open && target) {
      setMode('reassign');
      setReplacement(null);
      setTypedName('');
      setTypedCount('');
      setVerifiedCheckbox(false);
      setCountMismatch(false);
      setSubmitError(null);
      setPreview(null);
      loadPreview(target.id);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [open, target?.id]);

  const expectedName = preview?.target.name ?? target?.name ?? '';
  const expectedCount = preview?.compound_count ?? 0;
  const nameMatches = typedName === expectedName;
  const countMatches = typedCount === String(expectedCount);

  const reassignReady = mode === 'reassign' && !!replacement && nameMatches;
  const cascadeReady =
    mode === 'cascade' && nameMatches && countMatches && verifiedCheckbox;
  const canSubmit = (reassignReady || cascadeReady) && !submitting && !!preview;

  const reassignOptions: ReassignableTarget[] = useMemo(
    () => preview?.reassignable_targets ?? [],
    [preview]
  );

  const handleSubmit = async () => {
    if (!target || !preview) return;
    setSubmitting(true);
    setSubmitError(null);
    setCountMismatch(false);

    try {
      if (mode === 'reassign') {
        if (!replacement) {
          setSubmitError('Select a replacement target');
          setSubmitting(false);
          return;
        }
        const result = await apiPost<ReassignResult>(
          `targets/${target.id}/delete_with_reassign/`,
          { replacement_target_id: replacement.id }
        );
        onDeleted(
          `Reassigned ${result.compounds_reassigned} compounds to ${result.replacement_target.name} and deleted ${expectedName}`
        );
        onClose();
      } else {
        try {
          const result = await apiPost<CascadeResult>(
            `targets/${target.id}/delete_with_cascade/`,
            { confirmed_compound_count: expectedCount }
          );
          onDeleted(
            `Deleted ${result.compounds_deleted} compounds, ${result.batches_deleted} batches, and target ${expectedName}`
          );
          onClose();
        } catch (err: unknown) {
          const status = (err as { status?: number })?.status;
          if (status === 409) {
            setCountMismatch(true);
            setSubmitError(
              'Compound count changed since preview loaded. Reloading…'
            );
            await loadPreview(target.id);
            setTypedCount('');
            return;
          }
          throw err;
        }
      }
    } catch (err) {
      setSubmitError(err instanceof Error ? err.message : 'Delete failed');
    } finally {
      setSubmitting(false);
    }
  };

  const renderImpactSummary = () => {
    if (!preview) return null;
    const lines: string[] = [];
    lines.push(
      `${preview.compound_count} compounds, ${preview.cascade_impact.batches} batches, ${preview.cascade_impact.batch_qc_files} QC files, ${preview.cascade_impact.compound_documents} compound documents`
    );
    if (preview.template_count > 0) {
      lines.push(
        `${preview.template_count} compound template${preview.template_count === 1 ? '' : 's'} will be deleted (belong to this target)`
      );
    }
    if (preview.assay_count > 0 || preview.protocol_count > 0) {
      const parts: string[] = [];
      if (preview.assay_count > 0)
        parts.push(`${preview.assay_count} assay${preview.assay_count === 1 ? '' : 's'}`);
      if (preview.protocol_count > 0)
        parts.push(
          `${preview.protocol_count} protocol${preview.protocol_count === 1 ? '' : 's'}`
        );
      lines.push(`${parts.join(' + ')} will have their target cleared`);
    }
    if (preview.hypothesis_count > 0) {
      lines.push(
        `${preview.hypothesis_count} hypothes${preview.hypothesis_count === 1 ? 'is' : 'es'} will be deleted`
      );
    }
    if (preview.child_target_count > 0) {
      lines.push(
        `${preview.child_target_count} child target${preview.child_target_count === 1 ? '' : 's'} will become root-level`
      );
    }
    if (preview.cascade_impact.assay_data_series > 0) {
      lines.push(
        `${preview.cascade_impact.assay_data_series} assay data series will be orphaned (compound reference cleared)`
      );
    }
    return (
      <Box sx={{ mb: 2 }}>
        {lines.map((line, i) => (
          <Typography key={i} variant="body2" sx={{ mb: 0.25 }}>
            • {line}
          </Typography>
        ))}
        {preview.affected_sample.compounds.length > 0 && (
          <Typography variant="caption" color="text.secondary" sx={{ display: 'block', mt: 1 }}>
            First compounds: {preview.affected_sample.compounds.join(', ')}
            {preview.compound_count > preview.affected_sample.compounds.length &&
              ` … and ${preview.compound_count - preview.affected_sample.compounds.length} more`}
          </Typography>
        )}
      </Box>
    );
  };

  return (
    <Dialog open={open} onClose={submitting ? undefined : onClose} maxWidth="sm" fullWidth>
      <DialogTitle sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Warning color="warning" />
          Delete target{expectedName ? ` "${expectedName}"` : ''}?
        </Box>
        <IconButton onClick={onClose} size="small" disabled={submitting}>
          <Close />
        </IconButton>
      </DialogTitle>

      <DialogContent dividers>
        {previewLoading && (
          <Box sx={{ display: 'flex', justifyContent: 'center', py: 4 }}>
            <CircularProgress />
          </Box>
        )}

        {previewError && (
          <Alert severity="error" sx={{ mb: 2 }}>
            {previewError}
          </Alert>
        )}

        {preview && (
          <>
            <ToggleButtonGroup
              value={mode}
              exclusive
              onChange={(_, next) => next && setMode(next as Mode)}
              fullWidth
              size="small"
              sx={{ mb: 2 }}
            >
              <ToggleButton value="reassign">
                <SwapHoriz fontSize="small" sx={{ mr: 0.5 }} />
                Reassign
              </ToggleButton>
              <ToggleButton value="cascade" color="error">
                <Delete fontSize="small" sx={{ mr: 0.5 }} />
                Delete everything
              </ToggleButton>
            </ToggleButtonGroup>

            {renderImpactSummary()}

            <Divider sx={{ my: 2 }} />

            {mode === 'reassign' && (
              <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
                <Autocomplete
                  options={reassignOptions}
                  value={replacement}
                  onChange={(_, v) => setReplacement(v)}
                  getOptionLabel={(opt) => opt.name}
                  isOptionEqualToValue={(a, b) => a.id === b.id}
                  renderInput={(params) => (
                    <TextField
                      {...params}
                      label="Replacement target"
                      required
                      size="small"
                      helperText="Compounds in the deleted target will be moved here"
                    />
                  )}
                />
                <TextField
                  label={`Type "${expectedName}" to confirm`}
                  value={typedName}
                  onChange={(e) => setTypedName(e.target.value)}
                  size="small"
                  fullWidth
                  error={typedName.length > 0 && !nameMatches}
                />
              </Box>
            )}

            {mode === 'cascade' && (
              <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
                <Alert severity="error" icon={<Warning />}>
                  This will permanently delete {preview.compound_count} compounds,{' '}
                  {preview.cascade_impact.batches} batches,{' '}
                  {preview.cascade_impact.batch_qc_files} QC files, and{' '}
                  {preview.cascade_impact.compound_documents} documents. This cannot be
                  undone.
                </Alert>
                <TextField
                  label={`Type "${expectedName}" to confirm`}
                  value={typedName}
                  onChange={(e) => setTypedName(e.target.value)}
                  size="small"
                  fullWidth
                  error={typedName.length > 0 && !nameMatches}
                />
                <TextField
                  label={`Type ${expectedCount} to confirm the number of compounds`}
                  value={typedCount}
                  onChange={(e) => setTypedCount(e.target.value)}
                  size="small"
                  fullWidth
                  error={(typedCount.length > 0 && !countMatches) || countMismatch}
                  helperText={
                    countMismatch
                      ? `Count changed to ${expectedCount}. Re-type to confirm.`
                      : undefined
                  }
                />
                <FormControlLabel
                  control={
                    <Checkbox
                      checked={verifiedCheckbox}
                      onChange={(e) => setVerifiedCheckbox(e.target.checked)}
                    />
                  }
                  label="I have verified this target's compounds cannot be reassigned"
                />
              </Box>
            )}

            {submitError && (
              <Alert severity="error" sx={{ mt: 2 }}>
                {submitError}
              </Alert>
            )}
          </>
        )}
      </DialogContent>

      <DialogActions>
        <Button onClick={onClose} disabled={submitting}>
          Cancel
        </Button>
        <Button
          variant="contained"
          color={mode === 'cascade' ? 'error' : 'primary'}
          onClick={handleSubmit}
          disabled={!canSubmit}
          startIcon={submitting ? <CircularProgress size={16} /> : <Delete />}
        >
          {submitting
            ? mode === 'cascade'
              ? 'Deleting…'
              : 'Reassigning…'
            : mode === 'cascade'
              ? 'Delete Everything'
              : 'Reassign & Delete'}
        </Button>
      </DialogActions>
    </Dialog>
  );
}
