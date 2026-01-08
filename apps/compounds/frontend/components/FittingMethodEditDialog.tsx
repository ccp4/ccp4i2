'use client';

import { useState, useEffect } from 'react';
import dynamic from 'next/dynamic';
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  TextField,
  Box,
  Typography,
  Alert,
  IconButton,
  Divider,
  CircularProgress,
  Chip,
} from '@mui/material';
import { Close, Code, Save } from '@mui/icons-material';
import { useCompoundsApi } from '@/lib/api';
import type { FittingMethodDetail } from '@/types/models';

// Dynamically import Monaco to avoid SSR issues
const Editor = dynamic(() => import('@monaco-editor/react'), { ssr: false });

interface FittingMethodEditDialogProps {
  open: boolean;
  onClose: () => void;
  fittingMethodId: string;
  onSave?: () => void;
}

export function FittingMethodEditDialog({
  open,
  onClose,
  fittingMethodId,
  onSave,
}: FittingMethodEditDialogProps) {
  const api = useCompoundsApi();
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Fetch the full fitting method details
  const { data: method, isLoading, mutate } = api.get<FittingMethodDetail>(
    open ? `fitting-methods/${fittingMethodId}/` : null
  );

  // Form state
  const [name, setName] = useState('');
  const [description, setDescription] = useState('');
  const [script, setScript] = useState('');
  const [version, setVersion] = useState('');

  // Reset form when method loads
  useEffect(() => {
    if (method) {
      setName(method.name);
      setDescription(method.description || '');
      setScript(method.script || '');
      setVersion(method.version);
      setError(null);
    }
  }, [method]);

  const handleSave = async () => {
    setSaving(true);
    setError(null);

    try {
      const response = await fetch(`/api/proxy/compounds/fitting-methods/${fittingMethodId}/`, {
        method: 'PATCH',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          name,
          description,
          script,
          version,
        }),
      });

      if (!response.ok) {
        const data = await response.json().catch(() => ({}));
        const errorMessages: string[] = [];
        for (const [field, errors] of Object.entries(data)) {
          if (Array.isArray(errors)) {
            errorMessages.push(`${field}: ${errors.join(', ')}`);
          } else if (typeof errors === 'string') {
            errorMessages.push(`${field}: ${errors}`);
          }
        }
        throw new Error(errorMessages.length > 0 ? errorMessages.join('; ') : 'Failed to save');
      }

      mutate();
      onSave?.();
      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to save');
    } finally {
      setSaving(false);
    }
  };

  const handlePasteFromClipboard = async () => {
    try {
      const text = await navigator.clipboard.readText();
      setScript(text);
    } catch (err) {
      setError('Failed to read clipboard. Please paste manually.');
    }
  };

  return (
    <Dialog open={open} onClose={onClose} maxWidth="lg" fullWidth>
      <DialogTitle sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Code color="primary" />
          Edit Fitting Method
          {method?.is_builtin && (
            <Chip label="Built-in" size="small" color="info" variant="outlined" />
          )}
        </Box>
        <IconButton onClick={onClose} size="small">
          <Close />
        </IconButton>
      </DialogTitle>

      <DialogContent dividers>
        {isLoading ? (
          <Box sx={{ display: 'flex', justifyContent: 'center', py: 4 }}>
            <CircularProgress />
          </Box>
        ) : method ? (
          <>
            {error && (
              <Alert severity="error" sx={{ mb: 2 }}>
                {error}
              </Alert>
            )}

            <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
              {/* Metadata */}
              <Box sx={{ display: 'flex', gap: 2 }}>
                <TextField
                  label="Name"
                  value={name}
                  onChange={(e) => setName(e.target.value)}
                  size="small"
                  sx={{ flex: 2 }}
                />
                <TextField
                  label="Version"
                  value={version}
                  onChange={(e) => setVersion(e.target.value)}
                  size="small"
                  sx={{ flex: 1 }}
                />
                <TextField
                  label="Slug"
                  value={method.slug}
                  size="small"
                  disabled
                  sx={{ flex: 1 }}
                />
              </Box>

              <TextField
                label="Description"
                value={description}
                onChange={(e) => setDescription(e.target.value)}
                multiline
                rows={2}
                size="small"
              />

              <Divider />

              {/* Script editor */}
              <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
                <Typography variant="subtitle2" color="text.secondary">
                  Python Script
                </Typography>
                <Button size="small" onClick={handlePasteFromClipboard}>
                  Paste from Clipboard
                </Button>
              </Box>

              <Box
                sx={{
                  border: 1,
                  borderColor: 'divider',
                  borderRadius: 1,
                  overflow: 'hidden',
                }}
              >
                <Editor
                  height="450px"
                  defaultLanguage="python"
                  value={script}
                  onChange={(value) => setScript(value || '')}
                  theme="vs-dark"
                  options={{
                    minimap: { enabled: false },
                    fontSize: 13,
                    lineNumbers: 'on',
                    scrollBeyondLastLine: false,
                    automaticLayout: true,
                    tabSize: 4,
                    wordWrap: 'on',
                  }}
                />
              </Box>

              <Typography variant="caption" color="text.secondary">
                The script must define a <code>fit(input_data: dict) -&gt; dict</code> function.
                Input contains: concentrations, responses, controls, parameters.
                Output should include: ic50/ki, top, bottom, r_squared, curve_points, flags, kpi.
              </Typography>
            </Box>
          </>
        ) : (
          <Typography color="error">Fitting method not found</Typography>
        )}
      </DialogContent>

      <DialogActions>
        <Button onClick={onClose} disabled={saving}>
          Cancel
        </Button>
        <Button
          variant="contained"
          onClick={handleSave}
          disabled={saving || !name.trim() || !script.trim()}
          startIcon={saving ? <CircularProgress size={16} /> : <Save />}
        >
          {saving ? 'Saving...' : 'Save Script'}
        </Button>
      </DialogActions>
    </Dialog>
  );
}
