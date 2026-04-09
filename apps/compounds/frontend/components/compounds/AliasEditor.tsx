'use client';

import { useState, useCallback } from 'react';
import {
  Box,
  Chip,
  TextField,
  IconButton,
  Typography,
  Tooltip,
  CircularProgress,
  Alert,
} from '@mui/material';
import { Add, Edit, Save, Close } from '@mui/icons-material';
import { useCompoundsApi } from '@/lib/compounds/api';

interface AliasEditorProps {
  compoundId: string;
  aliases: string[];
  canEdit: boolean;
  onUpdate?: (newAliases: string[]) => void;
}

/**
 * Component for viewing and editing compound aliases.
 * Aliases are alternative identifiers used during data import matching
 * (e.g., supplier codes, abbreviations, internal project names).
 */
export function AliasEditor({ compoundId, aliases, canEdit, onUpdate }: AliasEditorProps) {
  const api = useCompoundsApi();
  const [isEditing, setIsEditing] = useState(false);
  const [editedAliases, setEditedAliases] = useState<string[]>(aliases || []);
  const [newAlias, setNewAlias] = useState('');
  const [isSaving, setIsSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleStartEdit = useCallback(() => {
    setEditedAliases(aliases || []);
    setNewAlias('');
    setError(null);
    setIsEditing(true);
  }, [aliases]);

  const handleCancel = useCallback(() => {
    setEditedAliases(aliases || []);
    setNewAlias('');
    setError(null);
    setIsEditing(false);
  }, [aliases]);

  const handleAddAlias = useCallback(() => {
    const trimmed = newAlias.trim();
    if (!trimmed) return;

    // Check for duplicates (case-insensitive)
    const lowerTrimmed = trimmed.toLowerCase();
    if (editedAliases.some(a => a.toLowerCase() === lowerTrimmed)) {
      setError(`Alias "${trimmed}" already exists`);
      return;
    }

    setEditedAliases([...editedAliases, trimmed]);
    setNewAlias('');
    setError(null);
  }, [newAlias, editedAliases]);

  const handleRemoveAlias = useCallback((aliasToRemove: string) => {
    setEditedAliases(editedAliases.filter(a => a !== aliasToRemove));
  }, [editedAliases]);

  const handleKeyDown = useCallback((e: React.KeyboardEvent) => {
    if (e.key === 'Enter') {
      e.preventDefault();
      handleAddAlias();
    }
  }, [handleAddAlias]);

  const handleSave = useCallback(async () => {
    setIsSaving(true);
    setError(null);

    try {
      await api.patch(`compounds/${compoundId}/`, {
        aliases: editedAliases,
      });
      setIsEditing(false);
      onUpdate?.(editedAliases);
    } catch (err: any) {
      setError(err.message || 'Failed to save aliases');
    } finally {
      setIsSaving(false);
    }
  }, [api, compoundId, editedAliases, onUpdate]);

  // Read-only view
  if (!isEditing) {
    return (
      <Box>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, mb: 1 }}>
          <Typography variant="h6">
            Aliases
          </Typography>
          {canEdit && (
            <Tooltip title="Edit aliases">
              <IconButton size="small" onClick={handleStartEdit}>
                <Edit fontSize="small" />
              </IconButton>
            </Tooltip>
          )}
        </Box>

        {aliases && aliases.length > 0 ? (
          <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: 1 }}>
            {aliases.map((alias, idx) => (
              <Chip
                key={idx}
                label={alias}
                size="small"
                variant="outlined"
                sx={{ fontFamily: 'monospace' }}
              />
            ))}
          </Box>
        ) : (
          <Typography color="text.secondary" variant="body2">
            No aliases defined. Aliases help match this compound during data imports
            when alternative names are used (e.g., supplier codes, abbreviations).
          </Typography>
        )}
      </Box>
    );
  }

  // Edit mode
  return (
    <Box>
      <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', mb: 1 }}>
        <Typography variant="h6">
          Aliases
        </Typography>
        <Box sx={{ display: 'flex', gap: 1 }}>
          <Tooltip title="Save changes">
            <span>
              <IconButton
                size="small"
                onClick={handleSave}
                disabled={isSaving}
                color="primary"
              >
                {isSaving ? <CircularProgress size={20} /> : <Save fontSize="small" />}
              </IconButton>
            </span>
          </Tooltip>
          <Tooltip title="Cancel">
            <IconButton size="small" onClick={handleCancel} disabled={isSaving}>
              <Close fontSize="small" />
            </IconButton>
          </Tooltip>
        </Box>
      </Box>

      {error && (
        <Alert severity="error" sx={{ mb: 1 }} onClose={() => setError(null)}>
          {error}
        </Alert>
      )}

      {/* Current aliases */}
      <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: 1, mb: 2 }}>
        {editedAliases.map((alias, idx) => (
          <Chip
            key={idx}
            label={alias}
            size="small"
            variant="outlined"
            onDelete={() => handleRemoveAlias(alias)}
            sx={{ fontFamily: 'monospace' }}
          />
        ))}
        {editedAliases.length === 0 && (
          <Typography color="text.secondary" variant="body2" sx={{ fontStyle: 'italic' }}>
            No aliases yet
          </Typography>
        )}
      </Box>

      {/* Add new alias */}
      <Box sx={{ display: 'flex', gap: 1, alignItems: 'flex-start' }}>
        <TextField
          size="small"
          placeholder="Add alias (e.g., 'AST-001', 'compound-a')"
          value={newAlias}
          onChange={(e) => setNewAlias(e.target.value)}
          onKeyDown={handleKeyDown}
          sx={{ flex: 1 }}
          slotProps={{
            input: {
              sx: { fontFamily: 'monospace' }
            }
          }}
          helperText="Press Enter or click + to add"
        />
        <IconButton
          onClick={handleAddAlias}
          disabled={!newAlias.trim()}
          color="primary"
          sx={{ mt: 0.5 }}
        >
          <Add />
        </IconButton>
      </Box>
    </Box>
  );
}
