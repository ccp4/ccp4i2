'use client';

import { use, useState } from 'react';
import { useRouter } from 'next/navigation';
import {
  Container,
  Typography,
  Box,
  Paper,
  Grid2 as Grid,
  Chip,
  Skeleton,
  Divider,
  Button,
  TextField,
  Alert,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  List,
  ListItem,
  ListItemText,
  ListItemButton,
} from '@mui/material';
import { GridOn, Edit, Save, Cancel, Delete, ArrowBack } from '@mui/icons-material';
import Link from 'next/link';
import { useSWRConfig } from 'swr';
import { PageHeader } from '@/components/compounds/PageHeader';
import { PlateLayoutEditor } from '@/components/compounds/PlateLayoutEditor';
import { PlatePreview } from '@/components/compounds/PlatePreview';
import { useCompoundsApi } from '@/lib/compounds/api';
import { routes } from '@/lib/compounds/routes';
import type { PlateLayoutRecord, PlateLayout, Protocol } from '@/types/compounds/models';

interface PageProps {
  params: Promise<{ id: string }>;
}

function InfoRow({ label, value }: { label: string; value: React.ReactNode }) {
  return (
    <Box sx={{ display: 'flex', py: 0.5 }}>
      <Typography
        color="text.secondary"
        sx={{ minWidth: 140, fontWeight: 500 }}
      >
        {label}:
      </Typography>
      <Typography component="div">{value ?? '-'}</Typography>
    </Box>
  );
}

export default function PlateLayoutDetailPage({ params }: PageProps) {
  const { id } = use(params);
  const router = useRouter();
  const { mutate } = useSWRConfig();
  const api = useCompoundsApi();

  // Edit mode state
  const [isEditing, setIsEditing] = useState(false);
  const [editedName, setEditedName] = useState('');
  const [editedDescription, setEditedDescription] = useState('');
  const [editedConfig, setEditedConfig] = useState<PlateLayout | null>(null);
  const [saving, setSaving] = useState(false);
  const [saveError, setSaveError] = useState<string | null>(null);

  // Delete dialog state
  const [deleteDialogOpen, setDeleteDialogOpen] = useState(false);
  const [deleting, setDeleting] = useState(false);
  const [deleteError, setDeleteError] = useState<string | null>(null);

  // Fetch plate layout data
  const { data: layout, isLoading, error } = api.get<PlateLayoutRecord>(`plate-layouts/${id}/`);
  const { data: protocols } = api.get<Protocol[]>(`plate-layouts/${id}/protocols/`);

  const handleEdit = () => {
    if (!layout) return;
    setEditedName(layout.name);
    setEditedDescription(layout.description || '');
    setEditedConfig(layout.config);
    setSaveError(null);
    setIsEditing(true);
  };

  const handleCancel = () => {
    setIsEditing(false);
    setEditedName('');
    setEditedDescription('');
    setEditedConfig(null);
    setSaveError(null);
  };

  const handleSave = async () => {
    if (!layout || !editedConfig) return;

    setSaving(true);
    setSaveError(null);

    try {
      await api.patch(`plate-layouts/${id}/`, {
        name: editedName,
        description: editedDescription || null,
        config: editedConfig,
      });
      mutate(`/api/proxy/compounds/plate-layouts/${id}/`);
      mutate('/api/proxy/compounds/plate-layouts/');
      setIsEditing(false);
    } catch (err: unknown) {
      console.error('Save failed:', err);
      if (err instanceof Error) {
        setSaveError(err.message);
      } else {
        setSaveError('Failed to save plate layout');
      }
    } finally {
      setSaving(false);
    }
  };

  const handleDelete = async () => {
    if (!layout) return;

    setDeleting(true);
    setDeleteError(null);

    try {
      await api.delete(`plate-layouts/${id}/`);
      mutate('/api/proxy/compounds/plate-layouts/');
      router.push(routes.assays.plateLayouts());
    } catch (err: unknown) {
      console.error('Delete failed:', err);
      if (err instanceof Error) {
        setDeleteError(err.message);
      } else {
        setDeleteError('Failed to delete plate layout');
      }
    } finally {
      setDeleting(false);
    }
  };

  if (error) {
    return (
      <Container maxWidth="lg" sx={{ py: 3 }}>
        <Alert severity="error">
          Failed to load plate layout. It may have been deleted.
        </Alert>
        <Button
          sx={{ mt: 2 }}
          variant="outlined"
          onClick={() => router.push(routes.assays.plateLayouts())}
        >
          Back to Plate Layouts
        </Button>
      </Container>
    );
  }

  const canDelete = layout && !layout.is_builtin && (!protocols || protocols.length === 0);

  return (
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <PageHeader
        breadcrumbs={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Assays', href: routes.assays.protocols() },
          { label: 'Plate Layouts', href: routes.assays.plateLayouts() },
          { label: isLoading ? 'Loading...' : layout?.name || 'Not Found', icon: 'protocol' },
        ]}
      />

      {isLoading ? (
        <Box>
          <Skeleton variant="text" width={300} height={50} />
          <Skeleton variant="rectangular" height={400} sx={{ mt: 2 }} />
        </Box>
      ) : layout ? (
        <>
          {/* Header */}
          <Box sx={{ mb: 3, display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
            <Box>
              {isEditing ? (
                <TextField
                  value={editedName}
                  onChange={(e) => setEditedName(e.target.value)}
                  variant="standard"
                  sx={{ fontSize: '2rem', '& input': { fontSize: '2rem', fontWeight: 500 } }}
                  fullWidth
                />
              ) : (
                <Typography variant="h4" gutterBottom>
                  <GridOn sx={{ mr: 1, verticalAlign: 'middle' }} />
                  {layout.name}
                </Typography>
              )}
              <Box sx={{ display: 'flex', gap: 1, mt: 1 }}>
                {layout.plate_format && (
                  <Chip label={`${layout.plate_format}-well`} size="small" variant="outlined" />
                )}
                {layout.is_builtin && (
                  <Chip label="Built-in" size="small" color="secondary" variant="outlined" />
                )}
                {protocols && (
                  <Chip
                    label={`${protocols.length} protocol${protocols.length !== 1 ? 's' : ''}`}
                    size="small"
                    color="primary"
                    variant="outlined"
                  />
                )}
              </Box>
            </Box>
            <Box sx={{ display: 'flex', gap: 1 }}>
              <Button
                variant="outlined"
                startIcon={<ArrowBack />}
                onClick={() => router.push(routes.assays.plateLayouts())}
              >
                Back
              </Button>
              {isEditing ? (
                <>
                  <Button
                    variant="outlined"
                    startIcon={<Cancel />}
                    onClick={handleCancel}
                    disabled={saving}
                  >
                    Cancel
                  </Button>
                  <Button
                    variant="contained"
                    startIcon={<Save />}
                    onClick={handleSave}
                    disabled={saving || !editedName.trim()}
                  >
                    {saving ? 'Saving...' : 'Save'}
                  </Button>
                </>
              ) : (
                <>
                  {!layout.is_builtin && (
                    <Button
                      variant="outlined"
                      startIcon={<Edit />}
                      onClick={handleEdit}
                    >
                      Edit
                    </Button>
                  )}
                  <Button
                    variant="outlined"
                    color="error"
                    startIcon={<Delete />}
                    onClick={() => setDeleteDialogOpen(true)}
                    disabled={!canDelete}
                  >
                    Delete
                  </Button>
                </>
              )}
            </Box>
          </Box>

          {saveError && (
            <Alert severity="error" sx={{ mb: 2 }}>
              {saveError}
            </Alert>
          )}

          <Grid container spacing={3}>
            {/* Left column: Layout info */}
            <Grid size={{ xs: 12, md: 4 }}>
              <Paper sx={{ p: 2 }}>
                <Typography variant="h6" gutterBottom>
                  Layout Information
                </Typography>
                <Divider sx={{ my: 1 }} />

                <InfoRow
                  label="Format"
                  value={layout.plate_format ? `${layout.plate_format}-well plate` : '-'}
                />
                <InfoRow
                  label="Type"
                  value={layout.is_builtin ? 'Built-in' : 'Custom'}
                />
                <InfoRow label="Created By" value={layout.created_by_email || '-'} />
                <InfoRow
                  label="Created"
                  value={new Date(layout.created_at).toLocaleDateString()}
                />
                <InfoRow
                  label="Updated"
                  value={new Date(layout.updated_at).toLocaleDateString()}
                />

                <Box sx={{ mt: 2 }}>
                  <Typography variant="subtitle2" color="text.secondary" gutterBottom>
                    Description
                  </Typography>
                  {isEditing ? (
                    <TextField
                      value={editedDescription}
                      onChange={(e) => setEditedDescription(e.target.value)}
                      multiline
                      rows={3}
                      fullWidth
                      placeholder="Add a description..."
                    />
                  ) : (
                    <Typography variant="body2">
                      {layout.description || 'No description'}
                    </Typography>
                  )}
                </Box>

                {/* Protocols using this layout */}
                {protocols && protocols.length > 0 && (
                  <Box sx={{ mt: 3 }}>
                    <Typography variant="subtitle2" color="text.secondary" gutterBottom>
                      Protocols Using This Layout
                    </Typography>
                    <List dense disablePadding>
                      {protocols.map((protocol) => (
                        <ListItem key={protocol.id} disablePadding>
                          <ListItemButton
                            component={Link}
                            href={routes.assays.protocol(protocol.id)}
                          >
                            <ListItemText
                              primary={protocol.name}
                              secondary={protocol.import_type}
                            />
                          </ListItemButton>
                        </ListItem>
                      ))}
                    </List>
                  </Box>
                )}
              </Paper>
            </Grid>

            {/* Right column: Layout editor/preview */}
            <Grid size={{ xs: 12, md: 8 }}>
              <Paper sx={{ p: 2 }}>
                <Typography variant="h6" gutterBottom>
                  Plate Configuration
                </Typography>
                <Divider sx={{ my: 1 }} />

                {isEditing ? (
                  <PlateLayoutEditor
                    value={editedConfig || layout.config}
                    onChange={setEditedConfig}
                  />
                ) : (
                  <Box sx={{ display: 'flex', justifyContent: 'center', py: 2 }}>
                    <PlatePreview layout={layout.config} />
                  </Box>
                )}
              </Paper>
            </Grid>
          </Grid>

          {/* Delete confirmation dialog */}
          <Dialog open={deleteDialogOpen} onClose={() => setDeleteDialogOpen(false)}>
            <DialogTitle>Delete Plate Layout?</DialogTitle>
            <DialogContent>
              {deleteError && (
                <Alert severity="error" sx={{ mb: 2 }}>
                  {deleteError}
                </Alert>
              )}
              <Typography>
                Are you sure you want to delete this plate layout?
              </Typography>
              <Typography sx={{ fontWeight: 'medium', mt: 1 }}>
                {layout.name}
              </Typography>
              {protocols && protocols.length > 0 && (
                <Alert severity="warning" sx={{ mt: 2 }}>
                  This layout is used by {protocols.length} protocol(s) and cannot be deleted.
                </Alert>
              )}
              {layout.is_builtin && (
                <Alert severity="warning" sx={{ mt: 2 }}>
                  Built-in layouts cannot be deleted.
                </Alert>
              )}
            </DialogContent>
            <DialogActions>
              <Button onClick={() => setDeleteDialogOpen(false)} disabled={deleting}>
                Cancel
              </Button>
              <Button
                onClick={handleDelete}
                color="error"
                variant="contained"
                disabled={deleting || !canDelete}
              >
                {deleting ? 'Deleting...' : 'Delete'}
              </Button>
            </DialogActions>
          </Dialog>
        </>
      ) : (
        <Alert severity="warning">Plate layout not found</Alert>
      )}
    </Container>
  );
}
