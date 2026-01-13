'use client';

import { useState, useEffect, useCallback } from 'react';
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  TextField,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  Box,
  Alert,
  IconButton,
  CircularProgress,
  Typography,
  Paper,
  Chip,
  Divider,
  List,
  ListItem,
  ListItemText,
  Collapse,
} from '@mui/material';
import {
  Close,
  Science,
  Add,
  CloudUpload,
  Description,
  Delete,
  ExpandMore,
  ExpandLess,
  CheckCircle,
} from '@mui/icons-material';
import { useCompoundsApi } from '@/lib/compounds/api';
import type { Plasmid, ConstructProject } from '@/types/compounds/constructs';

interface GenbankFeature {
  type: string;
  name: string;
  start: number;
  end: number;
  strand: number;
}

interface GenbankPreview {
  name: string;
  length: number;
  features: GenbankFeature[];
  definition?: string;
  organism?: string;
}

/**
 * Simple GenBank parser to extract key metadata for preview.
 * Parses LOCUS, DEFINITION, ORGANISM, and FEATURES sections.
 */
function parseGenbankFile(content: string): GenbankPreview | null {
  try {
    const lines = content.split('\n');
    const preview: GenbankPreview = {
      name: '',
      length: 0,
      features: [],
    };

    let inFeatures = false;
    let currentFeature: Partial<GenbankFeature> | null = null;

    for (let i = 0; i < lines.length; i++) {
      const line = lines[i];

      // Parse LOCUS line
      if (line.startsWith('LOCUS')) {
        const parts = line.split(/\s+/);
        preview.name = parts[1] || '';
        // Try to extract length (usually "1234 bp")
        const bpMatch = line.match(/(\d+)\s+bp/);
        if (bpMatch) {
          preview.length = parseInt(bpMatch[1], 10);
        }
      }

      // Parse DEFINITION
      if (line.startsWith('DEFINITION')) {
        preview.definition = line.substring(12).trim();
        // Continue reading if definition spans multiple lines
        let j = i + 1;
        while (j < lines.length && lines[j].startsWith('            ')) {
          preview.definition += ' ' + lines[j].trim();
          j++;
        }
      }

      // Parse ORGANISM
      if (line.startsWith('  ORGANISM')) {
        preview.organism = line.substring(12).trim();
      }

      // Track FEATURES section
      if (line.startsWith('FEATURES')) {
        inFeatures = true;
        continue;
      }

      // End of FEATURES section
      if (line.startsWith('ORIGIN') || line.startsWith('CONTIG')) {
        inFeatures = false;
        if (currentFeature && currentFeature.type && currentFeature.name) {
          preview.features.push(currentFeature as GenbankFeature);
        }
        continue;
      }

      // Parse features
      if (inFeatures) {
        // New feature line (starts with 5 spaces then feature type)
        const featureMatch = line.match(/^     (\w+)\s+(.+)/);
        if (featureMatch) {
          // Save previous feature
          if (currentFeature && currentFeature.type && currentFeature.name) {
            preview.features.push(currentFeature as GenbankFeature);
          }

          currentFeature = {
            type: featureMatch[1],
            name: '',
            start: 0,
            end: 0,
            strand: 1,
          };

          // Parse location
          const locationStr = featureMatch[2];
          const complementMatch = locationStr.match(/complement\((\d+)\.\.(\d+)\)/);
          const simpleMatch = locationStr.match(/(\d+)\.\.(\d+)/);

          if (complementMatch) {
            currentFeature.start = parseInt(complementMatch[1], 10);
            currentFeature.end = parseInt(complementMatch[2], 10);
            currentFeature.strand = -1;
          } else if (simpleMatch) {
            currentFeature.start = parseInt(simpleMatch[1], 10);
            currentFeature.end = parseInt(simpleMatch[2], 10);
          }
        }

        // Parse qualifier lines (starts with 21 spaces then /)
        if (currentFeature && line.match(/^\s{21}\//)) {
          const labelMatch = line.match(/\/label="?([^"]+)"?/);
          const geneMatch = line.match(/\/gene="?([^"]+)"?/);
          const productMatch = line.match(/\/product="?([^"]+)"?/);
          const noteMatch = line.match(/\/note="?([^"]+)"?/);

          if (labelMatch && !currentFeature.name) {
            currentFeature.name = labelMatch[1];
          } else if (geneMatch && !currentFeature.name) {
            currentFeature.name = geneMatch[1];
          } else if (productMatch && !currentFeature.name) {
            currentFeature.name = productMatch[1];
          } else if (noteMatch && !currentFeature.name) {
            currentFeature.name = noteMatch[1].substring(0, 50);
          }
        }
      }
    }

    // Assign default names to unnamed features
    preview.features = preview.features.map((f, i) => ({
      ...f,
      name: f.name || `${f.type}_${i + 1}`,
    }));

    return preview;
  } catch (e) {
    console.error('Failed to parse GenBank file:', e);
    return null;
  }
}

interface PlasmidCreateDialogProps {
  open: boolean;
  onClose: () => void;
  onCreated: (plasmid: Plasmid) => void;
}

export function PlasmidCreateDialog({
  open,
  onClose,
  onCreated,
}: PlasmidCreateDialogProps) {
  const api = useCompoundsApi();
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Form state
  const [name, setName] = useState('');
  const [projectId, setProjectId] = useState<string | null>(null);
  const [parentId, setParentId] = useState<string | null>(null);
  const [genbankFile, setGenbankFile] = useState<File | null>(null);
  const [genbankPreview, setGenbankPreview] = useState<GenbankPreview | null>(null);
  const [showFeatures, setShowFeatures] = useState(false);

  // Fetch projects and plasmids for dropdowns
  const { data: projects, isLoading: projectsLoading } = api.get<ConstructProject[]>(
    open ? 'construct-projects/' : null
  );
  const { data: plasmids, isLoading: plasmidsLoading } = api.get<Plasmid[]>(
    open ? 'plasmids/' : null
  );

  // Reset form when dialog opens
  useEffect(() => {
    if (open) {
      setName('');
      setProjectId(null);
      setParentId(null);
      setGenbankFile(null);
      setGenbankPreview(null);
      setShowFeatures(false);
      setError(null);
    }
  }, [open]);

  const handleDrop = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    const file = e.dataTransfer.files[0];
    if (file && isValidGenbankFile(file)) {
      processGenbankFile(file);
    } else {
      setError('Please upload a GenBank file (.gb or .gbk)');
    }
  }, []);

  const handleDragOver = useCallback((e: React.DragEvent) => {
    e.preventDefault();
  }, []);

  const handleFileSelect = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (file) {
      if (isValidGenbankFile(file)) {
        processGenbankFile(file);
      } else {
        setError('Please upload a GenBank file (.gb or .gbk)');
      }
    }
  };

  const isValidGenbankFile = (file: File): boolean => {
    const fileName = file.name.toLowerCase();
    return fileName.endsWith('.gb') || fileName.endsWith('.gbk');
  };

  const processGenbankFile = async (file: File) => {
    setGenbankFile(file);
    setError(null);

    try {
      const content = await file.text();
      const preview = parseGenbankFile(content);

      if (preview) {
        setGenbankPreview(preview);
        // Auto-fill name from GenBank if name is empty
        if (!name && preview.name) {
          setName(preview.name);
        }
      } else {
        setError('Could not parse GenBank file. The file may be malformed.');
      }
    } catch (e) {
      setError('Failed to read file');
    }
  };

  const handleSave = async () => {
    if (!name.trim()) {
      setError('Plasmid name is required');
      return;
    }

    setSaving(true);
    setError(null);

    try {
      const formData = new FormData();
      formData.append('name', name.trim());
      if (projectId) {
        formData.append('project', projectId);
      }
      if (parentId) {
        formData.append('parent', parentId);
      }
      if (genbankFile) {
        formData.append('genbank_file', genbankFile);
      }

      const newPlasmid = await api.upload<Plasmid>('plasmids/', formData);
      onCreated(newPlasmid);
      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to create plasmid');
    } finally {
      setSaving(false);
    }
  };

  const featureTypeColors: Record<string, string> = {
    CDS: 'primary',
    gene: 'success',
    promoter: 'warning',
    terminator: 'error',
    misc_feature: 'info',
    rep_origin: 'secondary',
  };

  return (
    <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
      <DialogTitle
        sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}
      >
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Science color="primary" />
          New Plasmid Construct
        </Box>
        <IconButton onClick={onClose} size="small" disabled={saving}>
          <Close />
        </IconButton>
      </DialogTitle>

      <DialogContent dividers>
        {error && (
          <Alert severity="error" sx={{ mb: 2 }}>
            {error}
          </Alert>
        )}

        <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2.5, pt: 1 }}>
          {/* GenBank file upload - moved to top for workflow */}
          <Box>
            <Typography variant="subtitle2" gutterBottom>
              GenBank File
            </Typography>
            {genbankFile ? (
              <Paper variant="outlined" sx={{ p: 2 }}>
                <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', mb: 1 }}>
                  <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                    <CheckCircle color="success" fontSize="small" />
                    <Typography variant="body2" fontWeight={500}>
                      {genbankFile.name}
                    </Typography>
                    <Typography variant="caption" color="text.secondary">
                      ({(genbankFile.size / 1024).toFixed(1)} KB)
                    </Typography>
                  </Box>
                  <IconButton
                    size="small"
                    onClick={() => {
                      setGenbankFile(null);
                      setGenbankPreview(null);
                    }}
                    disabled={saving}
                  >
                    <Delete fontSize="small" />
                  </IconButton>
                </Box>

                {/* GenBank Preview */}
                {genbankPreview && (
                  <Box sx={{ bgcolor: 'grey.50', p: 1.5, borderRadius: 1 }}>
                    <Typography variant="caption" color="text.secondary" display="block">
                      <strong>Locus:</strong> {genbankPreview.name}
                    </Typography>
                    <Typography variant="caption" color="text.secondary" display="block">
                      <strong>Length:</strong> {genbankPreview.length.toLocaleString()} bp
                    </Typography>
                    {genbankPreview.definition && (
                      <Typography variant="caption" color="text.secondary" display="block" noWrap>
                        <strong>Definition:</strong> {genbankPreview.definition}
                      </Typography>
                    )}
                    {genbankPreview.organism && (
                      <Typography variant="caption" color="text.secondary" display="block">
                        <strong>Organism:</strong> {genbankPreview.organism}
                      </Typography>
                    )}

                    {/* Features summary */}
                    {genbankPreview.features.length > 0 && (
                      <Box sx={{ mt: 1 }}>
                        <Button
                          size="small"
                          onClick={() => setShowFeatures(!showFeatures)}
                          endIcon={showFeatures ? <ExpandLess /> : <ExpandMore />}
                          sx={{ p: 0, minWidth: 0 }}
                        >
                          {genbankPreview.features.length} features
                        </Button>
                        <Collapse in={showFeatures}>
                          <List dense sx={{ maxHeight: 150, overflow: 'auto', mt: 0.5 }}>
                            {genbankPreview.features.slice(0, 20).map((feature, idx) => (
                              <ListItem key={idx} sx={{ py: 0, px: 0.5 }}>
                                <Chip
                                  label={feature.type}
                                  size="small"
                                  color={(featureTypeColors[feature.type] as any) || 'default'}
                                  sx={{ mr: 1, minWidth: 70, fontSize: '0.65rem' }}
                                />
                                <ListItemText
                                  primary={feature.name}
                                  secondary={`${feature.start}..${feature.end}${feature.strand === -1 ? ' (-)' : ''}`}
                                  primaryTypographyProps={{ variant: 'caption' }}
                                  secondaryTypographyProps={{ variant: 'caption' }}
                                />
                              </ListItem>
                            ))}
                            {genbankPreview.features.length > 20 && (
                              <ListItem sx={{ py: 0 }}>
                                <Typography variant="caption" color="text.secondary">
                                  ... and {genbankPreview.features.length - 20} more
                                </Typography>
                              </ListItem>
                            )}
                          </List>
                        </Collapse>
                      </Box>
                    )}
                  </Box>
                )}
              </Paper>
            ) : (
              <Paper
                variant="outlined"
                onDrop={handleDrop}
                onDragOver={handleDragOver}
                sx={{
                  p: 3,
                  textAlign: 'center',
                  bgcolor: 'grey.50',
                  borderStyle: 'dashed',
                  cursor: 'pointer',
                  '&:hover': { bgcolor: 'grey.100' },
                }}
                onClick={() => document.getElementById('genbank-file-input')?.click()}
              >
                <CloudUpload sx={{ fontSize: 36, color: 'grey.400', mb: 1 }} />
                <Typography color="text.secondary" variant="body2">
                  Drag and drop a GenBank file, or click to select
                </Typography>
                <Box sx={{ mt: 1 }}>
                  <Chip label=".gb" size="small" sx={{ mr: 0.5 }} />
                  <Chip label=".gbk" size="small" />
                </Box>
                <input
                  id="genbank-file-input"
                  type="file"
                  accept=".gb,.gbk"
                  onChange={handleFileSelect}
                  style={{ display: 'none' }}
                />
              </Paper>
            )}
          </Box>

          <Divider />

          {/* Name field */}
          <TextField
            label="Plasmid Name"
            value={name}
            onChange={(e) => setName(e.target.value)}
            fullWidth
            required
            size="small"
            placeholder="e.g., pET28a-EGFP, pGEX-6P-1"
            helperText={
              genbankPreview?.name && name === genbankPreview.name
                ? 'Auto-filled from GenBank file'
                : 'A descriptive name for this plasmid construct'
            }
          />

          {/* Project selection */}
          <FormControl fullWidth size="small">
            <InputLabel>Project (optional)</InputLabel>
            <Select
              value={projectId || ''}
              onChange={(e) => setProjectId(e.target.value || null)}
              label="Project (optional)"
            >
              <MenuItem value="">
                <em>None</em>
              </MenuItem>
              {projectsLoading ? (
                <MenuItem disabled>Loading...</MenuItem>
              ) : (
                projects?.map((project) => (
                  <MenuItem key={project.id} value={project.id}>
                    {project.name}
                  </MenuItem>
                ))
              )}
            </Select>
          </FormControl>

          {/* Parent plasmid selection */}
          <FormControl fullWidth size="small">
            <InputLabel>Parent Plasmid (optional)</InputLabel>
            <Select
              value={parentId || ''}
              onChange={(e) => setParentId(e.target.value || null)}
              label="Parent Plasmid (optional)"
            >
              <MenuItem value="">
                <em>None</em>
              </MenuItem>
              {plasmidsLoading ? (
                <MenuItem disabled>Loading...</MenuItem>
              ) : (
                plasmids?.map((plasmid) => (
                  <MenuItem key={plasmid.id} value={plasmid.id}>
                    {plasmid.formatted_id} - {plasmid.name}
                  </MenuItem>
                ))
              )}
            </Select>
          </FormControl>
        </Box>
      </DialogContent>

      <DialogActions>
        <Button onClick={onClose} disabled={saving}>
          Cancel
        </Button>
        <Button
          variant="contained"
          onClick={handleSave}
          disabled={saving || !name.trim()}
          startIcon={saving ? <CircularProgress size={16} /> : <Add />}
        >
          {saving ? 'Creating...' : 'Create Plasmid'}
        </Button>
      </DialogActions>
    </Dialog>
  );
}
