'use client';

import { useState, useCallback, useEffect, useMemo, Suspense } from 'react';
import { useRouter, useSearchParams } from 'next/navigation';
import Link from 'next/link';
import {
  Container,
  Paper,
  Typography,
  Box,
  TextField,
  Button,
  Alert,
  CircularProgress,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  Autocomplete,
  Grid,
  Divider,
  Chip,
  Accordion,
  AccordionSummary,
  AccordionDetails,
  InputAdornment,
  IconButton,
} from '@mui/material';
import {
  Add,
  Science,
  Medication,
  ArrowBack,
  Save,
  Upload,
  ExpandMore,
  Draw,
  Clear,
} from '@mui/icons-material';
import { Breadcrumbs } from '@/components/Breadcrumbs';
import { JSMEEditor } from '@/components/JSMEEditor';
import { MoleculeChip } from '@/components/MoleculeView';
import { useCompoundsApi, apiPost } from '@/lib/api';

interface Target {
  id: string;
  name: string;
}

interface Supplier {
  id: string;
  name: string;
  initials: string | null;
}

interface CompoundFormData {
  target: string;
  smiles: string;
  stereo_comment: string;
  supplier: string | null;
  supplier_ref: string;
  comments: string;
  labbook_number: number | null;
  page_number: number | null;
}

const STEREO_OPTIONS = [
  { value: 'unset', label: 'Unset' },
  { value: 'achiral', label: 'Achiral' },
  { value: 'racemic', label: 'Racemic mixture' },
  { value: 'single_unknown', label: 'Single enantiomer, configuration unknown' },
  { value: 'r_enantiomer', label: 'R enantiomer' },
  { value: 's_enantiomer', label: 'S enantiomer' },
  { value: 'non_racemic_mixture', label: 'Non-racemic stereoisomer mixture' },
  { value: 'four_diastereomers', label: 'Mixture of 4 diastereoisomers' },
  { value: 'two_diastereomers', label: 'Mixture of 2 diastereoisomers' },
  { value: 'single_diastereomer_unknown', label: 'Single diastereoisomer, configuration unknown' },
  { value: 'rr_diastereomer', label: 'RR diastereoisomer' },
  { value: 'rs_diastereomer', label: 'RS diastereoisomer' },
  { value: 'sr_diastereomer', label: 'SR diastereoisomer' },
  { value: 'ss_diastereomer', label: 'SS diastereoisomer' },
  { value: 'epimer_mixture', label: 'Mixture of epimers' },
  { value: 'ez_mixture', label: 'Mixture of E and Z isomers' },
  { value: 'e_isomer', label: 'E isomer' },
  { value: 'z_isomer', label: 'Z isomer' },
];

export default function NewCompoundPage() {
  return (
    <Suspense fallback={<Container maxWidth="lg" sx={{ py: 3 }}><CircularProgress /></Container>}>
      <NewCompoundPageContent />
    </Suspense>
  );
}

function NewCompoundPageContent() {
  const router = useRouter();
  const searchParams = useSearchParams();
  const api = useCompoundsApi();

  // Get pre-populated target from URL (e.g., from target detail page)
  const preselectedTargetId = searchParams.get('target');

  const [formData, setFormData] = useState<CompoundFormData>({
    target: preselectedTargetId || '',
    smiles: '',
    stereo_comment: 'unset',
    supplier: null,
    supplier_ref: '',
    comments: '',
    labbook_number: null,
    page_number: null,
  });

  const [submitting, setSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [success, setSuccess] = useState<{ formatted_id: string; id: string } | null>(null);
  const [sketcherExpanded, setSketcherExpanded] = useState(false);
  const [jsmeInitialSmiles, setJsmeInitialSmiles] = useState('');

  // Fetch targets and suppliers
  const { data: targetsData } = api.get<Target[]>('targets/');
  const { data: suppliersData } = api.get<Supplier[]>('suppliers/');

  const targets = targetsData || [];
  const suppliers = suppliersData || [];

  // Get pre-selected target info for display
  const preselectedTarget = useMemo(() => {
    if (!preselectedTargetId || !targets.length) return null;
    return targets.find((t) => t.id === preselectedTargetId) || null;
  }, [preselectedTargetId, targets]);

  // Update form data when targets load and we have a preselected target
  useEffect(() => {
    if (preselectedTargetId && targets.length > 0 && !formData.target) {
      setFormData((prev) => ({ ...prev, target: preselectedTargetId }));
    }
  }, [preselectedTargetId, targets, formData.target]);

  const handleFieldChange = useCallback((field: keyof CompoundFormData, value: any) => {
    setFormData((prev) => ({ ...prev, [field]: value }));
    setError(null);
  }, []);

  const handleSmilesChange = useCallback((smiles: string) => {
    setFormData((prev) => ({ ...prev, smiles }));
    setError(null);
  }, []);

  // Handle opening the sketcher - load current SMILES into JSME
  const handleSketcherExpand = useCallback((expanded: boolean) => {
    setSketcherExpanded(expanded);
    if (expanded && formData.smiles) {
      setJsmeInitialSmiles(formData.smiles);
    }
  }, [formData.smiles]);

  // Clear the SMILES field
  const handleClearSmiles = useCallback(() => {
    setFormData((prev) => ({ ...prev, smiles: '' }));
    setJsmeInitialSmiles('');
  }, []);

  const handleSubmit = useCallback(async (e: React.FormEvent) => {
    e.preventDefault();
    setError(null);

    // Validate required fields
    if (!formData.target) {
      setError('Please select a target');
      return;
    }
    if (!formData.smiles || formData.smiles.trim().length === 0) {
      setError('Please draw or enter a structure');
      return;
    }

    setSubmitting(true);

    try {
      // Prepare data - remove null/empty values
      const submitData: Record<string, any> = {
        target: formData.target,
        smiles: formData.smiles.trim(),
      };

      if (formData.stereo_comment && formData.stereo_comment !== 'unset') {
        submitData.stereo_comment = formData.stereo_comment;
      }
      if (formData.supplier) {
        submitData.supplier = formData.supplier;
      }
      if (formData.supplier_ref.trim()) {
        submitData.supplier_ref = formData.supplier_ref.trim();
      }
      if (formData.comments.trim()) {
        submitData.comments = formData.comments.trim();
      }
      if (formData.labbook_number !== null) {
        submitData.labbook_number = formData.labbook_number;
      }
      if (formData.page_number !== null) {
        submitData.page_number = formData.page_number;
      }

      const result = await apiPost<{ id: string; formatted_id: string }>(
        'registry/compounds/',
        submitData
      );

      setSuccess(result);

    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to register compound');
    } finally {
      setSubmitting(false);
    }
  }, [formData]);

  const handleRegisterAnother = useCallback(() => {
    setSuccess(null);
    setFormData({
      target: formData.target, // Keep the target
      smiles: '',
      stereo_comment: 'unset',
      supplier: formData.supplier, // Keep the supplier
      supplier_ref: '',
      comments: '',
      labbook_number: formData.labbook_number, // Keep lab notebook info
      page_number: null,
    });
  }, [formData.target, formData.supplier, formData.labbook_number]);

  // Show success screen
  if (success) {
    return (
      <Container maxWidth="md" sx={{ py: 4 }}>
        <Paper sx={{ p: 4, textAlign: 'center' }}>
          <Medication sx={{ fontSize: 64, color: 'success.main', mb: 2 }} />
          <Typography variant="h5" gutterBottom>
            Compound Registered Successfully
          </Typography>
          <Chip
            label={success.formatted_id}
            color="success"
            size="medium"
            sx={{ fontSize: '1.2rem', py: 2, mb: 3 }}
          />
          <Box sx={{ display: 'flex', justifyContent: 'center', gap: 2 }}>
            <Button
              variant="outlined"
              onClick={handleRegisterAnother}
              startIcon={<Add />}
            >
              Register Another
            </Button>
            <Button
              variant="contained"
              component={Link}
              href={`/registry/compounds/${success.id}`}
            >
              View Compound
            </Button>
          </Box>
        </Paper>
      </Container>
    );
  }

  return (
    <Container maxWidth="lg" sx={{ py: 4 }}>
      <Breadcrumbs
        items={
          preselectedTargetId && preselectedTarget
            ? [
                { label: 'Registry', href: '/registry/targets' },
                { label: preselectedTarget.name, href: `/registry/targets/${preselectedTargetId}`, icon: 'target' },
                { label: 'Register Compound' },
              ]
            : [
                { label: 'Registry', href: '/registry' },
                { label: 'Register Compound' },
              ]
        }
      />

      <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 3 }}>
        <Button
          component={Link}
          href={preselectedTargetId ? `/registry/targets/${preselectedTargetId}` : '/registry/search'}
          startIcon={<ArrowBack />}
          size="small"
        >
          {preselectedTargetId ? 'Back to Target' : 'Back to Search'}
        </Button>
        <Typography variant="h4" sx={{ flex: 1 }}>
          Register New Compound
        </Typography>
        <Button
          component={Link}
          href="/registry/import"
          variant="outlined"
          startIcon={<Upload />}
        >
          Bulk Import
        </Button>
      </Box>

      <form onSubmit={handleSubmit}>
        <Grid container spacing={3}>
          {/* Left column - Structure input */}
          <Grid item xs={12} md={6}>
            <Paper sx={{ p: 3 }}>
              <Typography variant="h6" gutterBottom sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                <Science color="primary" />
                Structure
              </Typography>

              {/* SMILES text input */}
              <TextField
                label="SMILES *"
                value={formData.smiles}
                onChange={(e) => handleSmilesChange(e.target.value)}
                fullWidth
                multiline
                rows={2}
                placeholder="Enter SMILES string, e.g., CC(=O)Oc1ccccc1C(=O)O"
                InputProps={{
                  endAdornment: formData.smiles && (
                    <InputAdornment position="end">
                      <IconButton onClick={handleClearSmiles} size="small" edge="end">
                        <Clear fontSize="small" />
                      </IconButton>
                    </InputAdornment>
                  ),
                  sx: { fontFamily: 'monospace' },
                }}
                sx={{ mb: 2 }}
              />

              {/* Structure preview */}
              <Box
                sx={{
                  display: 'flex',
                  justifyContent: 'center',
                  alignItems: 'center',
                  minHeight: 180,
                  bgcolor: 'grey.50',
                  borderRadius: 1,
                  border: '1px solid',
                  borderColor: 'divider',
                  mb: 2,
                }}
              >
                {formData.smiles ? (
                  <MoleculeChip smiles={formData.smiles} size={160} />
                ) : (
                  <Typography color="text.secondary">
                    Structure preview will appear here
                  </Typography>
                )}
              </Box>

              {/* JSME Sketcher in accordion */}
              <Accordion
                expanded={sketcherExpanded}
                onChange={(_, expanded) => handleSketcherExpand(expanded)}
                sx={{ bgcolor: 'grey.50' }}
              >
                <AccordionSummary expandIcon={<ExpandMore />}>
                  <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                    <Draw fontSize="small" />
                    <Typography>Draw Structure</Typography>
                  </Box>
                </AccordionSummary>
                <AccordionDetails sx={{ bgcolor: 'white' }}>
                  <Typography variant="caption" color="text.secondary" paragraph>
                    Use the sketcher to draw a structure. The SMILES will be updated automatically.
                  </Typography>
                  <JSMEEditor
                    id="new-compound-jsme"
                    onChange={handleSmilesChange}
                    editable={true}
                    initialSmiles={jsmeInitialSmiles}
                    showPreview={false}
                    width={380}
                    height={320}
                  />
                </AccordionDetails>
              </Accordion>
            </Paper>
          </Grid>

          {/* Right column - Form fields */}
          <Grid item xs={12} md={6}>
            <Paper sx={{ p: 3 }}>
              <Typography variant="h6" gutterBottom sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                <Medication color="secondary" />
                Compound Details
              </Typography>

              {/* Target - Required, read-only if pre-selected from URL */}
              {preselectedTargetId ? (
                <TextField
                  label="Target *"
                  value={preselectedTarget?.name || 'Loading...'}
                  margin="normal"
                  fullWidth
                  disabled
                  InputProps={{
                    startAdornment: <Science sx={{ mr: 1, color: 'primary.main' }} />,
                  }}
                  helperText="Target is locked when registering from target page"
                />
              ) : (
                <Autocomplete
                  options={targets}
                  getOptionLabel={(option) => option.name}
                  value={targets.find((t) => t.id === formData.target) || null}
                  onChange={(_, newValue) => handleFieldChange('target', newValue?.id || '')}
                  renderInput={(params) => (
                    <TextField
                      {...params}
                      label="Target *"
                      margin="normal"
                      fullWidth
                      required
                    />
                  )}
                />
              )}

              {/* Stereo Comment */}
              <FormControl fullWidth margin="normal">
                <InputLabel>Stereochemistry</InputLabel>
                <Select
                  value={formData.stereo_comment}
                  label="Stereochemistry"
                  onChange={(e) => handleFieldChange('stereo_comment', e.target.value)}
                >
                  {STEREO_OPTIONS.map((option) => (
                    <MenuItem key={option.value} value={option.value}>
                      {option.label}
                    </MenuItem>
                  ))}
                </Select>
              </FormControl>

              <Divider sx={{ my: 2 }} />

              {/* Supplier */}
              <Autocomplete
                options={suppliers}
                getOptionLabel={(option) => option.name}
                value={suppliers.find((s) => s.id === formData.supplier) || null}
                onChange={(_, newValue) => handleFieldChange('supplier', newValue?.id || null)}
                renderInput={(params) => (
                  <TextField
                    {...params}
                    label="Supplier"
                    margin="normal"
                    fullWidth
                  />
                )}
              />

              {/* Supplier Reference */}
              <TextField
                label="Supplier Reference"
                value={formData.supplier_ref}
                onChange={(e) => handleFieldChange('supplier_ref', e.target.value)}
                margin="normal"
                fullWidth
                placeholder="Catalog number or order ID"
              />

              <Divider sx={{ my: 2 }} />

              {/* Lab notebook info */}
              <Grid container spacing={2}>
                <Grid item xs={6}>
                  <TextField
                    label="Lab Notebook #"
                    type="number"
                    value={formData.labbook_number ?? ''}
                    onChange={(e) => handleFieldChange('labbook_number', e.target.value ? parseInt(e.target.value) : null)}
                    margin="normal"
                    fullWidth
                  />
                </Grid>
                <Grid item xs={6}>
                  <TextField
                    label="Page #"
                    type="number"
                    value={formData.page_number ?? ''}
                    onChange={(e) => handleFieldChange('page_number', e.target.value ? parseInt(e.target.value) : null)}
                    margin="normal"
                    fullWidth
                  />
                </Grid>
              </Grid>

              {/* Comments */}
              <TextField
                label="Comments"
                value={formData.comments}
                onChange={(e) => handleFieldChange('comments', e.target.value)}
                margin="normal"
                fullWidth
                multiline
                rows={3}
              />

              {error && (
                <Alert severity="error" sx={{ mt: 2 }}>
                  {error}
                </Alert>
              )}

              {/* Submit button */}
              <Box sx={{ mt: 3, display: 'flex', justifyContent: 'flex-end', gap: 2 }}>
                <Button
                  component={Link}
                  href="/registry/search"
                  variant="outlined"
                >
                  Cancel
                </Button>
                <Button
                  type="submit"
                  variant="contained"
                  disabled={submitting || !formData.target || !formData.smiles}
                  startIcon={submitting ? <CircularProgress size={20} /> : <Save />}
                >
                  {submitting ? 'Registering...' : 'Register Compound'}
                </Button>
              </Box>
            </Paper>
          </Grid>
        </Grid>
      </form>
    </Container>
  );
}
