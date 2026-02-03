'use client';

import { Suspense, use, useState, useCallback, useEffect } from 'react';
import { useRouter, useSearchParams } from 'next/navigation';
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
  Dialog,
  DialogTitle,
  DialogContent,
  DialogContentText,
  DialogActions,
  IconButton,
  Alert,
  CircularProgress,
  LinearProgress,
  Tooltip,
  Accordion,
  AccordionSummary,
  AccordionDetails,
} from '@mui/material';
import { Description, Science, Assessment, Edit, GridOn, Close, Add, Delete, CloudUpload, Download, OpenInNew, TableChart, CheckCircle, Biotech, ExpandMore, FiberNew } from '@mui/icons-material';
import Link from 'next/link';
import { DetailPageLayout } from '@/components/compounds/DetailPageLayout';
import { DataTable, Column } from '@/components/data-table';
import { PlatePreview } from '@/components/compounds/PlatePreview';
import { PlateLayoutCreateDialog } from '@/components/compounds/PlateLayoutCreateDialog';
import { AssayUploadDrawer } from '@/components/compounds/AssayUploadDrawer';
import { ProtocolEditDialog } from '@/components/compounds/ProtocolEditDialog';
import { useCompoundsApi, apiUpload, getAuthenticatedDownloadUrl } from '@/lib/compounds/api';
import { useAuth } from '@/lib/compounds/auth-context';
import { routes } from '@/lib/compounds/routes';
import { Protocol, Assay, ProtocolDocument, PlateLayoutRecord, ImportType } from '@/types/compounds/models';

interface UploadingFile {
  name: string;
  status: 'uploading' | 'success' | 'error';
  error?: string;
}

interface PageProps {
  params: Promise<{ id: string }>;
}

const IMPORT_TYPE_LABELS: Record<ImportType, string> = {
  raw_data: 'Raw Data (Dose-Response)',
  ms_intact: 'MS-Intact',
  table_of_values: 'Table of Values',
  pharmaron_adme: 'Pharmaron ADME',
};

function InfoRow({ label, value }: { label: string; value: React.ReactNode }) {
  return (
    <Box sx={{ display: 'flex', py: 0.5, minWidth: 0 }}>
      <Typography
        color="text.secondary"
        sx={{ minWidth: 140, flexShrink: 0, fontWeight: 500 }}
      >
        {label}:
      </Typography>
      <Typography component="div" sx={{ minWidth: 0, wordBreak: 'break-word' }}>
        {value ?? '-'}
      </Typography>
    </Box>
  );
}

// Helper to check if import type is Pharmaron ADME
function isAdmeProtocol(importType: ImportType | undefined): boolean {
  return importType === 'pharmaron_adme';
}

// Helper to check if import type is table of values
function isTableOfValuesProtocol(importType: ImportType | undefined): boolean {
  return importType === 'table_of_values';
}

// Helper to check if protocol is plate-based (not ADME or table of values)
function isPlateBasedProtocol(importType: ImportType | undefined): boolean {
  return !isAdmeProtocol(importType) && !isTableOfValuesProtocol(importType);
}

function ProtocolDetailPageContent({ params }: PageProps) {
  const { id } = use(params);
  const router = useRouter();
  const searchParams = useSearchParams();
  const api = useCompoundsApi();
  const { canContribute } = useAuth();

  // Check for openUpload and target query params (from import page redirect)
  const shouldOpenUpload = searchParams.get('openUpload') === 'true';
  const targetIdFromUrl = searchParams.get('target');

  const [layoutSelectOpen, setLayoutSelectOpen] = useState(false);
  const [createLayoutOpen, setCreateLayoutOpen] = useState(false);
  const [saveError, setSaveError] = useState<string | null>(null);
  const [saving, setSaving] = useState(false);
  const [uploadDrawerOpen, setUploadDrawerOpen] = useState(false);
  const [hasAutoOpened, setHasAutoOpened] = useState(false);
  const [deleteDialogOpen, setDeleteDialogOpen] = useState(false);
  const [assayToDelete, setAssayToDelete] = useState<Assay | null>(null);
  const [deleting, setDeleting] = useState(false);
  const [editDialogOpen, setEditDialogOpen] = useState(false);

  // Document drop zone state
  const [isDragOver, setIsDragOver] = useState(false);
  const [uploadingFiles, setUploadingFiles] = useState<UploadingFile[]>([]);
  const [uploadError, setUploadError] = useState<string | null>(null);

  // Document delete state
  const [deleteDocDialogOpen, setDeleteDocDialogOpen] = useState(false);
  const [docToDelete, setDocToDelete] = useState<ProtocolDocument | null>(null);
  const [deletingDoc, setDeletingDoc] = useState(false);

  const { data: protocol, isLoading: protocolLoading, mutate } = api.get<Protocol>(
    `protocols/${id}/`
  );
  const { data: assays, isLoading: assaysLoading, mutate: mutateAssays } = api.get<Assay[]>(
    `assays/?protocol=${id}`
  );
  const { data: documents, isLoading: documentsLoading, mutate: mutateDocuments } = api.get<ProtocolDocument[]>(
    `protocol-documents/?protocol=${id}`
  );

  // Fetch available plate layouts for selection
  const { data: plateLayouts } = api.get<PlateLayoutRecord[]>('plate-layouts/');

  const loading = protocolLoading;

  // Auto-open upload drawer when redirected from import page with openUpload=true
  useEffect(() => {
    if (shouldOpenUpload && !hasAutoOpened && protocol && isPlateBasedProtocol(protocol.import_type)) {
      setUploadDrawerOpen(true);
      setHasAutoOpened(true);
      // Clean up URL params without triggering navigation
      const url = new URL(window.location.href);
      url.searchParams.delete('openUpload');
      url.searchParams.delete('target');
      window.history.replaceState({}, '', url.toString());
    }
  }, [shouldOpenUpload, hasAutoOpened, protocol]);

  const handleLayoutChange = async (layoutId: string | null) => {
    setSaving(true);
    setSaveError(null);

    try {
      await api.patch(`protocols/${id}/`, { plate_layout: layoutId });
      await mutate();
      setLayoutSelectOpen(false);
    } catch (err) {
      setSaveError(err instanceof Error ? err.message : 'Failed to save');
    } finally {
      setSaving(false);
    }
  };

  const handleLayoutCreated = async (newLayout: PlateLayoutRecord) => {
    // Set the newly created layout on this protocol
    await handleLayoutChange(newLayout.id);
    setCreateLayoutOpen(false);
  };

  const handleDeleteClick = (assay: Assay, e: React.MouseEvent) => {
    e.stopPropagation();
    setAssayToDelete(assay);
    setDeleteDialogOpen(true);
  };

  const handleDeleteConfirm = async () => {
    if (!assayToDelete) return;

    setDeleting(true);
    try {
      await api.delete(`assays/${assayToDelete.id}/`);
      mutateAssays();
    } catch (err) {
      console.error('Delete error:', err);
    } finally {
      setDeleting(false);
      setDeleteDialogOpen(false);
      setAssayToDelete(null);
    }
  };

  const handleDeleteDocClick = (doc: ProtocolDocument) => {
    setDocToDelete(doc);
    setDeleteDocDialogOpen(true);
  };

  const handleDeleteDocConfirm = async () => {
    if (!docToDelete) return;

    setDeletingDoc(true);
    try {
      await api.delete(`protocol-documents/${docToDelete.id}/`);
      mutateDocuments();
    } catch (err) {
      console.error('Delete error:', err);
      setUploadError('Failed to delete document');
    } finally {
      setDeletingDoc(false);
      setDeleteDocDialogOpen(false);
      setDocToDelete(null);
    }
  };

  // Handle drag and drop for document uploads
  const handleDragOver = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    setIsDragOver(true);
  }, []);

  const handleDragLeave = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    setIsDragOver(false);
  }, []);

  const handleDrop = useCallback(async (e: React.DragEvent) => {
    e.preventDefault();
    setIsDragOver(false);
    setUploadError(null);

    const droppedFiles = Array.from(e.dataTransfer.files);
    if (droppedFiles.length === 0) return;

    // Initialize upload status for all files
    const initialStatus: UploadingFile[] = droppedFiles.map((f) => ({
      name: f.name,
      status: 'uploading',
    }));
    setUploadingFiles(initialStatus);

    let successCount = 0;
    let errorCount = 0;

    // Upload files serially
    for (let i = 0; i < droppedFiles.length; i++) {
      const file = droppedFiles[i];
      try {
        const formData = new FormData();
        formData.append('file', file);
        formData.append('protocol', id);

        await apiUpload('protocol-documents/', formData);

        setUploadingFiles((prev) =>
          prev.map((f, idx) =>
            idx === i ? { ...f, status: 'success' } : f
          )
        );
        successCount++;
      } catch (err) {
        setUploadingFiles((prev) =>
          prev.map((f, idx) =>
            idx === i
              ? { ...f, status: 'error', error: err instanceof Error ? err.message : 'Upload failed' }
              : f
          )
        );
        errorCount++;
      }
    }

    // Refresh the documents list
    if (successCount > 0) {
      mutateDocuments();
    }

    // Clear upload status after a delay
    setTimeout(() => {
      setUploadingFiles([]);
      if (errorCount > 0) {
        setUploadError(`${errorCount} file(s) failed to upload`);
      }
    }, 3000);
  }, [id, mutateDocuments]);

  const handleFileSelect = useCallback(async (e: React.ChangeEvent<HTMLInputElement>) => {
    if (!e.target.files || e.target.files.length === 0) return;

    const selectedFiles = Array.from(e.target.files);
    setUploadError(null);

    // Initialize upload status
    const initialStatus: UploadingFile[] = selectedFiles.map((f) => ({
      name: f.name,
      status: 'uploading',
    }));
    setUploadingFiles(initialStatus);

    let successCount = 0;
    let errorCount = 0;

    // Upload files serially
    for (let i = 0; i < selectedFiles.length; i++) {
      const file = selectedFiles[i];
      try {
        const formData = new FormData();
        formData.append('file', file);
        formData.append('protocol', id);

        await apiUpload('protocol-documents/', formData);

        setUploadingFiles((prev) =>
          prev.map((f, idx) =>
            idx === i ? { ...f, status: 'success' } : f
          )
        );
        successCount++;
      } catch (err) {
        setUploadingFiles((prev) =>
          prev.map((f, idx) =>
            idx === i
              ? { ...f, status: 'error', error: err instanceof Error ? err.message : 'Upload failed' }
              : f
          )
        );
        errorCount++;
      }
    }

    // Refresh the documents list
    if (successCount > 0) {
      mutateDocuments();
    }

    // Clear upload status after a delay
    setTimeout(() => {
      setUploadingFiles([]);
      if (errorCount > 0) {
        setUploadError(`${errorCount} file(s) failed to upload`);
      }
    }, 3000);

    // Reset input
    e.target.value = '';
  }, [id, mutateDocuments]);

  const documentColumns: Column<ProtocolDocument>[] = [
    {
      key: 'filename',
      label: 'File',
      sortable: true,
      searchable: true,
      render: (value) => (
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Description fontSize="small" color="action" />
          <Typography fontWeight={500}>{value || 'Unnamed file'}</Typography>
        </Box>
      ),
    },
    {
      key: 'created_at',
      label: 'Uploaded',
      sortable: true,
      width: 120,
      render: (value) =>
        value ? new Date(value).toLocaleDateString() : '-',
    },
    {
      key: 'file',
      label: 'Actions',
      width: 140,
      render: (value, row) => (
        <Box sx={{ display: 'flex', gap: 0.5 }}>
          {value && (
            <>
              <Tooltip title="Download">
                <IconButton
                  size="small"
                  onClick={async (e) => {
                    e.stopPropagation();
                    const url = await getAuthenticatedDownloadUrl(value);
                    window.open(url, '_blank');
                  }}
                >
                  <Download fontSize="small" />
                </IconButton>
              </Tooltip>
              <Tooltip title="Open in new tab">
                <IconButton
                  size="small"
                  onClick={async (e) => {
                    e.stopPropagation();
                    const url = await getAuthenticatedDownloadUrl(value);
                    window.open(url, '_blank');
                  }}
                >
                  <OpenInNew fontSize="small" />
                </IconButton>
              </Tooltip>
            </>
          )}
          <Tooltip title={canContribute ? 'Delete' : 'Requires Contributor or Admin operating level'}>
            <span>
              <IconButton
                size="small"
                color="error"
                disabled={!canContribute}
                onClick={(e) => {
                  e.stopPropagation();
                  handleDeleteDocClick(row);
                }}
              >
                <Delete fontSize="small" />
              </IconButton>
            </span>
          </Tooltip>
        </Box>
      ),
    },
  ];

  // Helper to check if assay is new (created in last 7 days)
  const isRecentAssay = (createdAt: string | undefined) => {
    if (!createdAt) return false;
    const sevenDaysAgo = new Date();
    sevenDaysAgo.setDate(sevenDaysAgo.getDate() - 7);
    return new Date(createdAt) >= sevenDaysAgo;
  };

  const columns: Column<Assay>[] = [
    {
      key: 'data_filename',
      label: 'Data File',
      sortable: true,
      searchable: true,
      render: (value, row) => (
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Assessment fontSize="small" color="info" />
          <Typography fontWeight={500}>{value || 'No file'}</Typography>
          {isRecentAssay(row.created_at) && (
            <Tooltip title="New in the last 7 days">
              <FiberNew fontSize="small" color="secondary" />
            </Tooltip>
          )}
        </Box>
      ),
    },
    {
      key: 'target_name',
      label: 'Target',
      sortable: true,
      searchable: true,
      hiddenOnMobile: true,
      render: (value, row) =>
        value ? (
          <Chip
            icon={<Science fontSize="small" />}
            label={value}
            size="small"
            variant="outlined"
            onClick={(e) => {
              e.stopPropagation();
              if (row.target) router.push(routes.registry.target(row.target));
            }}
          />
        ) : (
          '-'
        ),
    },
    {
      key: 'data_series_count',
      label: 'Series',
      sortable: true,
      width: 80,
      hiddenOnMobile: true,
      render: (value) =>
        value !== undefined ? (
          <Chip label={value} size="small" color="secondary" variant="outlined" />
        ) : (
          '-'
        ),
    },
    {
      key: 'created_by_email',
      label: 'Created By',
      searchable: true,
      width: 180,
      hiddenOnMobile: true,
      render: (value) => value || '-',
    },
    {
      key: 'created_at',
      label: 'Date',
      sortable: true,
      width: 100,
      render: (value) =>
        value ? new Date(value).toLocaleDateString() : '-',
    },
    {
      key: 'actions',
      label: '',
      width: 50,
      hiddenOnMobile: true,
      render: (_, row) => (
        <Tooltip title={canContribute ? 'Delete' : 'Requires Contributor or Admin operating level'}>
          <span>
            <IconButton
              size="small"
              color="error"
              disabled={!canContribute}
              onClick={(e) => handleDeleteClick(row, e)}
              sx={{ opacity: 0.5, '&:hover': { opacity: 1 } }}
            >
              <Delete fontSize="small" />
            </IconButton>
          </span>
        </Tooltip>
      ),
    },
  ];

  // Build action button based on protocol import type
  const getAddAssayButton = () => {
    if (protocol?.import_type === 'table_of_values') {
      return (
        <Tooltip title={canContribute ? '' : 'Requires Contributor or Admin operating level'} arrow>
          <span>
            <Button
              component={Link}
              href={routes.assays.importTableOfValues({ protocol: id })}
              variant="contained"
              size="small"
              startIcon={<TableChart />}
              disabled={!canContribute}
            >
              Import
            </Button>
          </span>
        </Tooltip>
      );
    } else if (isAdmeProtocol(protocol?.import_type)) {
      return (
        <Tooltip title={canContribute ? '' : 'Requires Contributor or Admin operating level'} arrow>
          <span>
            <Button
              component={Link}
              href={routes.assays.importAdme()}
              variant="contained"
              size="small"
              startIcon={<Biotech />}
              disabled={!canContribute}
            >
              Import ADME
            </Button>
          </span>
        </Tooltip>
      );
    } else {
      return (
        <Tooltip title={canContribute ? '' : 'Requires Contributor or Admin operating level'} arrow>
          <span>
            <Button
              variant="contained"
              size="small"
              startIcon={<Add />}
              onClick={() => setUploadDrawerOpen(true)}
              disabled={!protocol || !canContribute}
            >
              Add Assay
            </Button>
          </span>
        </Tooltip>
      );
    }
  };

  // Summary configuration for collapsed header
  const summary = {
    title: protocol?.name || 'Loading...',
    titleIcon: <Description sx={{ fontSize: 'inherit' }} />,
    fields: [
      { label: 'Type', value: protocol ? IMPORT_TYPE_LABELS[protocol.import_type] : '-' },
      { label: 'Assays', value: assays ? `${assays.length}` : '-' },
      { label: 'Documents', value: documents ? `${documents.length}` : '-' },
    ],
    chips: protocol ? (
      <Chip
        label={IMPORT_TYPE_LABELS[protocol.import_type] || protocol.import_type}
        size="small"
        color="primary"
        variant="outlined"
      />
    ) : undefined,
    actions: (
      <>
        <Tooltip title={canContribute ? 'Edit protocol' : 'Requires Contributor or Admin operating level'} arrow>
          <span>
            <IconButton
              size="small"
              onClick={() => setEditDialogOpen(true)}
              disabled={!canContribute}
            >
              <Edit fontSize="small" />
            </IconButton>
          </span>
        </Tooltip>
        {getAddAssayButton()}
      </>
    ),
  };

  // Detail content: protocol info, plate layout, documents
  const detailContent = protocol ? (
    <Box>
      <Grid container spacing={3}>
        <Grid size={{ xs: 12, md: 6 }}>
          <Typography variant="h6" gutterBottom>
            Configuration
          </Typography>
          <InfoRow label="Data Type" value={IMPORT_TYPE_LABELS[protocol.import_type]} />
          {!isAdmeProtocol(protocol.import_type) && (
            <InfoRow
              label="Dilutions"
              value={
                protocol.preferred_dilutions_display ? (
                  <Typography component="span" sx={{ fontFamily: 'monospace', fontSize: '0.85rem' }}>
                    {protocol.preferred_dilutions_display}
                  </Typography>
                ) : null
              }
            />
          )}
          {isAdmeProtocol(protocol.import_type) && (
            <>
              <InfoRow label="Data Source" value="Pharmaron/NCU Excel exports" />
              <InfoRow
                label="File Pattern"
                value={
                  <Typography component="span" sx={{ fontFamily: 'monospace', fontSize: '0.85rem' }}>
                    ADME-NCU-*.xlsx
                  </Typography>
                }
              />
            </>
          )}
        </Grid>
        <Grid size={{ xs: 12, md: 6 }}>
          <Typography variant="h6" gutterBottom>
            Metadata
          </Typography>
          <InfoRow label="Created By" value={protocol.created_by_email} />
          <InfoRow
            label="Created"
            value={protocol.created_at ? new Date(protocol.created_at).toLocaleString() : null}
          />
          <InfoRow label="Assays" value={assays ? `${assays.length} experiments` : 'Loading...'} />
        </Grid>
      </Grid>

      {protocol.comments && (
        <>
          <Divider sx={{ my: 2 }} />
          <Typography variant="h6" gutterBottom>Comments</Typography>
          <Typography>{protocol.comments}</Typography>
        </>
      )}

      {protocol.fitting_method_name && (
        <>
          <Divider sx={{ my: 2 }} />
          <Typography variant="h6" gutterBottom>Fitting Method</Typography>
          <InfoRow
            label="Method"
            value={<Chip label={protocol.fitting_method_name} color="secondary" variant="outlined" size="small" />}
          />
          {protocol.fitting_parameters?.protein_conc && (
            <>
              <InfoRow label="[Protein]" value={`${protocol.fitting_parameters.protein_conc} nM`} />
              <InfoRow label="[Ligand]" value={`${protocol.fitting_parameters.ligand_conc} nM`} />
              <InfoRow label="Ligand Kd" value={`${protocol.fitting_parameters.ligand_kd} nM`} />
            </>
          )}
        </>
      )}

      {isPlateBasedProtocol(protocol.import_type) && (
        <>
          <Divider sx={{ my: 2 }} />
          <Accordion
            defaultExpanded={!(protocol.plate_layout_config && Object.keys(protocol.plate_layout_config).length > 0)}
            sx={{ boxShadow: 'none', '&:before': { display: 'none' }, bgcolor: 'transparent' }}
          >
            <AccordionSummary expandIcon={<ExpandMore />} sx={{ px: 0, minHeight: 'auto', '& .MuiAccordionSummary-content': { my: 1 } }}>
              <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, flex: 1 }}>
                <Typography variant="h6">Plate Layout</Typography>
                {protocol.plate_layout_name && (
                  <Chip icon={<GridOn />} label={protocol.plate_layout_name} color="primary" variant="outlined" size="small" />
                )}
                {!protocol.plate_layout_config || Object.keys(protocol.plate_layout_config).length === 0 ? (
                  <Chip label="Not configured" color="warning" variant="outlined" size="small" />
                ) : null}
              </Box>
            </AccordionSummary>
            <AccordionDetails sx={{ px: 0 }}>
              <Box sx={{ display: 'flex', gap: 1, mb: 2 }}>
                <Tooltip title={canContribute ? '' : 'Requires Contributor or Admin operating level'} arrow>
                  <span>
                    <Button size="small" startIcon={<Edit />} onClick={() => setLayoutSelectOpen(true)} disabled={!canContribute}>
                      {protocol.plate_layout ? 'Change Layout' : 'Select Layout'}
                    </Button>
                  </span>
                </Tooltip>
                {protocol.plate_layout && (
                  <Button size="small" component={Link} href={routes.assays.plateLayout(protocol.plate_layout)} startIcon={<OpenInNew />}>
                    View Details
                  </Button>
                )}
              </Box>
              {protocol.plate_layout_config && Object.keys(protocol.plate_layout_config).length > 0 ? (
                <Box sx={{ display: 'flex', justifyContent: 'center' }}>
                  <PlatePreview layout={protocol.plate_layout_config} width={450} height={300} />
                </Box>
              ) : (
                <Paper sx={{ p: 3, bgcolor: 'grey.50', textAlign: 'center' }}>
                  <GridOn sx={{ fontSize: 48, color: 'grey.400', mb: 1 }} />
                  <Typography color="text.secondary">No plate layout configured</Typography>
                  <Typography variant="body2" color="text.secondary" sx={{ mt: 0.5 }}>
                    Select a plate layout to define control positions, sample regions, and dilution patterns
                  </Typography>
                </Paper>
              )}
            </AccordionDetails>
          </Accordion>
        </>
      )}

      {/* Documents section */}
      <Accordion defaultExpanded={false} sx={{ boxShadow: 'none', '&:before': { display: 'none' }, bgcolor: 'transparent', mt: 2 }}>
        <AccordionSummary expandIcon={<ExpandMore />} sx={{ px: 0, minHeight: 'auto', '& .MuiAccordionSummary-content': { my: 1 } }}>
          <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
            <Description color="action" />
            <Typography variant="h6">Documents</Typography>
            {documents && documents.length > 0 && <Chip label={documents.length} size="small" variant="outlined" />}
          </Box>
        </AccordionSummary>
        <AccordionDetails sx={{ px: 0, pt: 0 }}>
          {uploadError && (
            <Alert severity="error" sx={{ mb: 2 }} onClose={() => setUploadError(null)}>{uploadError}</Alert>
          )}
          <Tooltip title={canContribute ? '' : 'Requires Contributor or Admin operating level'} arrow>
            <Paper
              variant="outlined"
              onDrop={canContribute ? handleDrop : undefined}
              onDragOver={canContribute ? handleDragOver : undefined}
              onDragLeave={canContribute ? handleDragLeave : undefined}
              sx={{
                p: 2, mb: 2, textAlign: 'center',
                bgcolor: isDragOver ? 'primary.50' : 'grey.50',
                borderStyle: 'dashed', borderColor: isDragOver ? 'primary.main' : 'grey.300', borderWidth: 2,
                cursor: canContribute ? 'pointer' : 'not-allowed',
                opacity: canContribute ? 1 : 0.6,
                transition: 'all 0.2s ease',
                '&:hover': canContribute ? { bgcolor: 'grey.100', borderColor: 'grey.400' } : {},
              }}
              onClick={canContribute ? () => document.getElementById('protocol-doc-input')?.click() : undefined}
            >
              <CloudUpload sx={{ fontSize: 32, color: isDragOver ? 'primary.main' : 'grey.400', mb: 0.5 }} />
              <Typography variant="body2" color={isDragOver ? 'primary.main' : 'text.secondary'}>
                {isDragOver ? 'Drop files here' : 'Drag and drop documents, or click to select'}
              </Typography>
              <input
                id="protocol-doc-input"
                type="file"
                multiple
                accept=".pdf,.doc,.docx,.xlsx,.xls,.txt,.csv,.png,.jpg,.jpeg"
                onChange={handleFileSelect}
                disabled={!canContribute}
                style={{ display: 'none' }}
              />
            </Paper>
          </Tooltip>
          {uploadingFiles.length > 0 && (
            <Paper sx={{ p: 2, mb: 2 }}>
              <Typography variant="subtitle2" sx={{ mb: 1 }}>Uploading {uploadingFiles.length} file(s)...</Typography>
              {uploadingFiles.map((file, index) => (
                <Box key={index} sx={{ display: 'flex', alignItems: 'center', gap: 1, mb: 0.5 }}>
                  {file.status === 'uploading' && <LinearProgress sx={{ flex: 1, height: 6, borderRadius: 1 }} />}
                  {file.status === 'success' && <CheckCircle color="success" fontSize="small" />}
                  {file.status === 'error' && <Description color="error" fontSize="small" />}
                  <Typography variant="body2" color={file.status === 'error' ? 'error' : 'text.secondary'} sx={{ minWidth: 200 }}>
                    {file.name}{file.status === 'error' && file.error && ` - ${file.error}`}
                  </Typography>
                </Box>
              ))}
            </Paper>
          )}
          <DataTable
            data={documents}
            columns={documentColumns}
            loading={documentsLoading}
            getRowKey={(row) => row.id}
            title={documents ? `${documents.length} document(s)` : undefined}
            emptyMessage="No documents uploaded for this protocol"
          />
        </AccordionDetails>
      </Accordion>
    </Box>
  ) : null;

  return (
    <>
      <DetailPageLayout
        breadcrumbs={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Protocols', href: routes.assays.protocols(), icon: 'protocol' },
          { label: protocol?.name || 'Loading...', icon: 'protocol' },
        ]}
        summary={summary}
        detailContent={detailContent}
        loading={loading}
      >
        <DataTable
          data={assays}
          columns={columns}
          loading={assaysLoading}
          onRowClick={(assay) => router.push(routes.assays.detail(assay.id))}
          getRowKey={(row) => row.id}
          title={assays ? `${assays.length} assays` : undefined}
          emptyMessage="No assays using this protocol"
          fillHeight
        />
      </DetailPageLayout>

      {/* Plate Layout Selection Dialog */}
      <Dialog open={layoutSelectOpen} onClose={() => setLayoutSelectOpen(false)} maxWidth="sm" fullWidth>
        <DialogTitle sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
          <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
            <GridOn color="primary" />
            Select Plate Layout
          </Box>
          <IconButton onClick={() => setLayoutSelectOpen(false)} size="small"><Close /></IconButton>
        </DialogTitle>
        <DialogContent dividers>
          {saveError && <Alert severity="error" sx={{ mb: 2 }}>{saveError}</Alert>}
          <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
            Select a plate layout for this protocol, or create a new one.
          </Typography>
          {plateLayouts && plateLayouts.length > 0 ? (
            <Box sx={{ display: 'flex', flexDirection: 'column', gap: 1 }}>
              {plateLayouts.map((layout) => (
                <Paper
                  key={layout.id}
                  variant="outlined"
                  sx={{
                    p: 2, cursor: 'pointer',
                    bgcolor: protocol?.plate_layout === layout.id ? 'primary.50' : 'background.paper',
                    borderColor: protocol?.plate_layout === layout.id ? 'primary.main' : 'divider',
                    '&:hover': { bgcolor: 'action.hover' },
                  }}
                  onClick={() => handleLayoutChange(layout.id)}
                >
                  <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
                    <Box>
                      <Typography fontWeight="medium">{layout.name}</Typography>
                      {layout.description && <Typography variant="body2" color="text.secondary">{layout.description}</Typography>}
                    </Box>
                    <Box sx={{ display: 'flex', gap: 1 }}>
                      {layout.plate_format && <Chip label={`${layout.plate_format}-well`} size="small" variant="outlined" />}
                      {protocol?.plate_layout === layout.id && <Chip label="Current" size="small" color="primary" />}
                    </Box>
                  </Box>
                </Paper>
              ))}
            </Box>
          ) : (
            <Typography color="text.secondary">No plate layouts available. Create one to get started.</Typography>
          )}
        </DialogContent>
        <DialogActions>
          {protocol?.plate_layout && (
            <Button onClick={() => handleLayoutChange(null)} color="error" disabled={saving}>Remove Layout</Button>
          )}
          <Box sx={{ flex: 1 }} />
          <Button onClick={() => setLayoutSelectOpen(false)}>Cancel</Button>
          <Button variant="contained" startIcon={<Add />} onClick={() => { setLayoutSelectOpen(false); setCreateLayoutOpen(true); }}>
            Create New Layout
          </Button>
        </DialogActions>
      </Dialog>

      <PlateLayoutCreateDialog open={createLayoutOpen} onClose={() => setCreateLayoutOpen(false)} onCreated={handleLayoutCreated} />

      {protocol && (
        <AssayUploadDrawer
          open={uploadDrawerOpen}
          onClose={() => setUploadDrawerOpen(false)}
          protocolId={id}
          protocolName={protocol.name}
          defaultTargetId={targetIdFromUrl || undefined}
          onAssayCreated={() => mutateAssays()}
        />
      )}

      <Dialog open={deleteDialogOpen} onClose={() => setDeleteDialogOpen(false)}>
        <DialogTitle>Delete Assay?</DialogTitle>
        <DialogContent>
          <DialogContentText>
            This will permanently delete &quot;{assayToDelete?.data_filename}&quot; and all its data series with analysis results.
          </DialogContentText>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setDeleteDialogOpen(false)} disabled={deleting}>Cancel</Button>
          <Button onClick={handleDeleteConfirm} color="error" variant="contained" disabled={deleting} startIcon={deleting ? <CircularProgress size={16} /> : <Delete />}>
            {deleting ? 'Deleting...' : 'Delete'}
          </Button>
        </DialogActions>
      </Dialog>

      <Dialog open={deleteDocDialogOpen} onClose={() => setDeleteDocDialogOpen(false)}>
        <DialogTitle>Delete Document</DialogTitle>
        <DialogContent>
          <DialogContentText>
            Are you sure you want to delete &quot;{docToDelete?.filename || 'this document'}&quot;? This action cannot be undone.
          </DialogContentText>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setDeleteDocDialogOpen(false)} disabled={deletingDoc}>Cancel</Button>
          <Button onClick={handleDeleteDocConfirm} color="error" variant="contained" disabled={deletingDoc} startIcon={deletingDoc ? <CircularProgress size={16} /> : <Delete />}>
            {deletingDoc ? 'Deleting...' : 'Delete'}
          </Button>
        </DialogActions>
      </Dialog>

      {protocol && (
        <ProtocolEditDialog open={editDialogOpen} onClose={() => setEditDialogOpen(false)} protocol={protocol} onSave={() => mutate()} />
      )}
    </>
  );
}

function ProtocolDetailPageFallback() {
  return (
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <Skeleton variant="rectangular" height={40} sx={{ mb: 2 }} />
      <Skeleton variant="rectangular" height={200} sx={{ mb: 3 }} />
      <Skeleton variant="rectangular" height={400} />
    </Container>
  );
}

export default function ProtocolDetailPage({ params }: PageProps) {
  return (
    <Suspense fallback={<ProtocolDetailPageFallback />}>
      <ProtocolDetailPageContent params={params} />
    </Suspense>
  );
}
