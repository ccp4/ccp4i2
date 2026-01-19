'use client';

import { use, useState, useCallback } from 'react';
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
  IconButton,
  Tooltip,
  LinearProgress,
  Alert,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogContentText,
  DialogActions,
  Button,
  CircularProgress,
} from '@mui/material';
import {
  Inventory,
  Medication,
  Science,
  Description,
  Download,
  OpenInNew,
  CloudUpload,
  CheckCircle,
  Delete,
} from '@mui/icons-material';
import { PageHeader } from '@/components/compounds/PageHeader';
import { DataTable, Column } from '@/components/compounds/DataTable';
import { useCompoundsApi, apiUpload } from '@/lib/compounds/api';
import { routes } from '@/lib/compounds/routes';
import { Batch, BatchQCFile, Compound, Target } from '@/types/compounds/models';

interface UploadingFile {
  name: string;
  status: 'uploading' | 'success' | 'error';
  error?: string;
}

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

export default function BatchDetailPage({ params }: PageProps) {
  const { id } = use(params);
  const router = useRouter();
  const api = useCompoundsApi();

  // Drop zone state
  const [isDragOver, setIsDragOver] = useState(false);
  const [uploadingFiles, setUploadingFiles] = useState<UploadingFile[]>([]);
  const [uploadError, setUploadError] = useState<string | null>(null);

  // Delete state
  const [deleteDialogOpen, setDeleteDialogOpen] = useState(false);
  const [fileToDelete, setFileToDelete] = useState<BatchQCFile | null>(null);
  const [deleting, setDeleting] = useState(false);

  const { data: batch, isLoading: batchLoading } = api.get<Batch>(
    `batches/${id}/`
  );
  const { data: qcFiles, isLoading: qcFilesLoading, mutate: mutateQcFiles } = api.get<BatchQCFile[]>(
    `batch-qc-files/?batch=${id}`
  );
  const { data: compound } = api.get<Compound>(
    batch?.compound ? `compounds/${batch.compound}/` : null
  );
  const { data: target } = api.get<Target>(
    compound?.target ? `targets/${compound.target}/` : null
  );

  // Handle drag and drop for QC files
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
        formData.append('batch', id);

        await apiUpload('batch-qc-files/', formData);

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

    // Refresh the QC files list
    if (successCount > 0) {
      mutateQcFiles();
    }

    // Clear upload status after a delay
    setTimeout(() => {
      setUploadingFiles([]);
      if (errorCount > 0) {
        setUploadError(`${errorCount} file(s) failed to upload`);
      }
    }, 3000);
  }, [id, mutateQcFiles]);

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
        formData.append('batch', id);

        await apiUpload('batch-qc-files/', formData);

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

    // Refresh the QC files list
    if (successCount > 0) {
      mutateQcFiles();
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
  }, [id, mutateQcFiles]);

  const handleDeleteClick = (file: BatchQCFile) => {
    setFileToDelete(file);
    setDeleteDialogOpen(true);
  };

  const handleDeleteConfirm = async () => {
    if (!fileToDelete) return;

    setDeleting(true);
    try {
      await api.delete(`batch-qc-files/${fileToDelete.id}/`);
      mutateQcFiles();
    } catch (err) {
      console.error('Delete error:', err);
      setUploadError('Failed to delete file');
    } finally {
      setDeleting(false);
      setDeleteDialogOpen(false);
      setFileToDelete(null);
    }
  };

  const columns: Column<BatchQCFile>[] = [
    {
      key: 'filename',
      label: 'File',
      sortable: true,
      searchable: true,
      render: (value, row) => (
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <Description fontSize="small" color="action" />
          <Typography fontWeight={500}>{value || 'Unnamed file'}</Typography>
        </Box>
      ),
    },
    {
      key: 'comments',
      label: 'Comments',
      searchable: true,
      render: (value) =>
        value ? (
          <Typography
            sx={{
              maxWidth: 300,
              overflow: 'hidden',
              textOverflow: 'ellipsis',
              whiteSpace: 'nowrap',
            }}
            title={value}
          >
            {value}
          </Typography>
        ) : (
          '-'
        ),
    },
    {
      key: 'uploaded_at',
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
                  href={value}
                  target="_blank"
                  onClick={(e) => e.stopPropagation()}
                >
                  <Download fontSize="small" />
                </IconButton>
              </Tooltip>
              <Tooltip title="Open in new tab">
                <IconButton
                  size="small"
                  href={value}
                  target="_blank"
                  onClick={(e) => e.stopPropagation()}
                >
                  <OpenInNew fontSize="small" />
                </IconButton>
              </Tooltip>
            </>
          )}
          <Tooltip title="Delete">
            <IconButton
              size="small"
              color="error"
              onClick={(e) => {
                e.stopPropagation();
                handleDeleteClick(row);
              }}
            >
              <Delete fontSize="small" />
            </IconButton>
          </Tooltip>
        </Box>
      ),
    },
  ];

  return (
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <PageHeader
        breadcrumbs={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Targets', href: routes.registry.targets(), icon: 'target' },
          {
            label: target?.name || 'Target',
            href: compound?.target
              ? routes.registry.target(compound.target)
              : undefined,
            icon: 'target',
          },
          {
            label: compound?.formatted_id || 'Compound',
            href: batch?.compound
              ? routes.registry.compound(batch.compound)
              : undefined,
            icon: 'compound',
          },
          {
            label: batch ? `Batch #${batch.batch_number}` : 'Loading...',
            icon: 'batch',
          },
        ]}
      />

      {/* Batch header */}
      <Paper sx={{ p: 3, mb: 3 }}>
        {batchLoading ? (
          <>
            <Skeleton variant="text" width={300} height={40} />
            <Skeleton variant="text" width={200} />
          </>
        ) : batch ? (
          <>
            <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 2 }}>
              <Inventory sx={{ fontSize: 48, color: 'info.main' }} />
              <Box>
                <Typography variant="h4">Batch #{batch.batch_number}</Typography>
                <Box sx={{ display: 'flex', gap: 1, mt: 0.5 }}>
                  {compound && (
                    <Chip
                      icon={<Medication fontSize="small" />}
                      label={compound.formatted_id}
                      size="small"
                      onClick={() =>
                        router.push(routes.registry.compound(compound.id))
                      }
                    />
                  )}
                  {target && (
                    <Chip
                      icon={<Science fontSize="small" />}
                      label={target.name}
                      size="small"
                      variant="outlined"
                      onClick={() =>
                        router.push(routes.registry.target(target.id))
                      }
                    />
                  )}
                </Box>
              </Box>
            </Box>

            <Divider sx={{ my: 2 }} />

            <Grid container spacing={3}>
              <Grid size={{ xs: 12, md: 6 }}>
                <Typography variant="h6" gutterBottom>
                  Properties
                </Typography>
                <InfoRow
                  label="Amount"
                  value={
                    batch.amount ? `${parseFloat(batch.amount).toFixed(2)} mg` : null
                  }
                />
                <InfoRow label="Salt Code" value={batch.salt_code} />
                <InfoRow
                  label="MW (salt)"
                  value={batch.molecular_weight?.toFixed(2)}
                />
              </Grid>
              <Grid size={{ xs: 12, md: 6 }}>
                <Typography variant="h6" gutterBottom>
                  Provenance
                </Typography>
                <InfoRow label="Supplier" value={batch.supplier_name} />
                <InfoRow label="Supplier Ref" value={batch.supplier_ref} />
                <InfoRow label="Lab Book" value={batch.labbook_number} />
                <InfoRow label="Page" value={batch.page_number} />
                <InfoRow
                  label="Registered"
                  value={
                    batch.registered_at
                      ? new Date(batch.registered_at).toLocaleString()
                      : null
                  }
                />
              </Grid>
            </Grid>

            {batch.comments && (
              <>
                <Divider sx={{ my: 2 }} />
                <Typography variant="h6" gutterBottom>
                  Comments
                </Typography>
                <Typography>{batch.comments}</Typography>
              </>
            )}
          </>
        ) : (
          <Typography color="error">Batch not found</Typography>
        )}
      </Paper>

      {/* QC Files section */}
      <Typography variant="h6" sx={{ mb: 2 }}>
        QC Files
      </Typography>

      {/* Upload error alert */}
      {uploadError && (
        <Alert severity="error" sx={{ mb: 2 }} onClose={() => setUploadError(null)}>
          {uploadError}
        </Alert>
      )}

      {/* Drop zone for QC file uploads */}
      <Paper
        variant="outlined"
        onDrop={handleDrop}
        onDragOver={handleDragOver}
        onDragLeave={handleDragLeave}
        sx={{
          p: 3,
          mb: 2,
          textAlign: 'center',
          bgcolor: isDragOver ? 'primary.50' : 'grey.50',
          borderStyle: 'dashed',
          borderColor: isDragOver ? 'primary.main' : 'grey.300',
          borderWidth: 2,
          cursor: 'pointer',
          transition: 'all 0.2s ease',
          '&:hover': { bgcolor: 'grey.100', borderColor: 'grey.400' },
        }}
        onClick={() => document.getElementById('qc-file-input')?.click()}
      >
        <CloudUpload sx={{ fontSize: 40, color: isDragOver ? 'primary.main' : 'grey.400', mb: 1 }} />
        <Typography color={isDragOver ? 'primary.main' : 'text.secondary'}>
          {isDragOver ? 'Drop files here' : 'Drag and drop QC files here, or click to select'}
        </Typography>
        <Typography variant="body2" color="text.secondary" sx={{ mt: 0.5 }}>
          Multiple files supported (PDF, images, Excel, etc.)
        </Typography>
        <input
          id="qc-file-input"
          type="file"
          multiple
          accept=".pdf,.doc,.docx,.xlsx,.xls,.txt,.csv,.png,.jpg,.jpeg"
          onChange={handleFileSelect}
          style={{ display: 'none' }}
        />
      </Paper>

      {/* Upload progress */}
      {uploadingFiles.length > 0 && (
        <Paper sx={{ p: 2, mb: 2 }}>
          <Typography variant="subtitle2" sx={{ mb: 1 }}>
            Uploading {uploadingFiles.length} file(s)...
          </Typography>
          {uploadingFiles.map((file, index) => (
            <Box key={index} sx={{ display: 'flex', alignItems: 'center', gap: 1, mb: 0.5 }}>
              {file.status === 'uploading' && <LinearProgress sx={{ flex: 1, height: 6, borderRadius: 1 }} />}
              {file.status === 'success' && <CheckCircle color="success" fontSize="small" />}
              {file.status === 'error' && <Description color="error" fontSize="small" />}
              <Typography
                variant="body2"
                color={file.status === 'error' ? 'error' : 'text.secondary'}
                sx={{ minWidth: 200 }}
              >
                {file.name}
                {file.status === 'error' && file.error && ` - ${file.error}`}
              </Typography>
            </Box>
          ))}
        </Paper>
      )}

      {/* QC Files table */}
      <DataTable
        data={qcFiles}
        columns={columns}
        loading={qcFilesLoading}
        getRowKey={(row) => row.id}
        title={qcFiles ? `${qcFiles.length} QC file(s)` : undefined}
        emptyMessage="No QC files uploaded for this batch"
      />

      {/* Delete confirmation dialog */}
      <Dialog
        open={deleteDialogOpen}
        onClose={() => setDeleteDialogOpen(false)}
      >
        <DialogTitle>Delete QC File</DialogTitle>
        <DialogContent>
          <DialogContentText>
            Are you sure you want to delete &quot;{fileToDelete?.filename || 'this file'}&quot;?
            This action cannot be undone.
          </DialogContentText>
        </DialogContent>
        <DialogActions>
          <Button
            onClick={() => setDeleteDialogOpen(false)}
            disabled={deleting}
          >
            Cancel
          </Button>
          <Button
            onClick={handleDeleteConfirm}
            color="error"
            variant="contained"
            disabled={deleting}
            startIcon={deleting ? <CircularProgress size={16} /> : <Delete />}
          >
            {deleting ? 'Deleting...' : 'Delete'}
          </Button>
        </DialogActions>
      </Dialog>
    </Container>
  );
}
