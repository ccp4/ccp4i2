'use client';

import { useState, useCallback } from 'react';
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  Box,
  Typography,
  Alert,
  IconButton,
  CircularProgress,
  List,
  ListItem,
  ListItemIcon,
  ListItemText,
  ListItemSecondaryAction,
  Paper,
  LinearProgress,
} from '@mui/material';
import {
  Close,
  CloudUpload,
  Description,
  Delete,
  CheckCircle,
  Error as ErrorIcon,
} from '@mui/icons-material';

interface UploadFile {
  file: File;
  status: 'pending' | 'uploading' | 'success' | 'error';
  error?: string;
  progress?: number;
}

interface DocumentUploadDialogProps {
  open: boolean;
  onClose: () => void;
  /** Title for the dialog */
  title: string;
  /** API endpoint to POST files to (e.g., 'batch-qc-files/' or 'protocol-documents/') */
  endpoint: string;
  /** Field name for the parent reference (e.g., 'batch' or 'protocol') */
  parentField: string;
  /** ID of the parent object */
  parentId: string;
  /** Callback after successful upload(s) */
  onUploaded: () => void;
  /** Accepted file types (e.g., '.pdf,.doc,.docx,.xlsx') */
  accept?: string;
}

export function DocumentUploadDialog({
  open,
  onClose,
  title,
  endpoint,
  parentField,
  parentId,
  onUploaded,
  accept = '.pdf,.doc,.docx,.xlsx,.xls,.txt,.csv,.png,.jpg,.jpeg',
}: DocumentUploadDialogProps) {
  const [files, setFiles] = useState<UploadFile[]>([]);
  const [uploading, setUploading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleDrop = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    const droppedFiles = Array.from(e.dataTransfer.files);
    addFiles(droppedFiles);
  }, []);

  const handleDragOver = useCallback((e: React.DragEvent) => {
    e.preventDefault();
  }, []);

  const handleFileSelect = (e: React.ChangeEvent<HTMLInputElement>) => {
    if (e.target.files) {
      addFiles(Array.from(e.target.files));
    }
  };

  const addFiles = (newFiles: File[]) => {
    const uploadFiles: UploadFile[] = newFiles.map((file) => ({
      file,
      status: 'pending',
    }));
    setFiles((prev) => [...prev, ...uploadFiles]);
    setError(null);
  };

  const removeFile = (index: number) => {
    setFiles((prev) => prev.filter((_, i) => i !== index));
  };

  const handleUpload = async () => {
    if (files.length === 0) return;

    setUploading(true);
    setError(null);

    let successCount = 0;
    let errorCount = 0;

    for (let i = 0; i < files.length; i++) {
      const uploadFile = files[i];
      if (uploadFile.status === 'success') {
        successCount++;
        continue;
      }

      // Update status to uploading
      setFiles((prev) =>
        prev.map((f, idx) =>
          idx === i ? { ...f, status: 'uploading' as const, progress: 0 } : f
        )
      );

      try {
        const formData = new FormData();
        formData.append('file', uploadFile.file);
        formData.append(parentField, parentId);

        const response = await fetch(`/api/proxy/compounds/${endpoint}`, {
          method: 'POST',
          body: formData,
        });

        if (!response.ok) {
          const errorData = await response.json().catch(() => ({}));
          const errorMsg = errorData.detail || errorData.file?.[0] || 'Upload failed';
          throw new Error(errorMsg);
        }

        setFiles((prev) =>
          prev.map((f, idx) =>
            idx === i ? { ...f, status: 'success' as const } : f
          )
        );
        successCount++;
      } catch (err) {
        setFiles((prev) =>
          prev.map((f, idx) =>
            idx === i
              ? {
                  ...f,
                  status: 'error' as const,
                  error: err instanceof Error ? err.message : 'Upload failed',
                }
              : f
          )
        );
        errorCount++;
      }
    }

    setUploading(false);

    if (successCount > 0) {
      onUploaded();
    }

    if (errorCount === 0 && successCount > 0) {
      // All successful - close dialog
      handleClose();
    } else if (errorCount > 0) {
      setError(`${errorCount} file(s) failed to upload`);
    }
  };

  const handleClose = () => {
    if (!uploading) {
      setFiles([]);
      setError(null);
      onClose();
    }
  };

  const pendingFiles = files.filter((f) => f.status === 'pending' || f.status === 'error');
  const canUpload = pendingFiles.length > 0 && !uploading;

  return (
    <Dialog open={open} onClose={handleClose} maxWidth="sm" fullWidth>
      <DialogTitle sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <CloudUpload color="primary" />
          {title}
        </Box>
        <IconButton onClick={handleClose} size="small" disabled={uploading}>
          <Close />
        </IconButton>
      </DialogTitle>

      <DialogContent dividers>
        {error && (
          <Alert severity="error" sx={{ mb: 2 }}>
            {error}
          </Alert>
        )}

        {/* Drop zone */}
        <Paper
          variant="outlined"
          onDrop={handleDrop}
          onDragOver={handleDragOver}
          sx={{
            p: 4,
            textAlign: 'center',
            bgcolor: 'grey.50',
            borderStyle: 'dashed',
            cursor: 'pointer',
            '&:hover': { bgcolor: 'grey.100' },
            mb: 2,
          }}
          onClick={() => document.getElementById('file-input')?.click()}
        >
          <CloudUpload sx={{ fontSize: 48, color: 'grey.400', mb: 1 }} />
          <Typography color="text.secondary">
            Drag and drop files here, or click to select
          </Typography>
          <Typography variant="body2" color="text.secondary" sx={{ mt: 0.5 }}>
            Multiple files supported
          </Typography>
          <input
            id="file-input"
            type="file"
            multiple
            accept={accept}
            onChange={handleFileSelect}
            style={{ display: 'none' }}
          />
        </Paper>

        {/* File list */}
        {files.length > 0 && (
          <List dense>
            {files.map((uploadFile, index) => (
              <ListItem
                key={index}
                sx={{
                  bgcolor:
                    uploadFile.status === 'error'
                      ? 'error.50'
                      : uploadFile.status === 'success'
                      ? 'success.50'
                      : 'transparent',
                  borderRadius: 1,
                  mb: 0.5,
                }}
              >
                <ListItemIcon>
                  {uploadFile.status === 'success' ? (
                    <CheckCircle color="success" />
                  ) : uploadFile.status === 'error' ? (
                    <ErrorIcon color="error" />
                  ) : uploadFile.status === 'uploading' ? (
                    <CircularProgress size={20} />
                  ) : (
                    <Description color="action" />
                  )}
                </ListItemIcon>
                <ListItemText
                  primary={uploadFile.file.name}
                  secondary={
                    uploadFile.status === 'error'
                      ? uploadFile.error
                      : uploadFile.status === 'uploading'
                      ? 'Uploading...'
                      : uploadFile.status === 'success'
                      ? 'Uploaded'
                      : `${(uploadFile.file.size / 1024).toFixed(1)} KB`
                  }
                  secondaryTypographyProps={{
                    color: uploadFile.status === 'error' ? 'error' : 'text.secondary',
                  }}
                />
                {uploadFile.status === 'pending' && !uploading && (
                  <ListItemSecondaryAction>
                    <IconButton
                      edge="end"
                      size="small"
                      onClick={() => removeFile(index)}
                    >
                      <Delete fontSize="small" />
                    </IconButton>
                  </ListItemSecondaryAction>
                )}
              </ListItem>
            ))}
          </List>
        )}

        {uploading && (
          <Box sx={{ mt: 2 }}>
            <LinearProgress />
          </Box>
        )}
      </DialogContent>

      <DialogActions>
        <Button onClick={handleClose} disabled={uploading}>
          {files.some((f) => f.status === 'success') ? 'Done' : 'Cancel'}
        </Button>
        <Button
          variant="contained"
          onClick={handleUpload}
          disabled={!canUpload}
          startIcon={uploading ? <CircularProgress size={16} /> : <CloudUpload />}
        >
          {uploading ? 'Uploading...' : `Upload ${pendingFiles.length} file(s)`}
        </Button>
      </DialogActions>
    </Dialog>
  );
}
