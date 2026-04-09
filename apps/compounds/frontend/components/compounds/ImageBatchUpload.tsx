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
  Paper,
  LinearProgress,
  Chip,
  Divider,
} from '@mui/material';
import {
  Close,
  CloudUpload,
  Image as ImageIcon,
  CheckCircle,
  Error as ErrorIcon,
  Warning,
  Delete,
} from '@mui/icons-material';
import { apiUpload } from '@/lib/compounds/api';

interface UploadFile {
  file: File;
  status: 'pending' | 'uploading' | 'matched' | 'unmatched';
  matchedTo?: string;  // compound_name if matched
}

interface ImageBatchUploadProps {
  open: boolean;
  onClose: () => void;
  /** Assay ID to upload images to */
  assayId: string;
  /** Callback after upload completes */
  onUploaded: () => void;
}

interface UploadResponse {
  status: string;
  matched: number;
  unmatched: number;
  matched_files: Array<{
    filename: string;
    data_series_id: string;
    compound_name: string;
  }>;
  unmatched_files: string[];
}

/**
 * Component for batch uploading plot images to a Table-Of-Values assay.
 *
 * Images are matched to data series by comparing the uploaded filename
 * to the 'Image File' field in each DataSeries analysis results.
 */
export function ImageBatchUpload({
  open,
  onClose,
  assayId,
  onUploaded,
}: ImageBatchUploadProps) {
  const [files, setFiles] = useState<UploadFile[]>([]);
  const [uploading, setUploading] = useState(false);
  const [uploadResult, setUploadResult] = useState<UploadResponse | null>(null);
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
    // Filter to only image files
    const imageFiles = newFiles.filter((f) =>
      f.type.startsWith('image/') ||
      /\.(png|jpg|jpeg|gif|svg|webp)$/i.test(f.name)
    );

    const uploadFiles: UploadFile[] = imageFiles.map((file) => ({
      file,
      status: 'pending',
    }));
    setFiles((prev) => [...prev, ...uploadFiles]);
    setError(null);
    setUploadResult(null);
  };

  const removeFile = (index: number) => {
    setFiles((prev) => prev.filter((_, i) => i !== index));
  };

  const handleUpload = async () => {
    if (files.length === 0) return;

    setUploading(true);
    setError(null);
    setUploadResult(null);

    // Mark all files as uploading
    setFiles((prev) =>
      prev.map((f) => ({ ...f, status: 'uploading' as const }))
    );

    try {
      const formData = new FormData();
      files.forEach((uploadFile) => {
        formData.append('files', uploadFile.file);
      });

      const response = await apiUpload<UploadResponse>(
        `assays/${assayId}/upload_images/`,
        formData
      );

      setUploadResult(response);

      // Update file statuses based on response
      setFiles((prev) =>
        prev.map((f) => {
          const matchedFile = response.matched_files.find(
            (m) => m.filename === f.file.name
          );
          if (matchedFile) {
            return {
              ...f,
              status: 'matched' as const,
              matchedTo: matchedFile.compound_name,
            };
          }
          return { ...f, status: 'unmatched' as const };
        })
      );

      if (response.matched > 0) {
        onUploaded();
      }
    } catch (err) {
      console.error('Upload failed:', err);
      setError(err instanceof Error ? err.message : 'Upload failed');
      // Reset status on error
      setFiles((prev) =>
        prev.map((f) => ({ ...f, status: 'pending' as const }))
      );
    } finally {
      setUploading(false);
    }
  };

  const handleClose = () => {
    if (!uploading) {
      setFiles([]);
      setError(null);
      setUploadResult(null);
      onClose();
    }
  };

  const pendingFiles = files.filter((f) => f.status === 'pending');
  const canUpload = pendingFiles.length > 0 && !uploading;
  const hasResults = uploadResult !== null;

  return (
    <Dialog open={open} onClose={handleClose} maxWidth="sm" fullWidth>
      <DialogTitle sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <ImageIcon color="primary" />
          Upload Plot Images
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

        {/* Upload result summary */}
        {uploadResult && (
          <Alert
            severity={uploadResult.unmatched > 0 ? 'warning' : 'success'}
            sx={{ mb: 2 }}
          >
            <Box sx={{ display: 'flex', gap: 1, alignItems: 'center', flexWrap: 'wrap' }}>
              <Chip
                icon={<CheckCircle />}
                label={`${uploadResult.matched} matched`}
                color="success"
                size="small"
              />
              {uploadResult.unmatched > 0 && (
                <Chip
                  icon={<Warning />}
                  label={`${uploadResult.unmatched} unmatched`}
                  color="warning"
                  size="small"
                />
              )}
            </Box>
            {uploadResult.unmatched > 0 && (
              <Typography variant="body2" sx={{ mt: 1 }}>
                Unmatched files were not saved. Check that filenames match the "Image File" column values.
              </Typography>
            )}
          </Alert>
        )}

        {/* Drop zone - only show if no results yet */}
        {!hasResults && (
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
            onClick={() => document.getElementById('image-file-input')?.click()}
          >
            <CloudUpload sx={{ fontSize: 48, color: 'grey.400', mb: 1 }} />
            <Typography color="text.secondary">
              Drag and drop plot images here, or click to select
            </Typography>
            <Typography variant="body2" color="text.secondary" sx={{ mt: 0.5 }}>
              PNG, JPG, SVG, GIF supported. Files will be matched by filename.
            </Typography>
            <input
              id="image-file-input"
              type="file"
              multiple
              accept="image/*,.svg"
              onChange={handleFileSelect}
              style={{ display: 'none' }}
            />
          </Paper>
        )}

        {/* File list */}
        {files.length > 0 && (
          <>
            {!hasResults && (
              <Typography variant="subtitle2" gutterBottom>
                {files.length} file(s) selected
              </Typography>
            )}
            <List dense sx={{ maxHeight: 300, overflow: 'auto' }}>
              {files.map((uploadFile, index) => (
                <ListItem
                  key={index}
                  sx={{
                    bgcolor:
                      uploadFile.status === 'unmatched'
                        ? 'warning.50'
                        : uploadFile.status === 'matched'
                        ? 'success.50'
                        : 'transparent',
                    borderRadius: 1,
                    mb: 0.5,
                  }}
                >
                  <ListItemIcon>
                    {uploadFile.status === 'matched' ? (
                      <CheckCircle color="success" />
                    ) : uploadFile.status === 'unmatched' ? (
                      <Warning color="warning" />
                    ) : uploadFile.status === 'uploading' ? (
                      <CircularProgress size={20} />
                    ) : (
                      <ImageIcon color="action" />
                    )}
                  </ListItemIcon>
                  <ListItemText
                    primary={uploadFile.file.name}
                    secondary={
                      uploadFile.status === 'matched'
                        ? `Matched to ${uploadFile.matchedTo}`
                        : uploadFile.status === 'unmatched'
                        ? 'No matching data series found'
                        : uploadFile.status === 'uploading'
                        ? 'Uploading...'
                        : `${(uploadFile.file.size / 1024).toFixed(1)} KB`
                    }
                    secondaryTypographyProps={{
                      color:
                        uploadFile.status === 'unmatched'
                          ? 'warning.main'
                          : uploadFile.status === 'matched'
                          ? 'success.main'
                          : 'text.secondary',
                    }}
                  />
                  {uploadFile.status === 'pending' && !uploading && (
                    <IconButton
                      edge="end"
                      size="small"
                      onClick={() => removeFile(index)}
                    >
                      <Delete fontSize="small" />
                    </IconButton>
                  )}
                </ListItem>
              ))}
            </List>
          </>
        )}

        {uploading && (
          <Box sx={{ mt: 2 }}>
            <LinearProgress />
          </Box>
        )}
      </DialogContent>

      <DialogActions>
        <Button onClick={handleClose} disabled={uploading}>
          {hasResults ? 'Done' : 'Cancel'}
        </Button>
        {!hasResults && (
          <Button
            variant="contained"
            onClick={handleUpload}
            disabled={!canUpload}
            startIcon={uploading ? <CircularProgress size={16} /> : <CloudUpload />}
          >
            {uploading ? 'Uploading...' : `Upload ${files.length} image(s)`}
          </Button>
        )}
      </DialogActions>
    </Dialog>
  );
}
