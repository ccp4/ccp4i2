import React, { useEffect, useState } from "react";
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  List,
  ListItem,
  ListItemButton,
  ListItemText,
  CircularProgress,
  Alert,
  Typography,
  Box,
  Divider,
  Chip,
} from "@mui/material";
import {
  Download,
  Close,
  InsertDriveFile,
  Description,
  Image,
  DataObject,
} from "@mui/icons-material";
import { apiFetch, apiGet } from "../api-fetch";

/**
 * Props interface for the ExportJobMenu component
 */
interface ExportJobMenuProps {
  /** The ID of the job to export, null/undefined to close dialog */
  jobId: number | null;
  /** Callback to set the job ID, call with null to close dialog */
  setJobId: (jobId: number | null) => void;
}

/**
 * File Menu Item structure - array of three strings:
 * [0] identifier (descriptive string)
 * [1] menu label (display name)
 * [2] mimetype (file type)
 */
type FileMenuItem = [string, string, string];

/**
 * Parsed file menu item for easier handling
 */
interface ParsedFileMenuItem {
  identifier: string;
  label: string;
  mimetype: string;
}

/**
 * API response structure for export job file menu
 * Supports both legacy format (status/result/reason) and new format (success/data/error)
 */
interface ExportJobFileMenuResponse {
  // Legacy format
  status?: "Success" | "Failed";
  result?: FileMenuItem[] | any; // Could be array of FileMenuItem or other structure
  reason?: string;
  // New standardized format
  success?: boolean;
  data?: {
    result?: FileMenuItem[] | any;
    [key: string]: any;
  };
  error?: string;
}

/**
 * Download utility function
 *
 * Handles file downloads with proper error handling and user feedback
 * optimized for cloud-based deployments.
 *
 * @param url - The URL to download from
 * @param filename - Suggested filename for the download
 */
const doDownload = async (url: string, filename: string): Promise<void> => {
  try {
    // Fetch the response to get headers
    const response = await apiFetch(url);

    // apiFetch already checks response.ok, so we don't need to check again

    // Try to get filename from X-Export-Info header
    let downloadFilename = filename; // fallback
    const exportInfoHeader = response.headers.get("X-Export-Info");
    if (exportInfoHeader) {
      try {
        const exportInfo = JSON.parse(exportInfoHeader);
        if (exportInfo.original_filename) {
          downloadFilename = exportInfo.original_filename;
        }
      } catch (parseError) {
        console.warn("Failed to parse X-Export-Info header:", parseError);
        // Continue with fallback filename
      }
    }

    // Get the blob data
    const blob = await response.blob();

    // Create object URL for the blob
    const blobUrl = window.URL.createObjectURL(blob);

    // Create a temporary anchor element for download
    const link = document.createElement("a");
    link.href = blobUrl;
    link.download = downloadFilename;

    // Add to DOM temporarily
    document.body.appendChild(link);

    // Trigger download
    link.click();

    // Clean up
    document.body.removeChild(link);
    window.URL.revokeObjectURL(blobUrl);

    console.log(`Download initiated for: ${downloadFilename}`);
  } catch (error) {
    console.error("Download failed:", error);
    throw new Error(
      `Download failed for ${filename}: ${
        error instanceof Error ? error.message : "Unknown error"
      }`
    );
  }
};

/**
 * ExportJobMenu Component
 *
 * A dialog component that displays exportable file menu items for a specific job.
 * Fetches available export options from the Django backend and presents them
 * in a user-friendly dialog interface.
 *
 * - Handles API calls to containerized Django backend
 * - Optimized for cloud-based file operations
 * - Includes proper error handling for distributed environments
 *
 * @param props - Component properties
 * @returns React functional component
 */
export const ExportJobMenu: React.FC<ExportJobMenuProps> = ({
  jobId,
  setJobId,
}) => {
  // State for managing file menu items and loading state
  const [fileMenuItems, setFileMenuItems] = useState<ParsedFileMenuItem[]>([]);
  const [loading, setLoading] = useState<boolean>(false);
  const [error, setError] = useState<string | null>(null);
  const [downloadingItems, setDownloadingItems] = useState<Set<string>>(
    new Set()
  );

  // Determine if dialog should be open
  const isOpen = jobId !== null && jobId !== undefined;

  /**
   * Parse raw file menu item array into structured object
   */
  const parseFileMenuItem = (item: FileMenuItem): ParsedFileMenuItem => {
    return {
      identifier: item[0] || "",
      label: item[1] || "",
      mimetype: item[2] || "",
    };
  };

  /**
   * Validate if an item is a valid FileMenuItem (array of 3 strings)
   */
  const isValidFileMenuItem = (item: any): item is FileMenuItem => {
    return (
      Array.isArray(item) &&
      item.length === 3 &&
      typeof item[0] === "string" &&
      typeof item[1] === "string" &&
      typeof item[2] === "string"
    );
  };

  /**
   * Fetch file menu items from the API
   */
  const fetchFileMenuItems = async (currentJobId: number) => {
    setLoading(true);
    setError(null);

    try {
      const data: ExportJobFileMenuResponse = await apiGet(
        `jobs/${currentJobId}/export_job_file_menu`
      );

      // Handle new API format: {success: true, data: {result: ...}}
      if (data.success && data.data?.result) {
        let menuItems: FileMenuItem[] = [];
        const resultData = data.data.result;

        // Handle case where result is directly an array of FileMenuItem
        if (Array.isArray(resultData)) {
          menuItems = resultData.filter(isValidFileMenuItem);
        } else if (typeof resultData === "object") {
          // Handle case where result is an object containing arrays
          const resultKeys = Object.keys(resultData);

          // Look for arrays in the result object
          for (const key of resultKeys) {
            const value = resultData[key];
            if (Array.isArray(value)) {
              const validItems = value.filter(isValidFileMenuItem);
              menuItems.push(...validItems);
            }
          }
        }

        // Parse the menu items into structured objects
        const parsedItems = menuItems.map(parseFileMenuItem);
        setFileMenuItems(parsedItems);
      } else {
        throw new Error(data.error || "Failed to fetch file menu items");
      }
    } catch (err) {
      console.error("Error fetching file menu items:", err);
      setError(err instanceof Error ? err.message : "Unknown error occurred");
      setFileMenuItems([]);
    } finally {
      setLoading(false);
    }
  };

  /**
   * Effect to fetch data when jobId changes
   */
  useEffect(() => {
    if (jobId !== null && jobId !== undefined) {
      fetchFileMenuItems(jobId);
    } else {
      // Clear state when dialog closes
      setFileMenuItems([]);
      setError(null);
      setLoading(false);
      setDownloadingItems(new Set());
    }
  }, [jobId]);

  /**
   * Handle dialog close
   */
  const handleClose = () => {
    setJobId(null);
  };

  /**
   * Handle file menu item selection and download
   *
   * Constructs the download URL using the job ID and item identifier,
   * then triggers the download using the doDownload utility function.
   *
   * - Uses proxy endpoint
   * - Handles cloud-specific download patterns
   * - Includes proper error handling for distributed environments
   */
  const handleFileMenuItemClick = async (item: ParsedFileMenuItem) => {
    if (!jobId) {
      console.error("No job ID available for download");
      return;
    }

    // Use identifier as the mode parameter
    const mode = encodeURIComponent(item.identifier);
    const composite_path = `/api/proxy/jobs/${jobId}/export_job_file?mode=${mode}`;

    // Create a descriptive filename from the label and identifier
    const description = item.label || item.identifier;

    // Add to downloading set to show loading state
    setDownloadingItems((prev) => new Set(prev).add(item.identifier));

    try {
      console.log("Initiating download for:", {
        identifier: item.identifier,
        label: item.label,
        mimetype: item.mimetype,
        url: composite_path,
      });

      await doDownload(composite_path, description);

      console.log(`Successfully initiated download for: ${description}`);
    } catch (error) {
      console.error("Download failed:", error);
      setError(
        `Download failed for ${description}: ${
          error instanceof Error ? error.message : "Unknown error"
        }`
      );
    } finally {
      // Remove from downloading set
      setDownloadingItems((prev) => {
        const newSet = new Set(prev);
        newSet.delete(item.identifier);
        return newSet;
      });
    }
  };

  /**
   * Handle bulk download of all available files
   */
  const handleExportAll = async () => {
    if (!jobId || fileMenuItems.length === 0) return;

    setError(null);

    try {
      // Download all files sequentially to avoid overwhelming the server
      for (const item of fileMenuItems) {
        await handleFileMenuItemClick(item);
        // Small delay between downloads to prevent server overload
        await new Promise((resolve) => setTimeout(resolve, 500));
      }
    } catch (error) {
      console.error("Bulk export failed:", error);
      setError(
        `Bulk export failed: ${
          error instanceof Error ? error.message : "Unknown error"
        }`
      );
    }
  };

  /**
   * Get icon based on mimetype
   */
  const getFileIcon = (mimetype: string) => {
    if (mimetype.startsWith("image/")) {
      return <Image sx={{ mr: 2 }} />;
    } else if (mimetype.includes("json") || mimetype.includes("xml")) {
      return <DataObject sx={{ mr: 2 }} />;
    } else if (mimetype.includes("text/")) {
      return <Description sx={{ mr: 2 }} />;
    } else {
      return <InsertDriveFile sx={{ mr: 2 }} />;
    }
  };

  /**
   * Get chip color based on mimetype category
   */
  const getMimetypeChipColor = (mimetype: string) => {
    if (mimetype.startsWith("image/")) return "primary";
    if (mimetype.includes("json") || mimetype.includes("xml"))
      return "secondary";
    if (mimetype.includes("text/")) return "default";
    return "default";
  };

  /**
   * Render file menu item with download functionality
   */
  const renderFileMenuItem = (item: ParsedFileMenuItem, index: number) => {
    const displayName = item.label || item.identifier || `Item ${index + 1}`;
    const description = item.identifier !== item.label ? item.identifier : "";
    const isDownloading = downloadingItems.has(item.identifier);

    return (
      <ListItem key={`${item.identifier}-${index}`} disablePadding>
        <ListItemButton
          onClick={() => handleFileMenuItemClick(item)}
          disabled={isDownloading}
        >
          {isDownloading ? (
            <CircularProgress size={20} sx={{ mr: 2 }} />
          ) : (
            getFileIcon(item.mimetype)
          )}
          <ListItemText
            primary={displayName}
            secondary={description}
            sx={{ flexGrow: 1 }}
          />
          <Chip
            label={item.mimetype}
            size="small"
            variant="outlined"
            color={getMimetypeChipColor(item.mimetype)}
            sx={{ ml: 1 }}
          />
        </ListItemButton>
      </ListItem>
    );
  };

  return (
    <Dialog
      open={isOpen}
      onClose={handleClose}
      maxWidth="md"
      fullWidth
      PaperProps={{
        sx: {
          minHeight: "400px",
          maxHeight: "80vh",
        },
      }}
    >
      <DialogTitle>
        <Box display="flex" alignItems="center" justifyContent="space-between">
          <Box display="flex" alignItems="center">
            <Download sx={{ mr: 1 }} />
            <Typography variant="h6">
              Export Job Files {jobId && `(Job ${jobId})`}
            </Typography>
          </Box>
          <Button
            onClick={handleClose}
            size="small"
            sx={{ minWidth: "auto", p: 1 }}
          >
            <Close />
          </Button>
        </Box>
      </DialogTitle>

      <Divider />

      <DialogContent sx={{ p: 0 }}>
        {loading && (
          <Box
            display="flex"
            justifyContent="center"
            alignItems="center"
            minHeight="200px"
          >
            <CircularProgress />
            <Typography variant="body2" sx={{ ml: 2 }}>
              Loading exportable files...
            </Typography>
          </Box>
        )}

        {error && (
          <Box p={2}>
            <Alert severity="error">
              <Typography variant="body2">{error}</Typography>
            </Alert>
          </Box>
        )}

        {!loading && !error && fileMenuItems.length === 0 && (
          <Box
            display="flex"
            justifyContent="center"
            alignItems="center"
            minHeight="200px"
          >
            <Typography variant="body2" color="text.secondary">
              No exportable files available for this job.
            </Typography>
          </Box>
        )}

        {!loading && !error && fileMenuItems.length > 0 && (
          <List sx={{ py: 0 }}>
            {fileMenuItems.map((item, index) =>
              renderFileMenuItem(item, index)
            )}
          </List>
        )}
      </DialogContent>

      <Divider />

      <DialogActions>
        <Button onClick={handleClose} variant="outlined">
          Close
        </Button>
        {fileMenuItems.length > 0 && (
          <Button
            variant="contained"
            startIcon={<Download />}
            onClick={handleExportAll}
            disabled={downloadingItems.size > 0}
          >
            {downloadingItems.size > 0 ? "Downloading..." : "Export All"}
          </Button>
        )}
      </DialogActions>
    </Dialog>
  );
};

export default ExportJobMenu;
