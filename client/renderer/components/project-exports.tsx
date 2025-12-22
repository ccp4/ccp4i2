import React from "react";
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  IconButton,
  Tooltip,
  LinearProgress,
  Button,
} from "@mui/material";
import { Download, Close, CheckCircle, Delete } from "@mui/icons-material";
import { ProjectExport } from "../types/models";
import { useCCP4i2Window } from "../app-context";
import { useProject } from "../utils";
import { useApi } from "../api";

interface ProjectExportsDialogProps {
  open: boolean;
  onClose: () => void;
}

// Helper function to slugify project name (equivalent to Python's slugify)
const slugify = (text: string): string => {
  return text
    .toString()
    .toLowerCase()
    .trim()
    .replace(/\s+/g, "-") // Replace spaces with -
    .replace(/[^\w\-]+/g, "") // Remove all non-word chars
    .replace(/\-\-+/g, "-") // Replace multiple - with single -
    .replace(/^-+/, "") // Trim - from start of text
    .replace(/-+$/, ""); // Trim - from end of text
};

// Helper function to format timestamp like Python's strftime("%Y%m%d_%H%M%S")
// Uses UTC to match Django server's UTC configuration
const formatTimestamp = (dateString: string): string => {
  const date = new Date(dateString);
  const year = date.getUTCFullYear();
  const month = String(date.getUTCMonth() + 1).padStart(2, "0");
  const day = String(date.getUTCDate()).padStart(2, "0");
  const hours = String(date.getUTCHours()).padStart(2, "0");
  const minutes = String(date.getUTCMinutes()).padStart(2, "0");
  const seconds = String(date.getUTCSeconds()).padStart(2, "0");

  return `${year}${month}${day}_${hours}${minutes}${seconds}`;
};

// Recursive function to search through hierarchical directory structure
const findFileInDirectory = (container: any[], fileName: string): boolean => {
  if (!container || !Array.isArray(container)) {
    return false;
  }

  for (const item of container) {
    // Check if current item matches the filename
    console.log("Checking item:", item.name, fileName);
    if (item.name === fileName) {
      return true;
    }

    // If this item has contents, search recursively
    if (item.contents && Array.isArray(item.contents)) {
      if (findFileInDirectory(item.contents, fileName)) {
        return true;
      }
    }
  }

  return false;
};

export const ProjectExportsDialog: React.FC<ProjectExportsDialogProps> = ({
  open,
  onClose,
}) => {
  const api = useApi();
  const { projectId } = useCCP4i2Window();

  const { project } = useProject(projectId);

  // Use centralized API hook for directory fetching
  const { data: directory, mutate: mutateDirectory } = api.projectDirectory(
    project && open ? project.id : null,
    open
  );

  // Only poll exports when the dialog is open
  const {
    data: exports,
    error,
    isLoading,
    mutate: mutateExports,
  } = api.get<ProjectExport[]>(
    open && projectId ? `projects/${projectId}/exports/` : null,
    open ? 5000 : 0
  );

  //console.log("Project exports:", exports);
  // Force fresh data when dialog opens
  React.useEffect(() => {
    if (open && projectId) {
      // Clear cache and force fresh fetch
      mutateExports(undefined, { revalidate: true });
      if (project) {
        mutateDirectory(undefined, { revalidate: true });
      }
    }
  }, [open, projectId, project?.id, mutateExports, mutateDirectory]);

  // State for delete confirmation dialog
  const [deleteDialogOpen, setDeleteDialogOpen] = React.useState(false);
  const [exportToDelete, setExportToDelete] =
    React.useState<ProjectExport | null>(null);

  const formatDateTime = (dateString: string) => {
    const date = new Date(dateString);
    return date.toLocaleString();
  };

  const handleDownload = (exportItem: ProjectExport) => {
    // Create a download link for the export file
    const downloadUrl = `/api/proxy/projectexports/${exportItem.id}/download/`;
    const link = document.createElement("a");
    const projectName =
      typeof exportItem.project === "object"
        ? exportItem.project.name
        : `project_${exportItem.project}`;
    link.href = downloadUrl;
    link.download = `${projectName}_export_${new Date(exportItem.time).toISOString().slice(0, 19).replace(/:/g, "")}.ccp4_project.zip`;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  };

  const handleDeleteClick = (exportItem: ProjectExport) => {
    setExportToDelete(exportItem);
    setDeleteDialogOpen(true);
  };

  const handleDeleteConfirm = async () => {
    if (!exportToDelete) return;

    try {
      await api.delete(`projectexports/${exportToDelete.id}/`);
      // Refresh both the exports list and project directory
      mutateExports();
      // Note: Project directory will auto-refresh due to the 5-second interval
      setDeleteDialogOpen(false);
      setExportToDelete(null);
    } catch (error) {
      console.error("Failed to delete export:", error);
      // You might want to show an error message to the user here
    }
  };

  const handleDeleteCancel = () => {
    setDeleteDialogOpen(false);
    setExportToDelete(null);
  };

  const inferredNames = React.useMemo(() => {
    if (!project || !exports) return {};
    const names: { [key: number]: string } = {};
    exports.forEach((exportItem) => {
      const projectName = slugify(project.name);
      const timestamp = formatTimestamp(exportItem.time);
      const inferredName = `${projectName}_export_${timestamp}.ccp4_project.zip`;
      names[exportItem.id] = inferredName;
    });
    return names;
  }, [exports, project]);

  // Memoized function to check if files exist in the hierarchical directory structure
  const fileExistence = React.useMemo(() => {
    if (!directory?.container || !exports || !project) return {};

    const existence: { [key: number]: boolean } = {};
    exports.forEach((exportItem) => {
      const inferredName = inferredNames[exportItem.id];
      if (inferredName) {
        existence[exportItem.id] = findFileInDirectory(
          directory.container,
          inferredName
        );
      }
    });
    return existence;
  }, [directory?.container, exports, project, inferredNames]);

  return (
    <Dialog open={open} onClose={onClose} maxWidth="md" fullWidth>
      <DialogTitle>
        Project Exports - {project ? project.name : "Loading..."}
        <IconButton
          onClick={onClose}
          sx={{ position: "absolute", right: 8, top: 8 }}
        >
          <Close />
        </IconButton>
      </DialogTitle>
      <DialogContent>
        {isLoading && <LinearProgress />}
        {error && <div>Error loading exports</div>}
        {exports && exports.length === 0 && <div>No exports found</div>}
        {exports && exports.length > 0 && project && inferredNames && (
          <Table size="small">
            <TableHead>
              <TableRow>
                <TableCell>Exists</TableCell>
                <TableCell>File Name</TableCell>
                <TableCell>Export Time</TableCell>
                <TableCell>Actions</TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              {exports.map((exportItem) => (
                <TableRow key={exportItem.id}>
                  <TableCell>
                    {fileExistence[exportItem.id] ? (
                      <Tooltip title="File exists">
                        <CheckCircle color="success" />
                      </Tooltip>
                    ) : (
                      <Tooltip title="File missing">
                        <Close color="error" />
                      </Tooltip>
                    )}
                  </TableCell>
                  <TableCell>{inferredNames[exportItem.id]}</TableCell>
                  <TableCell>{formatDateTime(exportItem.time)}</TableCell>
                  <TableCell>
                    <Tooltip title="Download export">
                      <IconButton onClick={() => handleDownload(exportItem)}>
                        <Download />
                      </IconButton>
                    </Tooltip>
                    <Tooltip title="Delete export">
                      <IconButton
                        onClick={() => handleDeleteClick(exportItem)}
                        color="error"
                      >
                        <Delete />
                      </IconButton>
                    </Tooltip>
                  </TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>
        )}
      </DialogContent>

      {/* Delete Confirmation Dialog */}
      <Dialog
        open={deleteDialogOpen}
        onClose={handleDeleteCancel}
        maxWidth="sm"
        fullWidth
      >
        <DialogTitle>Confirm Delete</DialogTitle>
        <DialogContent>
          Are you sure you want to delete this export? This action cannot be
          undone and will also remove the associated files from disk.
        </DialogContent>
        <DialogActions>
          <Button onClick={handleDeleteCancel} color="inherit">
            Cancel
          </Button>
          <Button
            onClick={handleDeleteConfirm}
            color="error"
            variant="contained"
          >
            Delete
          </Button>
        </DialogActions>
      </Dialog>
    </Dialog>
  );
};
