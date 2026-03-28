/*
 * Copyright (C) 2025-2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
import {
  SyntheticEvent,
  useCallback,
  useContext,
  createContext,
  useState,
  useMemo,
  PropsWithChildren,
  useRef,
  useEffect,
} from "react";
import { doDownload, useApi } from "../api";
import {
  Divider,
  ListItemIcon,
  ListItemText,
  Menu,
  MenuItem,
  Popper,
  Paper,
  TextField,
  ClickAwayListener,
  Box,
} from "@mui/material";
import { ContentCopy, Download, SaveAlt, Preview, Terminal, Edit, Image as ImageIcon, Science } from "@mui/icons-material";
import { useFilePreviewContext } from "./file-preview-context";
import { File as DjangoFile } from "../types/models";
import { useRouter } from "next/navigation";
import { CCP4i2MoorhenIcon } from "../components/General/CCP4i2Icons";
import { useCCP4i2Window } from "../app-context";
import { TableChart } from "@mui/icons-material";

export interface FileMenuExtraItem {
  key: string;
  label: string;
  icon?: React.ReactNode;
  onClick: () => void;
  disabled?: boolean;
  divider?: boolean;
}

interface FileMenuContextProps {
  fileMenuAnchorEl: HTMLElement | null;
  setFileMenuAnchorEl: (element: HTMLElement | null) => void;
  file: DjangoFile | null;
  setFile: (file: DjangoFile | null) => void;
  extraMenuItems: FileMenuExtraItem[];
  setExtraMenuItems: (items: FileMenuExtraItem[]) => void;
}

export const FileMenuContext = createContext<FileMenuContextProps>({
  fileMenuAnchorEl: null,
  setFileMenuAnchorEl: () => {},
  file: null,
  setFile: () => {},
  extraMenuItems: [],
  setExtraMenuItems: () => {},
});

export const FileMenuProvider: React.FC<PropsWithChildren> = ({ children }) => {
  const [fileMenuAnchorEl, setFileMenuAnchorEl] = useState<HTMLElement | null>(
    null
  );
  const [file, setFile] = useState<DjangoFile | null>(null);
  const [extraMenuItems, setExtraMenuItems] = useState<FileMenuExtraItem[]>([]);

  const contextValue = useMemo(
    () => ({
      fileMenuAnchorEl,
      setFileMenuAnchorEl,
      file,
      setFile,
      extraMenuItems,
      setExtraMenuItems,
    }),
    [fileMenuAnchorEl, file, extraMenuItems]
  );

  return (
    <FileMenuContext.Provider value={contextValue}>
      {children}
      <FileMenu />
    </FileMenuContext.Provider>
  );
};

/** Check if a file type represents an MTZ file (including unmerged variants) */
function isMtzFileType(type: string): boolean {
  return type.startsWith("application/CCP4-mtz") ||
    type.startsWith("application/CCP4-unmerged");
}

export const FileMenu: React.FC = () => {
  const { fileMenuAnchorEl, setFileMenuAnchorEl, file, extraMenuItems, setExtraMenuItems } = useFileMenu();
  const { projectId } = useCCP4i2Window();
  const api = useApi();
  const { setContentSpecification } = useFilePreviewContext();
  const router = useRouter();
  const { mutate: mutateFiles } = api.get_endpoint<DjangoFile[]>({
    type: "projects",
    id: projectId || 0,
    endpoint: "files",
  });

  // State for annotation editor popper
  const [annotationPopperAnchorEl, setAnnotationPopperAnchorEl] =
    useState<HTMLElement | null>(null);
  const [annotationValue, setAnnotationValue] = useState<string>("");

  // Store the original annotation value when popper opens
  const originalAnnotationValue = useRef<string>("");

  // Ref to store the current timeout ID for cleanup
  const timeoutRef = useRef<NodeJS.Timeout | null>(null);

  // Function to save annotation immediately
  const saveAnnotation = useCallback(
    async (value: string) => {
      if (!file) return;

      try {
        await api.patch(`files/${file.id}`, {
          annotation: value,
        });
        mutateFiles(); // Refresh file list after saving annotation
      } catch (error) {
        console.error("Failed to update annotation:", error);
      }
    },
    [file, api, mutateFiles]
  );

  // Handle Enter and Escape key presses in the text field
  const handleKeyDown = useCallback(
    async (event: React.KeyboardEvent<HTMLDivElement>) => {
      if (event.key === "Enter" && !event.shiftKey) {
        event.preventDefault();

        // Clear any pending timeout
        if (timeoutRef.current) {
          clearTimeout(timeoutRef.current);
          timeoutRef.current = null;
        }

        // Save immediately and close popper
        await saveAnnotation(annotationValue);
        setAnnotationPopperAnchorEl(null);
      } else if (event.key === "Escape") {
        event.preventDefault();

        // Clear any pending timeout
        if (timeoutRef.current) {
          clearTimeout(timeoutRef.current);
          timeoutRef.current = null;
        }

        // Restore original value and close popper without saving
        await saveAnnotation(originalAnnotationValue.current);
        setAnnotationPopperAnchorEl(null);
      }
    },
    [annotationValue, saveAnnotation]
  );

  // Handle annotation text change with delay
  const handleAnnotationChange = useCallback(
    (newValue: string) => {
      setAnnotationValue(newValue);

      // Clear existing timeout
      if (timeoutRef.current) {
        clearTimeout(timeoutRef.current);
      }

      // Create a delayed function to save annotation
      timeoutRef.current = setTimeout(() => {
        saveAnnotation(newValue);
        timeoutRef.current = null;
      }, 500); // 500ms delay
    },
    [saveAnnotation]
  );

  // Handle popper close
  const handleAnnotationPopperClose = useCallback(() => {
    // Clear any pending timeout when closing
    if (timeoutRef.current) {
      clearTimeout(timeoutRef.current);
      timeoutRef.current = null;
    }
    setAnnotationPopperAnchorEl(null);
  }, []);

  // Handle click away from popper
  const handleClickAway = useCallback(
    (event: MouseEvent | TouchEvent) => {
      handleAnnotationPopperClose();
    },
    [handleAnnotationPopperClose]
  );

  // Cleanup timeout on unmount
  useEffect(() => {
    return () => {
      if (timeoutRef.current) {
        clearTimeout(timeoutRef.current);
      }
    };
  }, []);

  const handleDownloadFile = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (file) {
        const composite_path = api.noSlashUrl(`files/${file.id}/download/`);
        doDownload(composite_path, file.name);
        setFileMenuAnchorEl(null);
      }
    },
    [file, api, setFileMenuAnchorEl]
  );

  /** Save file to user-chosen location via native dialog (Electron only). */
  const handleSaveFileAs = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (file?.path && window.electronAPI) {
        window.electronAPI.sendMessage("save-file-as", {
          filePath: file.path,
          fileName: file.name,
        });
      }
      setFileMenuAnchorEl(null);
    },
    [file, setFileMenuAnchorEl]
  );

  /** Copy file reference to clipboard as JSON (for cross-window paste). */
  const handleCopyReference = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (file) {
        const ref = {
          ccp4i2_file: true,
          uuid: file.uuid,
          id: file.id,
          name: file.name,
          type: file.type,
          sub_type: file.sub_type,
          content: file.content,
          annotation: file.annotation,
          job: file.job,
          job_param_name: file.job_param_name,
        };
        await navigator.clipboard.writeText(JSON.stringify(ref));
        setFileMenuAnchorEl(null);
      }
    },
    [file, setFileMenuAnchorEl]
  );

  const handlePreviewFile = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (file) {
        let language = "text";
        let url = `/api/proxy/ccp4i2/files/${file.id}/download/`;

        if (isMtzFileType(file.type)) {
          language = "mtz";
        } else if (file.type === "application/CCP4-image") {
          language = "image";
        } else if (file.type === "application/refmac-dictionary") {
          language = "dict-preview";
          url = `/api/proxy/ccp4i2/files/${file.id}/digest/`;
        } else if (file.type === "chemical/x-mdl-molfile") {
          language = "molblock-raw";
        } else if (file.type === "application/CCP4-seqalign") {
          language = "clustalw";
        }

        setContentSpecification({ url, title: file.name, language });
        setFileMenuAnchorEl(null);
      }
    },
    [file, setContentSpecification, setFileMenuAnchorEl]
  );

  const handlePreviewFileDigest = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (file) {
        setContentSpecification({
          url: `/api/proxy/ccp4i2/files/${file.id}/digest/`,
          title: file.name,
          language: "json",
        });
        setFileMenuAnchorEl(null);
      }
    },
    [file, setContentSpecification, setFileMenuAnchorEl]
  );

  const handlePreviewDbInfo = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (file) {
        setContentSpecification({
          url: `/api/proxy/ccp4i2/files/${file.id}/`,
          title: file.name,
          language: "json",
        });
        setFileMenuAnchorEl(null);
      }
    },
    [file, setContentSpecification, setFileMenuAnchorEl]
  );

  const handlePreviewFileInCoot = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (file) {
        api.post<any>(`files/${file.id}/preview/`, { viewer: "coot" });
        setFileMenuAnchorEl(null);
      }
    },
    [file, api, setFileMenuAnchorEl]
  );

  const handlePreviewFileInViewHKL = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (file) {
        api.post<any>(`files/${file.id}/preview/`, { viewer: "viewhkl" });
        setFileMenuAnchorEl(null);
      }
    },
    [file, api, setFileMenuAnchorEl]
  );

  const handlePreviewFileInCCP4MG = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (file) {
        api.post<any>(`files/${file.id}/preview/`, { viewer: "ccp4mg" });
        setFileMenuAnchorEl(null);
      }
    },
    [file, api, setFileMenuAnchorEl]
  );

  const handleOpenInNewWindow = (path: string) => {
    window.open(path, "_blank", "noopener,noreferrer");
  };

  const handlePreviewFileInMoorhen = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (file) {
        handleOpenInNewWindow(`/ccp4i2/moorhen-page/file-by-id/${file.id}`);
        setFileMenuAnchorEl(null);
      }
    },
    [file, setFileMenuAnchorEl]
  );

  const handlePreviewFileInTerminal = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (file) {
        api.post<any>(`files/${file.id}/preview/`, { viewer: "terminal" });
        setFileMenuAnchorEl(null);
      }
    },
    [file, api, setFileMenuAnchorEl]
  );

  const handlePreviewFileInPostScript = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (file) {
        api.post<any>(`files/${file.id}/preview/`, { viewer: "postscript" });
        setFileMenuAnchorEl(null);
      }
    },
    [file, api, setFileMenuAnchorEl]
  );

  const handlePreviewImage = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (file) {
        setContentSpecification({
          url: `/api/proxy/ccp4i2/files/${file.id}/download/`,
          title: file.name,
          language: "image",
        });
        setFileMenuAnchorEl(null);
      }
    },
    [file, setContentSpecification, setFileMenuAnchorEl]
  );

  const handlePreviewMolblock = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (file) {
        setContentSpecification({
          url: `/api/proxy/ccp4i2/files/${file.id}/molblock/`,
          title: `${file.name} - 2D Structure`,
          language: "molblock",
        });
        setFileMenuAnchorEl(null);
      }
    },
    [file, setContentSpecification, setFileMenuAnchorEl]
  );

  const handlePreviewMolfileRaw = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (file) {
        setContentSpecification({
          url: `/api/proxy/ccp4i2/files/${file.id}/download/`,
          title: `${file.name} - 2D Structure`,
          language: "molblock-raw",
        });
        setFileMenuAnchorEl(null);
      }
    },
    [file, setContentSpecification, setFileMenuAnchorEl]
  );

  // Preview MTZ header using pure TypeScript parser
  const handlePreviewMtzHeader = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (file) {
        setContentSpecification({
          url: `/api/proxy/ccp4i2/files/${file.id}/download/`,
          title: file.name,
          language: "mtz",
        });
        setFileMenuAnchorEl(null);
      }
    },
    [file, setContentSpecification, setFileMenuAnchorEl]
  );

  // Handle edit annotation menu item click
  const handleEditAnnotation = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (file) {
        // Store the original annotation value for potential restoration
        const currentAnnotation = file.annotation || "";
        originalAnnotationValue.current = currentAnnotation;

        // Set current annotation value
        setAnnotationValue(currentAnnotation);
        // Use the menu anchor element as the popper anchor
        setAnnotationPopperAnchorEl(fileMenuAnchorEl);
        // Close the menu
        setFileMenuAnchorEl(null);
      }
    },
    [file, fileMenuAnchorEl, setFileMenuAnchorEl]
  );

  return (
    <>
      <Menu
        open={Boolean(fileMenuAnchorEl)}
        anchorEl={fileMenuAnchorEl}
        onClose={() => { setFileMenuAnchorEl(null); setExtraMenuItems([]); }}
      >
        {/* Context-specific items injected by the caller (e.g. Clear, Paste, subtype items) */}
        {extraMenuItems.length > 0 && [
          ...extraMenuItems.map((extra) => [
            extra.divider && <Divider key={`div-${extra.key}`} />,
            <MenuItem
              key={extra.key}
              onClick={() => { setFileMenuAnchorEl(null); setExtraMenuItems([]); extra.onClick(); }}
              disabled={extra.disabled}
            >
              {extra.icon && <ListItemIcon>{extra.icon}</ListItemIcon>}
              <ListItemText>{extra.label}</ListItemText>
            </MenuItem>,
          ]),
          <Divider key="extra-divider" />,
        ]}
        {file?.path && typeof window !== "undefined" && window.electronAPI ? (
          <MenuItem key="SaveAs" onClick={handleSaveFileAs}>
            <SaveAlt sx={{ mr: 1 }} /> Save to...
          </MenuItem>
        ) : (
          <MenuItem key="Download" onClick={handleDownloadFile}>
            <Download sx={{ mr: 1 }} /> Download
          </MenuItem>
        )}
        <MenuItem key="CopyRef" onClick={handleCopyReference}>
          <ContentCopy sx={{ mr: 1 }} /> Copy reference
        </MenuItem>
        {(!file || !isMtzFileType(file.type)) && (
          <MenuItem key="Preview" onClick={handlePreviewFile}>
            <Preview sx={{ mr: 1 }} /> Preview
          </MenuItem>
        )}
        <MenuItem key="Terminal" onClick={handlePreviewFileInTerminal}>
          <Terminal sx={{ mr: 1 }} /> Terminal
        </MenuItem>
        <MenuItem key="EditAnnotation" onClick={handleEditAnnotation}>
          <Edit sx={{ mr: 1 }} /> Edit annotation
        </MenuItem>
        {file && (
          <MenuItem key="DbInfo" onClick={handlePreviewDbInfo}>
            <Preview sx={{ mr: 1 }} /> DbInfo
          </MenuItem>
        )}
        {file &&
          ["chemical/x-pdb", "application/CCP4-mtz-map"].includes(
            file.type
          ) && (
            <MenuItem key="Coot" onClick={handlePreviewFileInCoot}>
              <Preview sx={{ mr: 1 }} /> Coot
            </MenuItem>
          )}
        {file &&
          ["chemical/x-pdb", "application/CCP4-mtz-map"].includes(
            file.type
          ) && (
            <MenuItem key="CCP4MG" onClick={handlePreviewFileInCCP4MG}>
              <Preview sx={{ mr: 1 }} /> CCP4MG
            </MenuItem>
          )}
        {file &&
          [
            "chemical/x-pdb",
            "application/CCP4-mtz-map",
            "application/refmac-dictionary",
          ].includes(file.type) && (
            <MenuItem key="Moorhen" onClick={handlePreviewFileInMoorhen}>
              <CCP4i2MoorhenIcon sx={{ mr: 1 }} /> Moorhen
            </MenuItem>
          )}
        {file && isMtzFileType(file.type) && (
          <MenuItem key="ViewHKL" onClick={handlePreviewFileInViewHKL}>
            <Preview sx={{ mr: 1 }} /> ViewHKL
          </MenuItem>
        )}
        {file && isMtzFileType(file.type) && (
          <MenuItem key="MtzHeader" onClick={handlePreviewMtzHeader}>
            <TableChart sx={{ mr: 1 }} /> MTZ Preview
          </MenuItem>
        )}
        {file && file.type === "application/CCP4-image" && (
          <MenuItem key="ViewImage" onClick={handlePreviewImage}>
            <ImageIcon sx={{ mr: 1 }} /> View Image
          </MenuItem>
        )}
        {file && file.type === "application/refmac-dictionary" && (
          <MenuItem key="ViewStructure" onClick={handlePreviewMolblock}>
            <Science sx={{ mr: 1 }} /> View 2D Structure
          </MenuItem>
        )}
        {file && file.type === "chemical/x-mdl-molfile" && (
          <MenuItem key="ViewMolfile" onClick={handlePreviewMolfileRaw}>
            <Science sx={{ mr: 1 }} /> View 2D Structure
          </MenuItem>
        )}
        {file &&
          ["application/postscript", "application/x-pdf"].includes(
            file.type
          ) && (
            <MenuItem key="PostScript" onClick={handlePreviewFileInPostScript}>
              <Preview sx={{ mr: 1 }} /> Document Viewer
            </MenuItem>
          )}
        {file && (
          <MenuItem key="DIGEST" onClick={handlePreviewFileDigest}>
            <Preview sx={{ mr: 1 }} /> DIGEST
          </MenuItem>
        )}
      </Menu>

      {/* Annotation Editor Popper */}
      <Popper
        open={Boolean(annotationPopperAnchorEl)}
        anchorEl={annotationPopperAnchorEl}
        placement="bottom-start"
        sx={{ zIndex: 1300 }}
      >
        <ClickAwayListener onClickAway={handleClickAway}>
          <Paper
            elevation={8}
            sx={{
              p: 2,
              minWidth: 300,
              maxWidth: 500,
            }}
          >
            <Box>
              <TextField
                fullWidth
                multiline
                rows={3}
                label="File Annotation"
                value={annotationValue}
                onChange={(e) => handleAnnotationChange(e.target.value)}
                onKeyDown={handleKeyDown}
                placeholder="Enter annotation for this file..."
                variant="outlined"
                size="small"
                autoFocus
                helperText="Press Enter to save and close, Escape to cancel, or wait 500ms for auto-save"
              />
            </Box>
          </Paper>
        </ClickAwayListener>
      </Popper>
    </>
  );
};

export const useFileMenu = () => {
  const context = useContext(FileMenuContext);
  if (!context) {
    throw new Error("useFileMenu must be used within a FileMenuProvider");
  }
  return context;
};
