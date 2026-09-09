"use client";
import React, {
  createContext,
  useCallback,
  useContext,
  useRef,
  useState,
  ReactNode,
} from "react";
import Dialog from "@mui/material/Dialog";
import DialogTitle from "@mui/material/DialogTitle";
import DialogContent from "@mui/material/DialogContent";
import DialogContentText from "@mui/material/DialogContentText";
import DialogActions from "@mui/material/DialogActions";
import TextField from "@mui/material/TextField";
import Button from "@mui/material/Button";
import { readUiPreference } from "@/lib/ui-preferences";

/**
 * Import-provenance capture.
 *
 * When the user imports a local file and the `captureImportProvenance`
 * preference is on, this prompts for a free-text note describing where the
 * file came from — the Qt-era "Describe source of this file" step. The note is
 * carried into the upload POST and stored on the file's import record
 * (FileImport.description), distinct from the auto-generated File.annotation
 * label.
 *
 * `requestImportProvenance(fileName)` resolves to:
 *   - null   → don't attach anything (preference off, or no provider mounted);
 *   - ""     → user chose Skip (imported, no note);
 *   - text   → the user's provenance note.
 *
 * The caller passes the resolved string (when non-null) as `description` to
 * uploadFileParam. Programmatic/derived imports simply never call this, so they
 * are never prompted.
 */
type ImportProvenanceContextType = {
  requestImportProvenance: (fileName: string) => Promise<string | null>;
};

const noop: ImportProvenanceContextType = {
  requestImportProvenance: () => Promise.resolve(null),
};

const ImportProvenanceContext =
  createContext<ImportProvenanceContextType>(noop);

export const useImportProvenance = () => useContext(ImportProvenanceContext);

export const ImportProvenanceProvider: React.FC<{ children: ReactNode }> = ({
  children,
}) => {
  const [open, setOpen] = useState(false);
  const [fileName, setFileName] = useState("");
  const [text, setText] = useState("");
  // The pending promise's resolver, set while the dialog is open.
  const resolverRef = useRef<((value: string | null) => void) | null>(null);

  const requestImportProvenance = useCallback(
    (name: string): Promise<string | null> => {
      // Preference read at call time (not render) so a mid-session toggle takes
      // effect on the next import without this provider re-rendering.
      if (!readUiPreference("captureImportProvenance")) {
        return Promise.resolve(null);
      }
      return new Promise<string | null>((resolve) => {
        resolverRef.current = resolve;
        setFileName(name);
        setText("");
        setOpen(true);
      });
    },
    [],
  );

  const finish = useCallback((value: string | null) => {
    setOpen(false);
    const resolve = resolverRef.current;
    resolverRef.current = null;
    resolve?.(value);
  }, []);

  return (
    <ImportProvenanceContext.Provider value={{ requestImportProvenance }}>
      {children}
      <Dialog
        open={open}
        onClose={() => finish("")}
        maxWidth="sm"
        fullWidth
        // Keep the click that opened the file picker from also dismissing
        // downstream menus oddly; standard modal behaviour is fine here.
      >
        <DialogTitle>Provenance of {fileName || "file"}</DialogTitle>
        <DialogContent>
          <DialogContentText sx={{ mb: 2 }}>
            Describe where this file came from — the experiment, processing
            software, or another project. This note is kept with the imported
            file. You can leave it blank.
          </DialogContentText>
          <TextField
            autoFocus
            fullWidth
            multiline
            minRows={3}
            label="Source of this file"
            value={text}
            onChange={(e) => setText(e.target.value)}
            onKeyDown={(e) => {
              // Cmd/Ctrl+Enter saves; plain Enter keeps making new lines.
              if ((e.metaKey || e.ctrlKey) && e.key === "Enter") {
                e.preventDefault();
                finish(text.trim());
              }
            }}
          />
        </DialogContent>
        <DialogActions>
          <Button onClick={() => finish("")}>Skip</Button>
          <Button variant="contained" onClick={() => finish(text.trim())}>
            Save
          </Button>
        </DialogActions>
      </Dialog>
    </ImportProvenanceContext.Provider>
  );
};
