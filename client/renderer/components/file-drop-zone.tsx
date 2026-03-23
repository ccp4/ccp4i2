"use client";

import { useCallback, useState } from "react";
import {
  Box,
  Chip,
  IconButton,
  Paper,
  Stack,
  Typography,
} from "@mui/material";
import {
  CloudUpload as UploadIcon,
  Clear as ClearIcon,
} from "@mui/icons-material";

/**
 * File type detection based on extension.
 * Returns the task name to run for importing this file type.
 */
export type DetectedFileType =
  | "reflections"    // .mtz -> splitMtz
  | "coordinates"    // .pdb, .ent -> coordinate_selector
  | "mmcif_coords"   // .cif (coordinates) -> coordinate_selector
  | "sequence"       // .fasta, .fa, .seq, .pir -> ProvideSequence
  | "alignment"      // .aln, .clw, .sto -> ProvideAlignment
  | "unknown";

export interface DroppedFile {
  file: File;
  detectedType: DetectedFileType;
}

const EXTENSION_MAP: Record<string, DetectedFileType> = {
  ".mtz": "reflections",
  ".pdb": "coordinates",
  ".ent": "coordinates",
  ".fasta": "sequence",
  ".fa": "sequence",
  ".seq": "sequence",
  ".pir": "sequence",
  ".aln": "alignment",
  ".clw": "alignment",
  ".sto": "alignment",
  ".stk": "alignment",
  ".phylip": "alignment",
  ".phy": "alignment",
};

/** Map detected file type to the task that will import it */
export const TASK_FOR_TYPE: Record<DetectedFileType, string | null> = {
  reflections: "splitMtz",
  coordinates: "coordinate_selector",
  mmcif_coords: "coordinate_selector",
  sequence: "ProvideSequence",
  alignment: "ProvideAlignment",
  unknown: null,
};

/** Map detected file type to the input parameter path for file upload */
export const PARAM_FOR_TYPE: Record<DetectedFileType, string | null> = {
  reflections: "inputData.HKLIN",
  coordinates: "inputData.XYZIN",
  mmcif_coords: "inputData.XYZIN",
  sequence: "controlParameters.SEQUENCETEXT",
  alignment: "inputData.ALIGNIN",
  unknown: null,
};

const TYPE_LABELS: Record<DetectedFileType, string> = {
  reflections: "Reflections (MTZ)",
  coordinates: "Coordinates",
  mmcif_coords: "Coordinates (mmCIF)",
  sequence: "Sequence",
  alignment: "Alignment",
  unknown: "Unknown",
};

const TYPE_COLORS: Record<DetectedFileType, "primary" | "secondary" | "success" | "warning" | "error" | "info" | "default"> = {
  reflections: "primary",
  coordinates: "success",
  mmcif_coords: "success",
  sequence: "info",
  alignment: "secondary",
  unknown: "default",
};

function detectFileType(file: File): DetectedFileType {
  const name = file.name.toLowerCase();
  const ext = "." + name.split(".").pop();

  // Direct extension match
  if (ext in EXTENSION_MAP) {
    return EXTENSION_MAP[ext];
  }

  // .cif is ambiguous: could be coordinates or reflections
  // Default to coordinates (more common for drag-and-drop)
  if (ext === ".cif" || ext === ".mmcif") {
    return "mmcif_coords";
  }

  return "unknown";
}

interface FileDropZoneProps {
  files: DroppedFile[];
  onChange: (files: DroppedFile[]) => void;
}

export const FileDropZone: React.FC<FileDropZoneProps> = ({
  files,
  onChange,
}) => {
  const [isDragOver, setIsDragOver] = useState(false);

  const handleDrop = useCallback(
    (e: React.DragEvent) => {
      e.preventDefault();
      setIsDragOver(false);
      const droppedFiles = e.dataTransfer.files;
      if (!droppedFiles) return;

      const newFiles: DroppedFile[] = Array.from(droppedFiles).map((file) => ({
        file,
        detectedType: detectFileType(file),
      }));

      onChange([...files, ...newFiles]);
    },
    [files, onChange]
  );

  const handleFileInput = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      const selectedFiles = e.target.files;
      if (!selectedFiles) return;

      const newFiles: DroppedFile[] = Array.from(selectedFiles).map((file) => ({
        file,
        detectedType: detectFileType(file),
      }));

      onChange([...files, ...newFiles]);
      // Reset the input so the same file can be selected again
      e.target.value = "";
    },
    [files, onChange]
  );

  const removeFile = useCallback(
    (index: number) => {
      onChange(files.filter((_, i) => i !== index));
    },
    [files, onChange]
  );

  return (
    <Stack spacing={1}>
      <Typography variant="subtitle2" color="text.secondary">
        Starting data (optional)
      </Typography>
      <Paper
        variant="outlined"
        sx={{
          p: files.length > 0 ? 1.5 : 4,
          border: "2px dashed",
          borderColor: isDragOver ? "primary.main" : "divider",
          backgroundColor: isDragOver ? "action.hover" : "transparent",
          textAlign: "center",
          cursor: "pointer",
          transition: "all 0.2s ease",
        }}
        onDragOver={(e) => {
          e.preventDefault();
          setIsDragOver(true);
        }}
        onDragLeave={() => setIsDragOver(false)}
        onDrop={handleDrop}
        onClick={() => {
          const input = document.createElement("input");
          input.type = "file";
          input.multiple = true;
          input.accept = ".mtz,.pdb,.ent,.cif,.mmcif,.fasta,.fa,.seq,.pir,.aln,.clw,.sto";
          input.onchange = (e) =>
            handleFileInput(e as unknown as React.ChangeEvent<HTMLInputElement>);
          input.click();
        }}
      >
        {files.length === 0 ? (
          <Stack spacing={1} alignItems="center">
            <UploadIcon sx={{ fontSize: 32, color: "text.secondary" }} />
            <Typography variant="body2" color="text.secondary">
              Drop files here to import after project creation
            </Typography>
            <Typography variant="caption" color="text.secondary">
              MTZ, PDB, CIF, FASTA, sequence alignments
            </Typography>
          </Stack>
        ) : (
          <Stack spacing={0.5} onClick={(e) => e.stopPropagation()}>
            {files.map((df, idx) => (
              <Box
                key={idx}
                sx={{
                  display: "flex",
                  alignItems: "center",
                  justifyContent: "space-between",
                }}
              >
                <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
                  <Chip
                    label={TYPE_LABELS[df.detectedType]}
                    color={TYPE_COLORS[df.detectedType]}
                    size="small"
                    variant="outlined"
                  />
                  <Typography variant="body2" noWrap sx={{ maxWidth: 300 }}>
                    {df.file.name}
                  </Typography>
                </Box>
                <IconButton
                  size="small"
                  onClick={(e) => {
                    e.stopPropagation();
                    removeFile(idx);
                  }}
                >
                  <ClearIcon fontSize="small" />
                </IconButton>
              </Box>
            ))}
            <Typography
              variant="caption"
              color="text.secondary"
              sx={{ cursor: "pointer", mt: 0.5 }}
              onClick={() => {
                const input = document.createElement("input");
                input.type = "file";
                input.multiple = true;
                input.accept =
                  ".mtz,.pdb,.ent,.cif,.mmcif,.fasta,.fa,.seq,.pir,.aln,.clw,.sto";
                input.onchange = (e) =>
                  handleFileInput(
                    e as unknown as React.ChangeEvent<HTMLInputElement>
                  );
                input.click();
              }}
            >
              + Drop or click to add more files
            </Typography>
          </Stack>
        )}
      </Paper>
    </Stack>
  );
};
