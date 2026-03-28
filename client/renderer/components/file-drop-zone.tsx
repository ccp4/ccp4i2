/*
 * Copyright (C) 2026 Newcastle University
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
"use client";

import { useCallback, useState } from "react";
import {
  Box,
  Chip,
  CircularProgress,
  IconButton,
  Paper,
  Stack,
  Typography,
} from "@mui/material";
import {
  CloudUpload as UploadIcon,
  Clear as ClearIcon,
} from "@mui/icons-material";
import { parseMtzFile } from "../lib/mtz-parser";

/**
 * Detected file types and their import handling.
 *
 * Types with autoRun=true are simple imports that run immediately.
 * Types with autoRun=false create the job with data loaded but leave
 * it for the user to review parameters before running.
 */
export type DetectedFileType =
  | "reflections"       // merged .mtz -> splitMtz (auto-run)
  | "unmerged"          // unmerged .mtz, .sca, .hkl, .refl -> aimless_pipe (no auto-run)
  | "coordinates"       // .pdb, .ent -> coordinate_selector (auto-run)
  | "mmcif_coords"      // .cif (coordinates) -> coordinate_selector (auto-run)
  | "sequence"          // .fasta, .fa, .seq, .pir -> ProvideSequence (auto-run)
  | "alignment"         // .aln, .clw, .sto -> ProvideAlignment (auto-run)
  | "ligand"            // .smi, .mol, .mol2 -> LidiaAcedrgNew (auto-run)
  | "unknown";

export interface DroppedFile {
  file: File;
  detectedType: DetectedFileType;
}

/** Map detected file type to the task that will import it */
export const TASK_FOR_TYPE: Record<DetectedFileType, string | null> = {
  reflections: "splitMtz",
  unmerged: "aimless_pipe",
  coordinates: "coordinate_selector",
  mmcif_coords: "coordinate_selector",
  sequence: "ProvideSequence",
  alignment: "ProvideAlignment",
  ligand: "LidiaAcedrgNew",
  unknown: null,
};

/** Map detected file type to the input parameter path for file upload */
export const PARAM_FOR_TYPE: Record<DetectedFileType, string | null> = {
  reflections: "inputData.HKLIN",
  unmerged: "inputData.UNMERGEDFILES",
  coordinates: "inputData.XYZIN",
  mmcif_coords: "inputData.XYZIN",
  sequence: "controlParameters.SEQUENCETEXT",
  alignment: "inputData.ALIGNIN",
  ligand: "inputData.MOLIN",
  unknown: null,
};

/** Whether the job should auto-run after creation */
export const AUTO_RUN_FOR_TYPE: Record<DetectedFileType, boolean> = {
  reflections: true,
  unmerged: false,   // User may need to set reference, space group, etc.
  coordinates: true,
  mmcif_coords: true,
  sequence: true,
  alignment: true,
  ligand: true,
  unknown: false,
};

const TYPE_LABELS: Record<DetectedFileType, string> = {
  reflections: "Reflections (merged)",
  unmerged: "Unmerged data",
  coordinates: "Coordinates",
  mmcif_coords: "Coordinates (mmCIF)",
  sequence: "Sequence",
  alignment: "Alignment",
  ligand: "Ligand",
  unknown: "Unknown",
};

const TYPE_COLORS: Record<
  DetectedFileType,
  "primary" | "secondary" | "success" | "warning" | "error" | "info" | "default"
> = {
  reflections: "primary",
  unmerged: "warning",
  coordinates: "success",
  mmcif_coords: "success",
  sequence: "info",
  alignment: "secondary",
  ligand: "secondary",
  unknown: "default",
};

/** Extensions that are unambiguously unmerged data */
const UNMERGED_EXTENSIONS: Record<string, boolean> = {
  ".sca": true,
  ".hkl": true,
  ".refl": true,
};

const EXTENSION_MAP: Record<string, DetectedFileType> = {
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
  ".smi": "ligand",
  ".mol": "ligand",
  ".mol2": "ligand",
};

/**
 * Detect file type from extension alone (synchronous, for non-MTZ files).
 * MTZ files need async sniffing — see detectFileTypeAsync.
 */
function detectFileTypeSync(file: File): DetectedFileType {
  const name = file.name.toLowerCase();
  const ext = "." + name.split(".").pop();

  if (ext in UNMERGED_EXTENSIONS) return "unmerged";
  if (ext in EXTENSION_MAP) return EXTENSION_MAP[ext];

  // .cif is ambiguous — default to coordinates
  if (ext === ".cif" || ext === ".mmcif") return "mmcif_coords";

  // .mtz needs header sniffing — mark as reflections initially,
  // will be refined asynchronously
  if (ext === ".mtz") return "reflections";

  return "unknown";
}

/**
 * For MTZ files, read the header to determine merged vs unmerged.
 * Returns the refined type, or the original type for non-MTZ files.
 */
async function sniffMtzFile(file: File): Promise<DetectedFileType> {
  try {
    const header = await parseMtzFile(file);
    return header.isMerged ? "reflections" : "unmerged";
  } catch {
    // If parsing fails, assume merged
    return "reflections";
  }
}

const ACCEPTED_EXTENSIONS =
  ".mtz,.pdb,.ent,.cif,.mmcif,.fasta,.fa,.seq,.pir,.aln,.clw,.sto,.sca,.hkl,.refl,.smi,.mol,.mol2";

interface FileDropZoneProps {
  files: DroppedFile[];
  onChange: (files: DroppedFile[]) => void;
}

export const FileDropZone: React.FC<FileDropZoneProps> = ({
  files,
  onChange,
}) => {
  const [isDragOver, setIsDragOver] = useState(false);
  const [isSniffing, setIsSniffing] = useState(false);

  const addFiles = useCallback(
    async (fileList: FileList | File[]) => {
      const rawFiles = Array.from(fileList);

      // Initial sync classification
      const newFiles: DroppedFile[] = rawFiles.map((file) => ({
        file,
        detectedType: detectFileTypeSync(file),
      }));

      // Immediately show files with sync-detected types
      const updatedFiles = [...files, ...newFiles];
      onChange(updatedFiles);

      // Async sniff MTZ files to distinguish merged/unmerged
      const mtzIndices = newFiles
        .map((df, i) => (df.file.name.toLowerCase().endsWith(".mtz") ? i : -1))
        .filter((i) => i >= 0);

      if (mtzIndices.length > 0) {
        setIsSniffing(true);
        const refined = [...updatedFiles];
        for (const localIdx of mtzIndices) {
          const globalIdx = files.length + localIdx;
          const refinedType = await sniffMtzFile(refined[globalIdx].file);
          refined[globalIdx] = {
            ...refined[globalIdx],
            detectedType: refinedType,
          };
        }
        onChange(refined);
        setIsSniffing(false);
      }
    },
    [files, onChange]
  );

  const handleDrop = useCallback(
    (e: React.DragEvent) => {
      e.preventDefault();
      setIsDragOver(false);
      if (e.dataTransfer.files) addFiles(e.dataTransfer.files);
    },
    [addFiles]
  );

  const handleFileInput = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      if (e.target.files) addFiles(e.target.files);
      e.target.value = "";
    },
    [addFiles]
  );

  const removeFile = useCallback(
    (index: number) => {
      onChange(files.filter((_, i) => i !== index));
    },
    [files, onChange]
  );

  const openFilePicker = () => {
    const input = document.createElement("input");
    input.type = "file";
    input.multiple = true;
    input.accept = ACCEPTED_EXTENSIONS;
    input.onchange = (e) =>
      handleFileInput(e as unknown as React.ChangeEvent<HTMLInputElement>);
    input.click();
  };

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
        onClick={openFilePicker}
      >
        {files.length === 0 ? (
          <Stack spacing={1} alignItems="center">
            <UploadIcon sx={{ fontSize: 32, color: "text.secondary" }} />
            <Typography variant="body2" color="text.secondary">
              Drop files here to import after project creation
            </Typography>
            <Typography variant="caption" color="text.secondary">
              MTZ, PDB, CIF, FASTA, alignments, SCA, ligands (SMI/MOL)
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
                  {isSniffing &&
                    df.file.name.toLowerCase().endsWith(".mtz") && (
                      <CircularProgress size={12} />
                    )}
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
              onClick={openFilePicker}
            >
              + Drop or click to add more files
            </Typography>
          </Stack>
        )}
      </Paper>
    </Stack>
  );
};
