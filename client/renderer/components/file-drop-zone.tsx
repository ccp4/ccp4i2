"use client";

import { useCallback, useState } from "react";
import {
  Alert,
  Box,
  Chip,
  CircularProgress,
  IconButton,
  Stack,
  Typography,
} from "@mui/material";
import {
  CloudUpload as UploadIcon,
  Clear as ClearIcon,
} from "@mui/icons-material";
import { parseMtzFile } from "../lib/mtz-parser";
import { sniffCifFile, sniffXmlFile } from "../lib/text-file-sniffer";
import { DropZone } from "./common/drop-zone";

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
  | "ligand"            // .smi, .mol, .mol2, .sdf -> LidiaAcedrgNew (auto-run)
  | "dictionary"        // restraint/monomer .cif -> ImportDictionary (auto-run)
  | "mmcif_reflections" // structure-factor .cif -> cif2mtz (no auto-run)
  | "asu_contents"      // .asu.xml -> ImportAsuContent (auto-run)
  | "map"               // .map, .mrc -> ImportMap (auto-run)
  | "tls"               // .tls -> ProvideTLS (auto-run)
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
  dictionary: "ImportDictionary",
  mmcif_reflections: "cif2mtz",
  asu_contents: "ImportAsuContent",
  map: "ImportMap",
  tls: "ProvideTLS",
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
  dictionary: "inputData.DICTIN",
  mmcif_reflections: "inputData.HKLIN",
  asu_contents: "inputData.ASUIN",
  map: "inputData.MAPIN",
  tls: "inputData.TLSIN",
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
  dictionary: true,
  mmcif_reflections: false, // User picks which data block / columns to convert
  asu_contents: true,
  map: true,
  tls: true,
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
  dictionary: "Restraint dictionary",
  mmcif_reflections: "Reflections (mmCIF)",
  asu_contents: "ASU contents",
  map: "Map",
  tls: "TLS definitions",
  unknown: "Unrecognised",
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
  dictionary: "secondary",
  mmcif_reflections: "primary",
  asu_contents: "info",
  map: "primary",
  tls: "default",
  unknown: "error",
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
  ".msf": "alignment",
  ".smi": "ligand",
  ".mol": "ligand",
  ".mol2": "ligand",
  ".sdf": "ligand",
  ".map": "map",
  ".mrc": "map",
  ".tls": "tls",
};

/**
 * Detect file type from the name alone.
 *
 * This is the provisional answer shown the instant a file is dropped. For
 * the ambiguous extensions it is only a guess; `sniffFile` reads the file
 * and corrects it. See `needsSniffing` for which those are.
 */
function detectFileTypeSync(file: File): DetectedFileType {
  const name = file.name.toLowerCase();
  const ext = "." + name.split(".").pop();

  // Compound extensions must be tested before the single-suffix ones:
  // splitting on "." would otherwise reduce "gamma.asu.xml" to ".xml".
  if (name.endsWith(".asu.xml")) return "asu_contents";

  if (ext in UNMERGED_EXTENSIONS) return "unmerged";
  if (ext in EXTENSION_MAP) return EXTENSION_MAP[ext];

  // .cif covers coordinates, structure factors and restraint dictionaries.
  // Guess coordinates now and correct it in sniffFile below.
  if (ext === ".cif" || ext === ".mmcif") return "mmcif_coords";

  // A bare .xml may still be ASU contents written under another name;
  // sniffFile decides. Anything else stays unrecognised.
  if (ext === ".xml") return "unknown";

  // .mtz needs header sniffing — mark as reflections initially,
  // will be refined asynchronously
  if (ext === ".mtz") return "reflections";

  return "unknown";
}

/** Extensions whose true type can only be settled by reading the file. */
function needsSniffing(file: File): boolean {
  const name = file.name.toLowerCase();
  // ".asu.xml" already names its type unambiguously; only a bare ".xml"
  // has to be opened to find out whether it is ASU contents.
  if (name.endsWith(".asu.xml")) return false;
  return (
    name.endsWith(".mtz") ||
    name.endsWith(".cif") ||
    name.endsWith(".mmcif") ||
    name.endsWith(".xml")
  );
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

/**
 * Refine a provisional classification by reading the file's contents.
 * Returns the given fallback when the content is not recognised, so a
 * sniff that tells us nothing never downgrades a confident guess.
 */
async function sniffFile(
  file: File,
  provisional: DetectedFileType
): Promise<DetectedFileType> {
  const name = file.name.toLowerCase();

  if (name.endsWith(".mtz")) return sniffMtzFile(file);

  if (name.endsWith(".cif") || name.endsWith(".mmcif")) {
    switch (await sniffCifFile(file)) {
      case "coordinates":
        return "mmcif_coords";
      case "reflections":
        return "mmcif_reflections";
      case "dictionary":
        return "dictionary";
      default:
        return provisional;
    }
  }

  if (name.endsWith(".xml")) {
    return (await sniffXmlFile(file)) === "asu_contents"
      ? "asu_contents"
      : provisional;
  }

  return provisional;
}

const ACCEPTED_EXTENSIONS =
  ".mtz,.pdb,.ent,.cif,.mmcif,.xml,.fasta,.fa,.seq,.pir,.aln,.clw,.sto,.stk," +
  ".phylip,.phy,.msf,.sca,.hkl,.refl,.smi,.mol,.mol2,.sdf,.map,.mrc,.tls";

interface FileDropZoneProps {
  files: DroppedFile[];
  onChange: (files: DroppedFile[]) => void;
}

export const FileDropZone: React.FC<FileDropZoneProps> = ({
  files,
  onChange,
}) => {
  const [isSniffing, setIsSniffing] = useState(false);

  const unrecognisedCount = files.filter(
    (df) => TASK_FOR_TYPE[df.detectedType] === null
  ).length;

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

      // Read the ambiguous ones (MTZ merged/unmerged, .cif flavour,
      // .xml ASU contents) to settle their type.
      const sniffIndices = newFiles
        .map((df, i) => (needsSniffing(df.file) ? i : -1))
        .filter((i) => i >= 0);

      if (sniffIndices.length > 0) {
        setIsSniffing(true);
        const refined = [...updatedFiles];
        for (const localIdx of sniffIndices) {
          const globalIdx = files.length + localIdx;
          const current = refined[globalIdx];
          refined[globalIdx] = {
            ...current,
            detectedType: await sniffFile(current.file, current.detectedType),
          };
        }
        onChange(refined);
        setIsSniffing(false);
      }
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
      <DropZone
        onFilesSelected={addFiles}
        accept={ACCEPTED_EXTENSIONS}
        multiple
        sx={{ p: files.length > 0 ? 1.5 : 4 }}
      >
        {files.length === 0 ? (
          <Stack spacing={1} alignItems="center">
            <UploadIcon sx={{ fontSize: 32, color: "text.secondary" }} />
            <Typography variant="body2" color="text.secondary">
              Drop files here to import after project creation
            </Typography>
            <Typography variant="caption" color="text.secondary">
              MTZ, PDB, CIF, ASU contents, maps, TLS, FASTA, alignments,
              SCA, ligands (SMI/MOL/SDF)
            </Typography>
          </Stack>
        ) : (
          <Stack spacing={0.5}>
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
                  {isSniffing && needsSniffing(df.file) && (
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
            <Typography variant="caption" color="text.secondary" sx={{ mt: 0.5 }}>
              + Drop or click to add more files
            </Typography>
            {unrecognisedCount > 0 && !isSniffing && (
              // The enclosing DropZone opens a file picker on click, which
              // is unhelpful while reading a warning: keep the click here.
              <Alert
                severity="warning"
                sx={{ mt: 1, textAlign: "left" }}
                onClick={(e) => e.stopPropagation()}
              >
                {unrecognisedCount === 1
                  ? "1 file was not recognised and will not be imported."
                  : `${unrecognisedCount} files were not recognised and ` +
                    "will not be imported."}{" "}
                Remove them, or create the project and import them by hand
                afterwards.
              </Alert>
            )}
          </Stack>
        )}
      </DropZone>
    </Stack>
  );
};
