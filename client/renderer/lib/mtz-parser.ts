/**
 * Pure TypeScript MTZ file header parser.
 *
 * Parses MTZ reflection file headers to extract column labels and types
 * without requiring any WASM modules (coot/gemmi).
 *
 * MTZ file format reference:
 * https://www.ccp4.ac.uk/html/mtzformat.html
 * https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/ccp4/html/mtzformat.html
 *
 * File structure:
 * - Bytes 0-3: "MTZ " identifier
 * - Bytes 4-7: 32-bit integer header offset (or -1 for 64-bit offset)
 * - Bytes 8-11: Machine stamp (byte order encoding)
 * - Bytes 12-19: 64-bit header offset (if 32-bit offset is -1)
 * - Byte 80+: Reflection data (20-word preamble × 4 bytes = 80)
 * - Header offset+: ASCII header records (80 chars each)
 *
 * Machine stamp encoding (byte 8, high nibble):
 * - 0x10: Big-endian IEEE
 * - 0x40: Little-endian IEEE
 */

/**
 * Column information extracted from MTZ header
 */
export interface MtzColumn {
  /** Column label (e.g., "FP", "SIGFP", "FREE") */
  label: string;
  /** Column type code (e.g., "F", "Q", "I", "H") */
  type: string;
  /** Minimum value in column */
  minValue: number;
  /** Maximum value in column */
  maxValue: number;
  /** Dataset ID this column belongs to */
  datasetId: number;
}

/**
 * Dataset information extracted from MTZ header
 */
export interface MtzDataset {
  id: number;
  projectName: string;
  crystalName: string;
  datasetName: string;
  cell?: {
    a: number;
    b: number;
    c: number;
    alpha: number;
    beta: number;
    gamma: number;
  };
  wavelength?: number;
}

/**
 * Complete MTZ header information
 */
export interface MtzHeader {
  /** MTZ format version */
  version: string;
  /** File title */
  title: string;
  /** Number of columns */
  nColumns: number;
  /** Number of reflections */
  nReflections: number;
  /** Number of batches (0 for merged data) */
  nBatches: number;
  /** Overall cell parameters */
  cell?: {
    a: number;
    b: number;
    c: number;
    alpha: number;
    beta: number;
    gamma: number;
  };
  /** Space group name */
  spaceGroup?: string;
  /** Space group number */
  spaceGroupNumber?: number;
  /** Resolution limits [min, max] in Angstroms */
  resolution?: [number, number];
  /** Column definitions */
  columns: MtzColumn[];
  /** Dataset definitions */
  datasets: MtzDataset[];
  /** Whether data is merged (nBatches === 0) */
  isMerged: boolean;
  /** Missing value marker from VALM record (NaN if not specified) */
  missingValue: number;
}

/**
 * Reflection data extracted from an MTZ file
 */
export interface MtzReflectionData {
  reflections: { h: number; k: number; l: number; intensity?: number }[];
  hRange: [number, number];
  kRange: [number, number];
  lRange: [number, number];
  intensityRange?: [number, number];
  intensityLabel?: string;
}

/**
 * Error thrown when MTZ parsing fails
 */
export class MtzParseError extends Error {
  constructor(message: string) {
    super(message);
    this.name = "MtzParseError";
  }
}

/**
 * Detect if the current system is little-endian
 */
function isLittleEndian(): boolean {
  const buffer = new ArrayBuffer(2);
  new DataView(buffer).setInt16(0, 256, true); // little-endian
  return new Int16Array(buffer)[0] === 256;
}

/**
 * Swap byte order of a 32-bit value
 */
function swap32(val: number): number {
  return (
    ((val & 0xff) << 24) |
    ((val & 0xff00) << 8) |
    ((val & 0xff0000) >> 8) |
    ((val >> 24) & 0xff)
  );
}

/**
 * Read a 32-bit integer from DataView, handling byte order
 */
function readInt32(view: DataView, offset: number, swapBytes: boolean): number {
  const val = view.getInt32(offset, true); // read as little-endian
  return swapBytes ? swap32(val) : val;
}

/**
 * Read a 32-bit float from DataView, handling byte order
 */
function readFloat32(view: DataView, offset: number, swapBytes: boolean): number {
  if (swapBytes) {
    // Need to swap bytes before interpreting as float
    const bytes = new Uint8Array(4);
    for (let i = 0; i < 4; i++) {
      bytes[i] = view.getUint8(offset + 3 - i);
    }
    return new DataView(bytes.buffer).getFloat32(0, true);
  }
  return view.getFloat32(offset, true);
}

/**
 * Parse the ASCII header section of an MTZ file.
 *
 * MTZ header records are 80 characters each, written in Fortran 20A4 format.
 * Records start with a keyword (e.g., VERS, TITL, NCOL, CELL) followed by data.
 */
function parseHeaderRecords(headerText: string): MtzHeader {
  const header: MtzHeader = {
    version: "",
    title: "",
    nColumns: 0,
    nReflections: 0,
    nBatches: 0,
    columns: [],
    datasets: [],
    isMerged: true,
    missingValue: NaN,
  };

  // Split into 80-character records
  const lines: string[] = [];
  for (let i = 0; i < headerText.length; i += 80) {
    lines.push(headerText.slice(i, i + 80).trimEnd());
  }

  // Temporary storage for dataset info
  const datasetMap = new Map<number, Partial<MtzDataset>>();
  let currentDatasetId = 0;

  for (const line of lines) {
    if (!line || line.startsWith("END")) break;

    // Keywords can be 4-6 characters, followed by space and arguments
    // Common keywords: VERS, TITLE, NCOL, CELL, SORT, SYMINF, SYMM, RESO, VALM, COLUMN, etc.
    const firstSpace = line.indexOf(" ");
    const keyword = firstSpace > 0 ? line.slice(0, firstSpace) : line.slice(0, 6);
    const args = firstSpace > 0 ? line.slice(firstSpace).trim() : "";

    switch (keyword) {
      case "VERS":
        header.version = args;
        break;

      case "TITLE":
        header.title = args;
        break;

      case "NCOL":
        // NCOL ncols nrefs nbatches
        const ncolParts = args.split(/\s+/);
        header.nColumns = parseInt(ncolParts[0], 10) || 0;
        header.nReflections = parseInt(ncolParts[1], 10) || 0;
        header.nBatches = parseInt(ncolParts[2], 10) || 0;
        header.isMerged = header.nBatches === 0;
        break;

      case "CELL":
        // CELL a b c alpha beta gamma
        const cellParts = args.split(/\s+/).map(parseFloat);
        if (cellParts.length >= 6) {
          header.cell = {
            a: cellParts[0],
            b: cellParts[1],
            c: cellParts[2],
            alpha: cellParts[3],
            beta: cellParts[4],
            gamma: cellParts[5],
          };
        }
        break;

      case "SORT":
        // Sort order - skip
        break;

      case "SYMINF":
        // SYMINF nops nprim lattype spgnum spgname pgname
        const symiParts = args.split(/\s+/);
        if (symiParts.length >= 4) {
          header.spaceGroupNumber = parseInt(symiParts[3], 10);
        }
        // Extract space group name (may be quoted or have apostrophe)
        const spgMatch = args.match(/'([^']+)'|(\S+)$/);
        if (spgMatch) {
          header.spaceGroup = spgMatch[1] || spgMatch[2];
        }
        break;

      case "SYMM":
        // Symmetry operators - skip for now
        break;

      case "RESO":
        // RESO minres maxres (stored as 1/d^2)
        const resoParts = args.split(/\s+/).map(parseFloat);
        if (resoParts.length >= 2) {
          // Convert from 1/d^2 to d (Angstroms)
          const d1 = resoParts[0] > 0 ? Math.sqrt(1 / resoParts[0]) : 999;
          const d2 = resoParts[1] > 0 ? Math.sqrt(1 / resoParts[1]) : 999;
          header.resolution = [Math.max(d1, d2), Math.min(d1, d2)];
        }
        break;

      case "VALM":
        // Missing value marker (MNF)
        const valmStr = args.trim().toUpperCase();
        if (valmStr === "NAN" || valmStr === "NAN(QUIET)") {
          header.missingValue = NaN;
        } else {
          const val = parseFloat(args);
          if (!isNaN(val)) header.missingValue = val;
        }
        break;

      case "COLUMN":
        // COLUMN label type min max dataset_id
        const colMatch = args.match(/^(\S+)\s+(\S)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+(\d+)/);
        if (colMatch) {
          header.columns.push({
            label: colMatch[1],
            type: colMatch[2],
            minValue: parseFloat(colMatch[3]),
            maxValue: parseFloat(colMatch[4]),
            datasetId: parseInt(colMatch[5], 10),
          });
        }
        break;

      case "NDIF":
        // Number of datasets - just informational
        break;

      case "PROJECT":
        // PROJECT dataset_id project_name
        const projMatch = args.match(/^(\d+)\s+(.+)/);
        if (projMatch) {
          const id = parseInt(projMatch[1], 10);
          if (!datasetMap.has(id)) datasetMap.set(id, { id });
          datasetMap.get(id)!.projectName = projMatch[2].trim();
        }
        break;

      case "CRYSTAL":
        // CRYSTAL dataset_id crystal_name
        const crysMatch = args.match(/^(\d+)\s+(.+)/);
        if (crysMatch) {
          const id = parseInt(crysMatch[1], 10);
          if (!datasetMap.has(id)) datasetMap.set(id, { id });
          datasetMap.get(id)!.crystalName = crysMatch[2].trim();
        }
        break;

      case "DATASET":
        // DATASET dataset_id dataset_name
        const dataMatch = args.match(/^(\d+)\s+(.+)/);
        if (dataMatch) {
          const id = parseInt(dataMatch[1], 10);
          if (!datasetMap.has(id)) datasetMap.set(id, { id });
          datasetMap.get(id)!.datasetName = dataMatch[2].trim();
          currentDatasetId = id;
        }
        break;

      case "DCELL":
        // DCELL dataset_id a b c alpha beta gamma
        const dcellMatch = args.match(/^(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)/);
        if (dcellMatch) {
          const id = parseInt(dcellMatch[1], 10);
          if (!datasetMap.has(id)) datasetMap.set(id, { id });
          datasetMap.get(id)!.cell = {
            a: parseFloat(dcellMatch[2]),
            b: parseFloat(dcellMatch[3]),
            c: parseFloat(dcellMatch[4]),
            alpha: parseFloat(dcellMatch[5]),
            beta: parseFloat(dcellMatch[6]),
            gamma: parseFloat(dcellMatch[7]),
          };
        }
        break;

      case "DWAVEL":
        // DWAVEL dataset_id wavelength
        const dwavMatch = args.match(/^(\d+)\s+([\d.]+)/);
        if (dwavMatch) {
          const id = parseInt(dwavMatch[1], 10);
          if (!datasetMap.has(id)) datasetMap.set(id, { id });
          datasetMap.get(id)!.wavelength = parseFloat(dwavMatch[2]);
        }
        break;

      case "BATCH":
        // Batch information - skip for column parsing
        break;

      default:
        // Unknown or extension keywords - skip
        break;
    }
  }

  // Convert dataset map to array
  header.datasets = Array.from(datasetMap.values()).map((d) => ({
    id: d.id ?? 0,
    projectName: d.projectName ?? "",
    crystalName: d.crystalName ?? "",
    datasetName: d.datasetName ?? "",
    cell: d.cell,
    wavelength: d.wavelength,
  }));

  return header;
}

/**
 * Internal: validate MTZ magic number and extract byte-order and header offset info.
 */
function getMtzFileInfo(data: ArrayBuffer): {
  view: DataView;
  swapBytes: boolean;
  headerByteOffset: number;
} {
  if (data.byteLength < 20) {
    throw new MtzParseError("File too small to be an MTZ file");
  }

  const view = new DataView(data);

  // Check magic number "MTZ "
  const magic = String.fromCharCode(
    view.getUint8(0),
    view.getUint8(1),
    view.getUint8(2),
    view.getUint8(3)
  );
  if (magic !== "MTZ ") {
    throw new MtzParseError(`Not an MTZ file - starts with "${magic}" instead of "MTZ "`);
  }

  // Detect byte order from machine stamp (byte 8, high nibble)
  // 0x10 = big-endian, 0x40 = little-endian
  const machineStamp = view.getUint8(8);
  const fileIsLittleEndian = (machineStamp & 0xf0) === 0x40;
  const systemIsLittleEndian = isLittleEndian();
  const swapBytes = fileIsLittleEndian !== systemIsLittleEndian;

  // Read header offset (32-bit at bytes 4-7, or 64-bit at 12-19 if 32-bit is -1)
  let headerOffset: number;
  const headerOffset32 = readInt32(view, 4, swapBytes);

  if (headerOffset32 === -1) {
    const low = readInt32(view, 12, swapBytes);
    const high = readInt32(view, 16, swapBytes);
    headerOffset = (high * 0x100000000) + (low >>> 0);
  } else {
    headerOffset = headerOffset32;
  }

  // Header offset is in "words" (4-byte units), 1-indexed (Fortran convention)
  const headerByteOffset = (headerOffset - 1) * 4;

  if (headerByteOffset >= data.byteLength || headerByteOffset < 20) {
    throw new MtzParseError(
      `Invalid header offset: ${headerOffset} words (${headerByteOffset} bytes) ` +
      `in file of ${data.byteLength} bytes`
    );
  }

  return { view, swapBytes, headerByteOffset };
}

/**
 * Parse an MTZ file and extract header information.
 *
 * This is a pure TypeScript implementation that doesn't require
 * any WASM modules (coot/gemmi). It handles both big-endian and
 * little-endian MTZ files.
 *
 * @param data - The MTZ file as an ArrayBuffer
 * @returns Parsed MTZ header information
 * @throws MtzParseError if the file is not a valid MTZ file
 */
export function parseMtzHeader(data: ArrayBuffer): MtzHeader {
  const { headerByteOffset } = getMtzFileInfo(data);

  // Read header as ASCII text
  const headerBytes = new Uint8Array(data, headerByteOffset);
  const decoder = new TextDecoder("ascii");
  const headerText = decoder.decode(headerBytes);

  return parseHeaderRecords(headerText);
}

/**
 * Parse an MTZ file from a File object.
 *
 * @param file - The File object to parse
 * @returns Promise resolving to parsed MTZ header information
 */
export async function parseMtzFile(file: File): Promise<MtzHeader> {
  const buffer = await file.arrayBuffer();
  return parseMtzHeader(buffer);
}

/**
 * Get column names organized by type, suitable for column selection UI.
 * Returns a map from column label to column type code.
 *
 * This matches the format expected by the existing mtz-column-dialog.tsx
 *
 * @param header - Parsed MTZ header
 * @returns Map of column label to type code
 */
export function getColumnNamesByType(header: MtzHeader): Record<string, string> {
  const result: Record<string, string> = {};
  for (const col of header.columns) {
    result[col.label] = col.type;
  }
  return result;
}

/**
 * Get column information with dataset names for UI display.
 *
 * @param header - Parsed MTZ header
 * @returns Array of columns with dataset information
 */
export function getColumnsWithDatasets(header: MtzHeader): Array<{
  label: string;
  type: string;
  dataset: string;
  groupIndex: number;
}> {
  // Build dataset name lookup
  const datasetNames = new Map<number, string>();
  for (const ds of header.datasets) {
    datasetNames.set(ds.id, ds.datasetName || ds.crystalName || `Dataset ${ds.id}`);
  }

  return header.columns.map((col, idx) => ({
    label: col.label,
    type: col.type,
    dataset: datasetNames.get(col.datasetId) || "",
    groupIndex: idx,
  }));
}

/**
 * Parse reflection data from an MTZ file, extracting h, k, l indices
 * and optionally an intensity/amplitude column.
 *
 * Reflection data starts at byte 20 in the file. Each reflection row
 * contains nColumns float32 values. We only read the H, K, L columns
 * (type "H") and an optional intensity column for efficiency.
 *
 * @param data - The MTZ file as an ArrayBuffer
 * @param header - Previously parsed header (from parseMtzHeader)
 * @param intensityColumn - Optional column label for intensity coloring
 * @returns Extracted reflection data with h,k,l indices and ranges
 */
export function parseMtzReflections(
  data: ArrayBuffer,
  header: MtzHeader,
  intensityColumn?: string
): MtzReflectionData {
  const { view, swapBytes } = getMtzFileInfo(data);

  // Find H, K, L column indices (type "H", always first 3 of this type)
  const hklIndices: number[] = [];
  for (let i = 0; i < header.columns.length; i++) {
    if (header.columns[i].type === "H") {
      hklIndices.push(i);
    }
  }
  if (hklIndices.length < 3) {
    throw new MtzParseError(
      `Expected at least 3 H-type columns (Miller indices), found ${hklIndices.length}`
    );
  }
  const hIdx = hklIndices[0];
  const kIdx = hklIndices[1];
  const lIdx = hklIndices[2];

  // Find optional intensity column
  let intensityIdx = -1;
  let intensityLabel: string | undefined;
  if (intensityColumn) {
    const idx = header.columns.findIndex((c) => c.label === intensityColumn);
    if (idx >= 0) {
      intensityIdx = idx;
      intensityLabel = intensityColumn;
    }
  }
  // Auto-detect: find the best column for coloring, in priority order:
  // J (mean intensity), F (mean amplitude), K (anom intensity), G (anom amplitude), E (normalised F)
  if (intensityIdx < 0) {
    const priorities = ["J", "F", "K", "G", "E"];
    for (const type of priorities) {
      const idx = header.columns.findIndex((c) => c.type === type);
      if (idx >= 0) {
        intensityIdx = idx;
        intensityLabel = header.columns[idx].label;
        break;
      }
    }
  }

  const dataStart = 80; // Reflection data starts at byte 80 (20 words × 4 bytes)
  const rowBytes = header.nColumns * 4;
  const nReflections = header.nReflections;

  // Validate data extent
  const dataEnd = dataStart + nReflections * rowBytes;
  if (dataEnd > data.byteLength) {
    throw new MtzParseError(
      `Reflection data extends beyond file: needs ${dataEnd} bytes, file is ${data.byteLength}`
    );
  }

  const reflections: { h: number; k: number; l: number; intensity?: number }[] = [];
  let hMin = Infinity, hMax = -Infinity;
  let kMin = Infinity, kMax = -Infinity;
  let lMin = Infinity, lMax = -Infinity;
  let intMin = Infinity, intMax = -Infinity;

  // Missing value check: uses VALM from header (NaN or a specific sentinel value)
  const mnf = header.missingValue;
  const isMissing = isNaN(mnf)
    ? (v: number) => isNaN(v)
    : (v: number) => isNaN(v) || v === mnf;

  for (let i = 0; i < nReflections; i++) {
    const rowStart = dataStart + i * rowBytes;

    const hRaw = readFloat32(view, rowStart + hIdx * 4, swapBytes);
    const kRaw = readFloat32(view, rowStart + kIdx * 4, swapBytes);
    const lRaw = readFloat32(view, rowStart + lIdx * 4, swapBytes);

    // Skip reflections with missing H, K, or L
    if (isMissing(hRaw) || isMissing(kRaw) || isMissing(lRaw)) {
      continue;
    }

    const h = Math.round(hRaw);
    const k = Math.round(kRaw);
    const l = Math.round(lRaw);

    // Read intensity column; undefined if missing
    let intensity: number | undefined;
    if (intensityIdx >= 0) {
      const val = readFloat32(view, rowStart + intensityIdx * 4, swapBytes);
      if (!isMissing(val)) {
        intensity = val;
        if (val < intMin) intMin = val;
        if (val > intMax) intMax = val;
      }
    }

    reflections.push({ h, k, l, intensity });

    if (h < hMin) hMin = h;
    if (h > hMax) hMax = h;
    if (k < kMin) kMin = k;
    if (k > kMax) kMax = k;
    if (l < lMin) lMin = l;
    if (l > lMax) lMax = l;
  }

  return {
    reflections,
    hRange: reflections.length > 0 ? [hMin, hMax] : [0, 0],
    kRange: reflections.length > 0 ? [kMin, kMax] : [0, 0],
    lRange: reflections.length > 0 ? [lMin, lMax] : [0, 0],
    intensityRange:
      intensityIdx >= 0 && intMin <= intMax ? [intMin, intMax] : undefined,
    intensityLabel,
  };
}
