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
 * - Byte 20+: Reflection data
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
        // Missing value marker - skip
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
    // 64-bit offset at bytes 12-19
    // For simplicity, assume files aren't that large and use lower 32 bits
    // (JavaScript numbers can safely represent integers up to 2^53)
    const low = readInt32(view, 12, swapBytes);
    const high = readInt32(view, 16, swapBytes);
    headerOffset = (high * 0x100000000) + (low >>> 0);
  } else {
    headerOffset = headerOffset32;
  }

  // Header offset is in "words" (4-byte units), 1-indexed (Fortran convention)
  // The header offset value N means header starts at word N, which is byte (N-1)*4
  // because Fortran arrays are 1-indexed
  const headerByteOffset = (headerOffset - 1) * 4;

  if (headerByteOffset >= data.byteLength || headerByteOffset < 20) {
    throw new MtzParseError(
      `Invalid header offset: ${headerOffset} words (${headerByteOffset} bytes) ` +
      `in file of ${data.byteLength} bytes`
    );
  }

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
