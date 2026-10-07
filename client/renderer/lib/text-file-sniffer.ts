/**
 * Content sniffing for text-based crystallographic files whose extension is
 * ambiguous or compound.
 *
 * Two extensions cannot be classified by name alone:
 *
 * - `.cif` is used by three unrelated datatypes — coordinates (mmCIF),
 *   merged reflection data (structure factors, as downloaded from the PDB),
 *   and monomer/restraint dictionaries. Routing all three to a coordinate
 *   import produces a job that fails or, worse, quietly does the wrong thing.
 * - `.xml` covers CCP4i2's own ASU-contents files, which carry the compound
 *   extension `.asu.xml` but are otherwise indistinguishable by suffix.
 *
 * Only a bounded prefix of the file is read: the discriminating category
 * markers appear in the header, and a structure-factor CIF can be tens of
 * megabytes of reflection data we have no reason to pull into the renderer.
 */

/** How much of the file to inspect. Generous enough for a long mmCIF preamble. */
const SNIFF_BYTES = 64 * 1024;

/** What a `.cif` file turned out to contain. */
export type CifFlavour = "coordinates" | "reflections" | "dictionary" | "unknown";

/** What an `.xml` file turned out to contain. */
export type XmlFlavour = "asu_contents" | "unknown";

/** Read the leading `SNIFF_BYTES` of a file as text. */
async function readPrefix(file: File): Promise<string> {
  const slice = file.slice(0, SNIFF_BYTES);
  const buffer = await slice.arrayBuffer();
  // CIF and CCP4i2 XML are both ASCII in practice; a lone truncated multi-byte
  // character at the cut is harmless because we only match ASCII markers.
  return new TextDecoder("utf-8", { fatal: false }).decode(buffer);
}

/**
 * Classify a CIF file by the mmCIF categories it declares.
 *
 * Order matters. A refined coordinate file commonly carries `_chem_comp`
 * (describing the entities it contains) without being a dictionary, so the
 * coordinate test must come first. Conversely a restraint dictionary has no
 * `_atom_site.` at all — its atoms live in `_chem_comp_atom.` — so the two
 * do not overlap once checked in this order.
 */
export function classifyCifText(text: string): CifFlavour {
  // Category markers are line-initial in every CIF writer we care about.
  // `marker` is a plain literal; escape it here rather than at each call site.
  const has = (marker: string) =>
    new RegExp(
      `^\\s*${marker.replace(/[.\\^$*+?()[\]{}|]/g, "\\$&")}`,
      "m"
    ).test(text);

  if (has("_atom_site.")) return "coordinates";
  if (has("_refln.") || has("_diffrn_refln.")) return "reflections";
  if (has("_chem_comp_atom.") || has("data_comp_")) return "dictionary";
  return "unknown";
}

/** Classify an XML file by its CCP4i2 header `<function>` element. */
export function classifyXmlText(text: string): XmlFlavour {
  if (/<function>\s*ASUCONTENT\s*<\/function>/i.test(text)) return "asu_contents";
  return "unknown";
}

/** Read and classify a `.cif` file. Unreadable files fall back to "unknown". */
export async function sniffCifFile(file: File): Promise<CifFlavour> {
  try {
    return classifyCifText(await readPrefix(file));
  } catch {
    return "unknown";
  }
}

/** Read and classify an `.xml` file. Unreadable files fall back to "unknown". */
export async function sniffXmlFile(file: File): Promise<XmlFlavour> {
  try {
    return classifyXmlText(await readPrefix(file));
  } catch {
    return "unknown";
  }
}

/** How much of a file `looksBinary` inspects. */
const BINARY_SNIFF_BYTES = 8 * 1024;

/**
 * Whether file content is binary rather than text, judged from its leading
 * bytes rather than its name: a file in a job directory called `.log`,
 * `.out` or `.txt` is not guaranteed to be text, and an MTZ or a map shown
 * as text is a screen of mojibake.
 *
 * A NUL byte is decisive (text files never carry one). Otherwise a prefix
 * more than a tenth control characters, other than whitespace, backspace
 * and the escape that starts a terminal colour code, is taken as binary.
 * Bytes above 0x7f are not counted, so a log in Latin-1 or UTF-8 reads as
 * text either way.
 */
export function looksBinary(content: ArrayBuffer | Uint8Array): boolean {
  const bytes =
    content instanceof Uint8Array
      ? content.subarray(0, BINARY_SNIFF_BYTES)
      : new Uint8Array(content, 0, Math.min(content.byteLength, BINARY_SNIFF_BYTES));
  if (bytes.length === 0) return false;
  let control = 0;
  for (const byte of bytes) {
    if (byte === 0) return true;
    const allowed = (byte >= 8 && byte <= 13) || byte === 27;
    if ((byte < 32 && !allowed) || byte === 127) control++;
  }
  return control / bytes.length > 0.1;
}
