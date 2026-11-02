/**
 * Which letters of an AU-contents sequence the server accepts, and where
 * the others are, so a table row or its editor can show them rather than
 * leave the reader to find an X in 300 residues (1H1S, 2026-10-31).
 *
 * The alphabets mirror CAsuContentSeq.validity (core/CCP4ModelData.py):
 * the twenty standard amino acids for protein, C A U G T for DNA and RNA.
 * Case does not matter there; whitespace is ignored.
 */

const PROTEIN = new Set("GALMFWKQESPVICYHRNDT");
const NUCLEIC = new Set("CAUGT");

export type PolymerType = "PROTEIN" | "DNA" | "RNA" | string;

export function canonicalLetters(polymerType: PolymerType): Set<string> {
  return (polymerType || "").toUpperCase() === "PROTEIN" ? PROTEIN : NUCLEIC;
}

export interface SequenceSegment {
  text: string;
  /** 1-based position of the segment's first letter in the whitespace-free sequence. */
  start: number;
  canonical: boolean;
}

/** The sequence as runs of accepted and unaccepted letters, in order. */
export function sequenceSegments(sequence: string, polymerType: PolymerType): SequenceSegment[] {
  const seq = (sequence || "").replace(/\s/g, "");
  const ok = canonicalLetters(polymerType);
  const out: SequenceSegment[] = [];
  for (let i = 0; i < seq.length; i++) {
    const canonical = ok.has(seq[i].toUpperCase());
    const last = out[out.length - 1];
    if (last && last.canonical === canonical) last.text += seq[i];
    else out.push({ text: seq[i], start: i + 1, canonical });
  }
  return out;
}

export interface OddLetter {
  letter: string;
  /** 1-based position in the whitespace-free sequence. */
  position: number;
}

/** Every letter the server will reject, with its position. */
export function nonCanonicalLetters(sequence: string, polymerType: PolymerType): OddLetter[] {
  return sequenceSegments(sequence, polymerType)
    .filter((s) => !s.canonical)
    .flatMap((s) => Array.from(s.text, (letter, i) => ({ letter, position: s.start + i })));
}

/** "X at 161" / "X at 161, B at 200" / "X at 161 and 3 more". */
export function describeOddLetters(odd: OddLetter[], max = 3): string {
  if (odd.length === 0) return "";
  const shown = odd.slice(0, max).map((o) => `${o.letter} at ${o.position}`).join(", ");
  const rest = odd.length - max;
  return rest > 0 ? `${shown} and ${rest} more` : shown;
}

/**
 * The part of a long sequence a one-line preview shows: its start, unless
 * it holds a letter the server rejects, when the window is centred on the
 * first one so the preview shows the problem. Returns the window's text
 * and its 1-based start.
 */
export function previewWindow(
  sequence: string,
  polymerType: PolymerType,
  maxLength = 40,
): { text: string; start: number; truncatedBefore: boolean; truncatedAfter: boolean } {
  const seq = (sequence || "").replace(/\s/g, "");
  if (seq.length <= maxLength) return { text: seq, start: 1, truncatedBefore: false, truncatedAfter: false };
  const first = nonCanonicalLetters(seq, polymerType)[0];
  let begin = 0;
  if (first && first.position > maxLength) {
    begin = Math.min(first.position - 1 - Math.floor(maxLength / 2), seq.length - maxLength);
  }
  return {
    text: seq.substring(begin, begin + maxLength),
    start: begin + 1,
    truncatedBefore: begin > 0,
    truncatedAfter: begin + maxLength < seq.length,
  };
}

export interface SequenceWindow {
  /** 1-based position of the first letter shown. */
  start: number;
  /** 1-based position of the last letter shown. */
  end: number;
  text: string;
}

/**
 * Numbered excerpts of the sequence around each rejected letter, `context`
 * residues either side, neighbours merged into one excerpt. What the row
 * editor shows instead of the whole sequence: a report of where the
 * problems are, clearly not a second place to type.
 */
export function errorWindows(
  sequence: string,
  polymerType: PolymerType,
  context = 7,
  max = 6,
): { windows: SequenceWindow[]; more: number } {
  const seq = (sequence || "").replace(/\s/g, "");
  const odd = nonCanonicalLetters(seq, polymerType);
  const windows: SequenceWindow[] = [];
  for (const o of odd) {
    const start = Math.max(1, o.position - context);
    const end = Math.min(seq.length, o.position + context);
    const last = windows[windows.length - 1];
    if (last && start <= last.end + 1) {
      last.end = Math.max(last.end, end);
      last.text = seq.substring(last.start - 1, last.end);
    } else {
      windows.push({ start, end, text: seq.substring(start - 1, end) });
    }
  }
  return { windows: windows.slice(0, max), more: Math.max(0, windows.length - max) };
}
