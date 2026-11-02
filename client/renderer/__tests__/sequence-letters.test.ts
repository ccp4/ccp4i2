import { describe, it, expect } from "vitest";
import {
  describeOddLetters,
  errorWindows,
  nonCanonicalLetters,
  previewWindow,
  sequenceSegments,
} from "../lib/sequence-letters";

describe("sequenceSegments", () => {
  it("splits a protein sequence into accepted and rejected runs with positions", () => {
    expect(sequenceSegments("MSTXVA", "PROTEIN")).toEqual([
      { text: "MST", start: 1, canonical: true },
      { text: "X", start: 4, canonical: false },
      { text: "VA", start: 5, canonical: true },
    ]);
  });
  it("ignores whitespace and case, and knows the nucleic alphabet", () => {
    expect(sequenceSegments("ac gu\nT", "RNA").every((s) => s.canonical)).toBe(true);
    expect(sequenceSegments("ACGTX", "DNA").filter((s) => !s.canonical)).toEqual([
      { text: "X", start: 5, canonical: false },
    ]);
    // B, Z, U and O are not in the server's protein alphabet
    expect(nonCanonicalLetters("MBZUO", "PROTEIN").map((o) => o.letter)).toEqual(["B", "Z", "U", "O"]);
  });
});

describe("describeOddLetters", () => {
  it("names letters with positions and summarises a long tail", () => {
    expect(describeOddLetters([])).toBe("");
    expect(describeOddLetters([{ letter: "X", position: 161 }])).toBe("X at 161");
    expect(
      describeOddLetters([
        { letter: "X", position: 1 }, { letter: "X", position: 2 },
        { letter: "X", position: 3 }, { letter: "X", position: 4 }, { letter: "X", position: 5 },
      ]),
    ).toBe("X at 1, X at 2, X at 3 and 2 more");
  });
});

describe("previewWindow", () => {
  const body = "A".repeat(300);
  it("shows the start of a clean sequence", () => {
    expect(previewWindow(body, "PROTEIN", 40)).toEqual({
      text: "A".repeat(40), start: 1, truncatedBefore: false, truncatedAfter: true,
    });
  });
  it("centres on the first rejected letter when it lies past the preview", () => {
    const seq = "A".repeat(160) + "X" + "A".repeat(139); // X at 161, as 1H1S's TPO
    const w = previewWindow(seq, "PROTEIN", 40);
    expect(w.text).toContain("X");
    expect(w.text).toHaveLength(40);
    expect(w.start).toBe(141);
    expect(w.truncatedBefore).toBe(true);
    expect(w.truncatedAfter).toBe(true);
  });
  it("keeps the start when the rejected letter is already in view", () => {
    const seq = "AAAXA" + "A".repeat(300);
    expect(previewWindow(seq, "PROTEIN", 40).start).toBe(1);
  });
  it("does not run past the end", () => {
    const seq = "A".repeat(295) + "X" + "AAAA";
    const w = previewWindow(seq, "PROTEIN", 40);
    expect(w.start + w.text.length - 1).toBe(300);
    expect(w.truncatedAfter).toBe(false);
  });
});

describe("errorWindows", () => {
  it("is empty for a clean sequence", () => {
    expect(errorWindows("MSTVA", "PROTEIN")).toEqual({ windows: [], more: 0 });
  });
  it("excerpts seven residues either side, clipped to the sequence", () => {
    const seq = "A".repeat(160) + "X" + "A".repeat(139);
    const { windows } = errorWindows(seq, "PROTEIN");
    expect(windows).toEqual([{ start: 154, end: 168, text: "A".repeat(7) + "X" + "A".repeat(7) }]);
    expect(errorWindows("XAAAAAAAAAAAA", "PROTEIN").windows[0]).toEqual({ start: 1, end: 8, text: "XAAAAAAA" });
  });
  it("merges neighbouring letters into one excerpt and caps the count", () => {
    const seq = "A".repeat(20) + "X" + "AAA" + "X" + "A".repeat(20);
    const { windows } = errorWindows(seq, "PROTEIN");
    expect(windows).toHaveLength(1);
    expect(windows[0]).toEqual({ start: 14, end: 32, text: seq.substring(13, 32) });
    const many = Array.from({ length: 10 }, () => "A".repeat(20) + "X").join("");
    expect(errorWindows(many, "PROTEIN", 7, 6)).toMatchObject({ more: 4 });
  });
});
