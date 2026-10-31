import { describe, it, expect } from "vitest";
import {
  bestMatthewsRow,
  describeCopies,
  namesOfRows,
  commonFactor,
  copiesOfRows,
  scaledCopies,
} from "../lib/asu-copies";

const row = (nCopies: unknown) => ({ _value: { nCopies: { _value: nCopies } } });

describe("copiesOfRows", () => {
  it("reads the copies and treats unset or bad values as one", () => {
    expect(copiesOfRows([row(2), row(null), row("3"), row(0), row(-1)])).toEqual([2, 1, 3, 1, 1]);
  });
});

describe("commonFactor", () => {
  it("is the gcd of the copies, and 1 for an empty table", () => {
    expect(commonFactor([2, 2])).toBe(2);
    expect(commonFactor([4, 6])).toBe(2);
    expect(commonFactor([2, 1])).toBe(1);
    expect(commonFactor([3])).toBe(3);
    expect(commonFactor([])).toBe(1);
  });
});

describe("scaledCopies", () => {
  it("multiplies the smallest whole set, so the table's ratios survive", () => {
    // 2A+2B is analysed as multiples of 1A+1B, so its current row is 2x
    expect(scaledCopies([2, 2], 2)).toEqual([2, 2]); // the current row: no change
    expect(scaledCopies([2, 2], 4)).toEqual([4, 4]);
    expect(scaledCopies([2, 2], 1)).toEqual([1, 1]); // back to one set
    expect(scaledCopies([2, 1], 2)).toEqual([4, 2]); // ratio kept
    expect(scaledCopies([1], 3)).toEqual([3]);
  });
});

describe("bestMatthewsRow", () => {
  it("picks the most probable row and null for none", () => {
    const rows = [
      { nmol_in_asu: 1, percent_solvent: 70, prob_matth: 0.1 },
      { nmol_in_asu: 2, percent_solvent: 40, prob_matth: 0.8 },
      { nmol_in_asu: 3, percent_solvent: 10, prob_matth: 0.1 },
    ];
    expect(bestMatthewsRow(rows)?.nmol_in_asu).toBe(2);
    expect(bestMatthewsRow([])).toBeNull();
    expect(bestMatthewsRow(undefined)).toBeNull();
  });
});

describe("namesOfRows", () => {
  it("reads the names and labels an unnamed row by its position", () => {
    const named = (name: unknown) => ({ _value: { name: { _value: name } } });
    expect(namesOfRows([named("Chain_A_C"), named(""), named(undefined)])).toEqual([
      "Chain_A_C",
      "Sequence 2",
      "Sequence 3",
    ]);
  });
});

describe("describeCopies", () => {
  it("spells out the copies each entry would hold for a multiple", () => {
    // 2A+2B: the smallest set is 1A+1B, so row 3 means 3 of each
    expect(describeCopies(["Chain_A_C", "Chain_B_D"], [2, 2], 3)).toBe(
      "Chain_A_C × 3, Chain_B_D × 3"
    );
    // 2A+1B keeps its ratio
    expect(describeCopies(["A", "B"], [2, 1], 2)).toBe("A × 4, B × 2");
    expect(describeCopies(["A"], [1], 1)).toBe("A × 1");
  });
});
