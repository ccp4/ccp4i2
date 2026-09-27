/**
 * The arithmetic behind the warning, pinned on the CDK4 campaign that
 * prompted it: 48 datasets between 2.04 and 3.92 A plus two at 6.48 and
 * 6.71 A, processed at 6.71 A because max_shell_datasets (60) exceeds the
 * dataset count, so every crystal is a comparator and max() governs.
 */
import { describe, expect, it } from "vitest";

import { shellResolution } from "../components/task/task-interfaces/pandda-resolution-notice";

const CDK4 = [
  ...Array.from({ length: 48 }, (_, i) => 2.04 + i * 0.04),
  6.48,
  6.71,
];

describe("shellResolution", () => {
  it("returns the worst crystal when the shell is bigger than the set", () => {
    expect(shellResolution(CDK4, 60)).toBeCloseTo(6.71, 2);
  });

  it("returns the nth best once the shell is smaller than the set", () => {
    const sorted = [...CDK4].sort((a, b) => a - b);
    expect(shellResolution(CDK4, 30)).toBeCloseTo(sorted[29], 6);
  });

  it("is unaffected by the order it is given", () => {
    const shuffled = [...CDK4].reverse();
    expect(shellResolution(shuffled, 30)).toBeCloseTo(shellResolution(CDK4, 30)!, 6);
  });

  it("has nothing to say about an empty set", () => {
    expect(shellResolution([], 60)).toBeNull();
  });

  it("handles a shell larger than the set without running off the end", () => {
    expect(shellResolution([2.0, 2.5], 60)).toBeCloseTo(2.5, 6);
  });
});
