/**
 * Opening a pending job must not change it. The atom-selection builder pushed
 * its own string (all chains, no filters) to the server on mount, whenever it
 * differed from the job's: a cloned or reopened job lost its selection,
 * "(NUT)" or "not (HOH)", the moment its page was displayed. The push waits
 * until the user has used the builder.
 */
import { readFileSync } from "fs";
import { join } from "path";
import { describe, expect, it } from "vitest";

const source = readFileSync(
  join(__dirname, "..", "components", "task", "task-elements", "cpdbdatafile.tsx"),
  "utf8"
);

describe("cpdbdatafile atom selection", () => {
  it("writes the builder's selection only after the builder is used", () => {
    expect(source).toMatch(/if \(!builderUsed \|\| manualEdit \|\|/);
  });

  it("marks the builder used wherever a builder control hands it control", () => {
    const handsOver = source.match(/setManualEdit\(false\);/g) ?? [];
    const marks = source.match(/setManualEdit\(false\);\s*setBuilderUsed\(true\);/g) ?? [];
    expect(handsOver.length).toBeGreaterThan(0);
    expect(marks.length).toBe(handsOver.length);
  });
});
