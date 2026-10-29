/**
 * The file preview's font-size buttons step along a fixed scale and stop at
 * its ends (#790).
 */
import { describe, it, expect } from "vitest";
import { stepFontSize } from "../lib/preview-font-size";
import { withArticle } from "../lib/with-article";

describe("stepFontSize", () => {
  it("steps up and down the scale", () => {
    expect(stepFontSize(14, 1)).toBe(16);
    expect(stepFontSize(14, -1)).toBe(13);
  });

  it("stays put at either end", () => {
    expect(stepFontSize(28, 1)).toBe(28);
    expect(stepFontSize(10, -1)).toBe(10);
  });

  it("steps from the default when the current size is off the scale", () => {
    expect(stepFontSize(15, 1)).toBe(16);
  });
});

describe("withArticle (list add button, #837)", () => {
  it("picks a or an by the first letter", () => {
    expect(withArticle("map")).toBe("a map");
    expect(withArticle("atomic model")).toBe("an atomic model");
  });
});
