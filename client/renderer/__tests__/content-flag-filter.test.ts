import { describe, expect, it } from "vitest";
import { canConvertToRequired } from "../lib/content-flags";

const OBSERVED = "application/CCP4-mtz-observed";
const COORDINATES = "chemical/x-pdb";

describe("content-flag filter in the file pickers", () => {
  it("converts observed data along the obs conversion table", () => {
    expect(canConvertToRequired(1, [4], OBSERVED)).toBe(true); // IPAIR -> FMEAN
    expect(canConvertToRequired(4, [1], OBSERVED)).toBe(false); // FMEAN -/-> IPAIR
  });

  it("requires an exact flag for any other class", () => {
    // Coordinates: 1 is PDB format, 2 is mmCIF. A PDB-format model must not
    // satisfy a field that asks for mmCIF.
    expect(canConvertToRequired(1, [2], COORDINATES)).toBe(false);
    expect(canConvertToRequired(2, [2], COORDINATES)).toBe(true);
    expect(canConvertToRequired(1, [2], undefined)).toBe(false);
  });

  it("rejects a file with no content flag", () => {
    expect(canConvertToRequired(null, [2], COORDINATES)).toBe(false);
  });
});
