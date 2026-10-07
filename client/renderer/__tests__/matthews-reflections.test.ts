// @vitest-environment node
import { describe, expect, it } from "vitest";
import { matthewsReflectionFile } from "../lib/matthews-reflections";

// What splitMtz's import of a dropped MTZ leaves (#609)
const files = [
  { uuid: "aaaa-1111", type: "application/CCP4-mtz-freerflag", job_param_name: "FREEOUT" },
  { uuid: "bbbb-2222", type: "application/CCP4-mtz-observed", job_param_name: "OBSOUT" },
  { uuid: "cccc-3333", type: "text/plain", job_param_name: "LOG" },
];

describe("matthewsReflectionFile", () => {
  it("prefers the observed data, as a dashless database id", () => {
    expect(matthewsReflectionFile(files)).toBe("bbbb2222");
  });

  it("takes any MTZ when there are no observed data", () => {
    expect(matthewsReflectionFile([files[0], files[2]])).toBe("aaaa1111");
  });

  it("gives nothing when the import wrote no MTZ", () => {
    expect(matthewsReflectionFile([files[2]])).toBeNull();
  });
});
