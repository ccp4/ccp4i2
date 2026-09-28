/**
 * The file elements render once per file in a task. Subscribing to the digest
 * or the file content from inside them fetched both for every file on mount:
 * 948 requests for a 158-dataset PanDDA job, enough to exhaust the database's
 * connection slots on DDU and starve the run check. They invalidate through
 * useJob's helpers after a change instead; only the code that displays a
 * digest subscribes, and the PDB element only when its builder can appear.
 */
import { readFileSync } from "fs";
import { join } from "path";
import { describe, expect, it } from "vitest";

const elements = join(__dirname, "..", "components", "task", "task-elements");
const read = (name: string) => readFileSync(join(elements, name), "utf8");

describe("file elements rendered once per file", () => {
  it.each(["cdatafile.tsx", "csimpledatafile.tsx", "cminimtzdatafile.tsx"])(
    "%s does not subscribe to a digest or the file content",
    (name) => {
      const source = read(name);
      expect(source).not.toMatch(/useFileDigest\s*\(/);
      expect(source).not.toMatch(/useFileContent\s*\(/);
      expect(source).toMatch(/mutateFileDigest\(/);
    }
  );

  it("cdatafile still invalidates both caches after a file change", () => {
    const source = read("cdatafile.tsx");
    expect(source).toMatch(/mutateContent = mutateFileContent/);
    expect(source).toMatch(/mutateContent\(\), mutateDigest\(\)/);
  });

  it("cpdbdatafile fetches a digest only when the atom-selection builder can appear", () => {
    const source = read("cpdbdatafile.tsx");
    expect(source).toMatch(/const digestPath =\s*hasFile && wantsComposition/);
    expect(source).toMatch(/wantsComposition = Boolean\(\s*qualifiers\?\.ifAtomSelection/);
  });
});

describe("useJob invalidation helpers", () => {
  const utils = readFileSync(join(__dirname, "..", "utils.ts"), "utf8");
  it("match the digest key with or without its cacheKey suffix, and every content key", () => {
    expect(utils).toMatch(/key === digestKey \|\| key\.startsWith\(`\$\{digestKey\}&`\)/);
    expect(utils).toMatch(/key\.startsWith\("files_by_uuid\/"\) &&\s*key\.endsWith\("\/download\/"\)/);
  });
});
