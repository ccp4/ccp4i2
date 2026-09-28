/**
 * CDataFile renders once per file in a task. Subscribing to the digest or the
 * file content from inside it fetched both for every file on mount: 948
 * requests for a 158-dataset PanDDA job, enough to exhaust the database's
 * connection slots on DDU and starve the run check. The element invalidates
 * those caches after a change instead; the elements that display a digest or
 * content subscribe for themselves.
 */
import { readFileSync } from "fs";
import { join } from "path";
import { describe, expect, it } from "vitest";

const source = readFileSync(
  join(__dirname, "..", "components", "task", "task-elements", "cdatafile.tsx"),
  "utf8"
);

describe("CDataFile", () => {
  it("does not subscribe to a digest or the file content for every file it renders", () => {
    expect(source).not.toMatch(/useFileDigest\s*\(/);
    expect(source).not.toMatch(/useFileContent\s*\(/);
  });

  it("still invalidates both caches after a file change", () => {
    expect(source).toMatch(/mutateSwr\(/);
    expect(source).toMatch(/digest\?object_path=/);
    expect(source).toMatch(/files_by_uuid\//);
    expect(source).toMatch(/mutateContent\(\), mutateDigest\(\)/);
  });
});
