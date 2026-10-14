import { describe, expect, it } from "vitest";
import { DirectoryEntry, directoryRows } from "../lib/directory-rows";

const dir = (name: string): DirectoryEntry => ({ path: `/p/${name}`, name, type: "directory" });
const file = (name: string): DirectoryEntry => ({ path: `/p/${name}`, name, type: "file", size: 1 });

const labels = (rows: ReturnType<typeof directoryRows>) =>
  rows.map((r) => (r.kind === "loading" ? `${"  ".repeat(r.depth)}…` : `${"  ".repeat(r.depth)}${r.entry.name}`));

describe("directoryRows", () => {
  it("draws only the root until a folder is opened", () => {
    const rows = directoryRows({ "": [file("b.log"), dir("job_10"), dir("job_2")] }, new Set());
    expect(labels(rows)).toEqual(["job_2", "job_10", "b.log"]);
  });

  it("puts an open folder's entries under it, and a loading row while it is unlisted", () => {
    const listings = {
      "": [dir("pandda2_out"), file("program.xml")],
      pandda2_out: [dir("processed_datasets"), file("input.yaml")],
    };
    const open = new Set(["pandda2_out", "pandda2_out/processed_datasets"]);
    expect(labels(directoryRows(listings, open))).toEqual([
      "pandda2_out",
      "  processed_datasets",
      "    …",
      "  input.yaml",
      "program.xml",
    ]);
  });

  it("filters files but keeps folders, whose contents are not known yet", () => {
    const rows = directoryRows(
      { "": [dir("sub"), file("x.log"), file("x.mtz")] },
      new Set(),
      (e) => e.name.endsWith(".log")
    );
    expect(labels(rows)).toEqual(["sub", "x.log"]);
  });

  it("gives each row its path relative to the root", () => {
    const rows = directoryRows({ "": [dir("a")], a: [file("f")] }, new Set(["a"]));
    expect(rows.map((r) => r.relative)).toEqual(["a", "a/f"]);
  });
});
