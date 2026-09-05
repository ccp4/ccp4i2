/**
 * Content sniffing for ambiguous crystallographic file extensions.
 *
 * The cases that matter are the ones where the extension lies: a `.cif` may
 * be coordinates, structure factors or a restraint dictionary, and routing
 * the wrong one to a coordinate import produces a job that fails or quietly
 * does the wrong thing. These run against the real files in `demo_data/`,
 * because handwritten fixtures would not have caught the case that broke
 * the first draft — a refined coordinate file carrying `_chem_comp`.
 */
import * as fs from "fs";
import * as path from "path";
import { fileURLToPath } from "url";
import { describe, expect, it } from "vitest";
import { classifyCifText, classifyXmlText } from "@/lib/text-file-sniffer";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const REPO_ROOT = path.resolve(__dirname, "../../..");
const DEMO = path.join(REPO_ROOT, "server/ccp4i2/demo_data");

/** Read the same bounded prefix the sniffer itself reads. */
function readPrefix(relativePath: string): string {
  const full = path.join(REPO_ROOT, relativePath);
  const fd = fs.openSync(full, "r");
  try {
    const buffer = Buffer.alloc(64 * 1024);
    const read = fs.readSync(fd, buffer, 0, buffer.length, 0);
    return buffer.subarray(0, read).toString("utf8");
  } finally {
    fs.closeSync(fd);
  }
}

describe("classifyCifText", () => {
  it.each([
    // A refined structure. Note this file also contains _chem_comp, which an
    // earlier draft misread as a dictionary — coordinates must win.
    ["server/ccp4i2/demo_data/ahir/ahir_model.cif", "coordinates"],
    ["server/ccp4i2/demo_data/mdm2/4qo4.cif", "coordinates"],
    // Structure factors as downloaded from the PDB.
    ["server/ccp4i2/demo_data/mdm2/4hg7-sf.cif", "reflections"],
    // eLBOW/acedrg-style restraint dictionaries.
    ["server/ccp4i2/demo_data/baz2b/BAZ2BA-x839-LIG.cif", "dictionary"],
    ["server/ccp4i2/demo_data/baz2b/BAZ2BA-x828-LIG.cif", "dictionary"],
  ])("classifies %s as %s", (file, expected) => {
    expect(classifyCifText(readPrefix(file))).toBe(expected);
  });

  it("returns unknown rather than guessing at unrelated text", () => {
    expect(classifyCifText("hello world\nnot a cif at all\n")).toBe("unknown");
    expect(classifyCifText("")).toBe("unknown");
  });

  it("does not match a category named only mid-line", () => {
    // "_atom_site." inside prose must not make this a coordinate file.
    expect(classifyCifText("# see the _atom_site. category for details\n")).toBe(
      "unknown"
    );
  });
});

describe("classifyXmlText", () => {
  it("recognises a CCP4i2 ASU contents file", () => {
    expect(
      classifyXmlText(readPrefix("server/ccp4i2/demo_data/gamma/gamma.asu.xml"))
    ).toBe("asu_contents");
  });

  it.each([
    "server/ccp4i2/wrappers/areaimol/areaimol.def.xml",
    "server/ccp4i2/wrappers/qtpisa/script/qtpisa.scene.xml",
  ])("does not claim %s as ASU contents", (file) => {
    expect(classifyXmlText(readPrefix(file))).toBe("unknown");
  });

  it("returns unknown for non-XML text", () => {
    expect(classifyXmlText("ATOM      1  N   MET A   1\n")).toBe("unknown");
  });
});

describe("demo_data corpus", () => {
  // A guard against silent regression: every .cif we ship should classify as
  // something. A new "unknown" here means a writer we no longer recognise.
  it("classifies every demo_data .cif file", () => {
    const cifs = fs
      .readdirSync(DEMO, { recursive: true, encoding: "utf8" })
      .filter((f) => typeof f === "string" && f.endsWith(".cif"));

    expect(cifs.length).toBeGreaterThan(0);
    const unknown = cifs.filter(
      (f) =>
        classifyCifText(readPrefix(path.join("server/ccp4i2/demo_data", f))) ===
        "unknown"
    );
    expect(unknown).toEqual([]);
  });
});
