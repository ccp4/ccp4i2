/**
 * Moorhen Scene parser / serialiser tests.
 *
 * The point of these tests is to lock down the format itself —
 * what's valid, what isn't, and that a scene survives round-trip
 * (parse → serialise → re-parse → deep-equal). The lifter (Moorhen state
 * → scene) and the resolver (scene → Moorhen session) are tested
 * separately and are not exercised here.
 */
import { describe, it, expect } from "vitest";

import {
  parseScene,
  serialiseScene,
  validateScene,
  SceneParseError,
} from "../lib/moorhen-scene";
import {
  MoorhenScene,
  SCENE_SCHEMA_VERSION,
} from "../types/moorhen-scene";

// --------------------------------------------------------------------------
// Worked example: matches the design conversation.
// --------------------------------------------------------------------------

const KINASE_YAML = `scene: kinase-domain-view
version: 1
authoredIn:
  projectId: 3f8a-aaaa-bbbb-cccc-uuid
  projectName: kinase-2026
  createdAt: 2026-05-22T14:30:00Z
  createdBy: martin.noble@ncl.ac.uk
files:
  - name: protein
    projectId: 3f8a-aaaa-bbbb-cccc-uuid
    projectName: kinase-2026
    job: 14
    param: XYZOUT
  - name: reference
    url: https://ddudatabase/api/ccp4i2/fileBy/abc
  - name: external
    path: /abs/path/to/foo.pdb
domains:
  - name: N-lobe
    chain: A
    range: 1-120
    color: "#4b8bbe"
  - name: hinge
    chain: A
    range: 121-130
    color: "#f1c40f"
  - name: C-lobe
    chain: A
    range: 131-300
    color: "#e74c3c"
elements:
  - file: protein
    representations:
      - style: CRs
        selection: //A
        colour: by-domain
      - style: ligands
        selection: //*/LIG
        colour: "#2ecc71"
view:
  origin: [10.5, 20.0, -5.25]
  quat: [0, 0, 0, -1]
  zoom: 1.5
resolver:
  onMissingResidues: clamp-and-log
`;

describe("Moorhen scene — worked kinase example", () => {
  it("parses the example without errors", () => {
    const scene = parseScene(KINASE_YAML);
    expect(scene.scene).toBe("kinase-domain-view");
    expect(scene.version).toBe(SCENE_SCHEMA_VERSION);
    expect(scene.authoredIn?.projectName).toBe("kinase-2026");
    expect(scene.files).toHaveLength(3);
    expect(scene.files![0]).toMatchObject({
      name: "protein",
      job: 14,
      param: "XYZOUT",
      projectId: "3f8a-aaaa-bbbb-cccc-uuid",
    });
    expect(scene.domains).toHaveLength(3);
    expect(scene.domains![1]).toEqual({
      name: "hinge",
      chain: "A",
      range: "121-130",
      color: "#f1c40f",
    });
    expect(scene.elements).toHaveLength(1);
    expect(scene.elements![0].representations).toHaveLength(2);
    expect(scene.elements![0].representations![0].colour).toBe("by-domain");
    expect(scene.elements![0].representations![1].colour).toBe("#2ecc71");
    expect(scene.view?.origin).toEqual([10.5, 20.0, -5.25]);
    expect(scene.resolver?.onMissingResidues).toBe("clamp-and-log");
  });

  it("round-trips: parse → serialise → re-parse is stable", () => {
    const first = parseScene(KINASE_YAML);
    const yaml = serialiseScene(first);
    const second = parseScene(yaml);
    expect(second).toEqual(first);
  });
});

// --------------------------------------------------------------------------
// Validation: each error class.
// --------------------------------------------------------------------------

describe("Moorhen scene — validation errors", () => {
  it("rejects non-mapping root", () => {
    expect(() => parseScene("- not a map")).toThrow(SceneParseError);
  });

  it("rejects wrong schema version", () => {
    const yaml = `scene: x\nversion: 99\n`;
    expect(() => parseScene(yaml)).toThrow(/unsupported scene version 99/);
  });

  it("rejects a file ref with no resolution method", () => {
    const yaml = `scene: x\nversion: 1\nfiles:\n  - name: orphan\n`;
    expect(() => parseScene(yaml)).toThrow(
      /must set one of: pdb, url, path, bundle, fileId/,
    );
  });

  it("accepts a pdb: file ref with a 4-char id", () => {
    const yaml = `scene: x\nversion: 1\nfiles:\n  - name: structure\n    pdb: 1m17\n`;
    const scene = parseScene(yaml);
    expect(scene.files![0].pdb).toBe("1m17");
  });

  it("accepts a pdb: file ref with an extended id", () => {
    const yaml = `scene: x\nversion: 1\nfiles:\n  - name: s\n    pdb: pdb_00001m17\n`;
    expect(() => parseScene(yaml)).not.toThrow();
  });

  it("rejects an obviously malformed pdb id", () => {
    const yaml = `scene: x\nversion: 1\nfiles:\n  - name: s\n    pdb: my-favourite-protein\n`;
    expect(() => parseScene(yaml)).toThrow(/does not look like a PDB ID/);
  });

  it("requires projectId when fileId is set", () => {
    const yaml = `scene: x\nversion: 1\nfiles:\n  - name: f\n    fileId: 42\n`;
    expect(() => parseScene(yaml)).toThrow(
      /projectId.*required when fileId or job\+param is set/,
    );
  });

  it("requires job and param to be set together", () => {
    const yaml = `scene: x\nversion: 1\nfiles:\n  - name: f\n    projectId: p\n    job: 14\n`;
    expect(() => parseScene(yaml)).toThrow(/job and param must be set together/);
  });

  it("rejects duplicate file names", () => {
    const yaml = `scene: x\nversion: 1\nfiles:\n  - name: a\n    url: u1\n  - name: a\n    url: u2\n`;
    expect(() => parseScene(yaml)).toThrow(/duplicate file name "a"/);
  });

  it("rejects malformed residue range", () => {
    const yaml = `scene: x\nversion: 1\ndomains:\n  - name: d\n    chain: A\n    range: not-a-range\n    color: "#ffffff"\n`;
    expect(() => parseScene(yaml)).toThrow(/must be "start-end"/);
  });

  it("rejects reversed residue range", () => {
    const yaml = `scene: x\nversion: 1\ndomains:\n  - name: d\n    chain: A\n    range: 100-50\n    color: "#ffffff"\n`;
    expect(() => parseScene(yaml)).toThrow(/end \(50\) must be >= start \(100\)/);
  });

  it("rejects non-hex domain colour", () => {
    const yaml = `scene: x\nversion: 1\ndomains:\n  - name: d\n    chain: A\n    range: 1-10\n    color: red\n`;
    expect(() => parseScene(yaml)).toThrow(/must be hex like "#rrggbb"/);
  });

  it("rejects an element referencing an unknown file", () => {
    const yaml = `scene: x\nversion: 1\nfiles:\n  - name: a\n    url: u\nelements:\n  - file: b\n`;
    expect(() => parseScene(yaml)).toThrow(/unknown file "b"/);
  });

  it("rejects unknown named colour", () => {
    const yaml = `scene: x\nversion: 1\nfiles:\n  - name: a\n    url: u\nelements:\n  - file: a\n    representations:\n      - style: CRs\n        colour: by-something-weird\n`;
    expect(() => parseScene(yaml)).toThrow(/unknown colour "by-something-weird"/);
  });

  it("accepts the raw colour escape hatch", () => {
    const yaml = `scene: x\nversion: 1\nfiles:\n  - name: a\n    url: u\nelements:\n  - file: a\n    representations:\n      - style: CRs\n        colour:\n          raw:\n            ruleType: custom\n            args: ["//A/1-50^#ff0000", 42]\n            isMultiColourRule: true\n`;
    const scene = parseScene(yaml);
    const c = scene.elements![0].representations![0].colour;
    expect(c).toEqual({
      raw: {
        ruleType: "custom",
        args: ["//A/1-50^#ff0000", 42],
        isMultiColourRule: true,
        applyColourToNonCarbonAtoms: undefined,
      },
    });
  });

  it("rejects bad resolver policy", () => {
    const yaml = `scene: x\nversion: 1\nresolver:\n  onMissingResidues: pretend-they-exist\n`;
    expect(() => parseScene(yaml)).toThrow(
      /must be "clamp-and-log" or "strict"/,
    );
  });

  it("collects multiple errors in a single pass", () => {
    const yaml = `scene: ""\nversion: 1\ndomains:\n  - name: d\n    chain: A\n    range: oops\n    color: red\n`;
    const { errors } = validateScene({
      scene: "x",
      version: 1,
      domains: [{ name: "d", chain: "A", range: "oops", color: "red" }],
    });
    expect(errors.length).toBeGreaterThanOrEqual(2);
    expect(errors.map((e) => e.path)).toEqual(
      expect.arrayContaining(["domains[0].range", "domains[0].color"]),
    );
  });
});

// --------------------------------------------------------------------------
// Serialisation: omits empty/undefined fields for stable round-trip.
// --------------------------------------------------------------------------

describe("Moorhen scene — serialiser canonicalisation", () => {
  it("omits empty top-level blocks", () => {
    const scene: MoorhenScene = { scene: "minimal", version: 1 };
    const yaml = serialiseScene(scene);
    expect(yaml).toContain("scene: minimal");
    expect(yaml).toContain("version: 1");
    expect(yaml).not.toContain("files:");
    expect(yaml).not.toContain("domains:");
    expect(yaml).not.toContain("elements:");
    expect(yaml).not.toContain("view:");
    expect(yaml).not.toContain("resolver:");
    expect(yaml).not.toContain("authoredIn:");
  });

  it("omits undefined fields inside authoredIn", () => {
    const scene: MoorhenScene = {
      scene: "x",
      version: 1,
      authoredIn: { projectName: "only-this" },
    };
    const yaml = serialiseScene(scene);
    expect(yaml).toContain("projectName: only-this");
    expect(yaml).not.toContain("projectId:");
    expect(yaml).not.toContain("createdAt:");
  });

  it("preserves field order: scene, version, authoredIn, files, domains, elements, view, resolver", () => {
    const scene = parseScene(KINASE_YAML);
    const yaml = serialiseScene(scene);
    const order = [
      "scene:",
      "version:",
      "authoredIn:",
      "files:",
      "domains:",
      "elements:",
      "view:",
      "resolver:",
    ].map((k) => yaml.indexOf(k));
    for (let i = 1; i < order.length; i++) {
      expect(order[i]).toBeGreaterThan(order[i - 1]);
    }
  });
});

// --------------------------------------------------------------------------
// chain: extended forms (string | "*" | array)
// --------------------------------------------------------------------------

describe("Moorhen scene — chain field forms", () => {
  it("accepts a single-string chain (the legacy form)", () => {
    const yaml = `scene: x\nversion: 1\ndomains:\n  - { name: d, chain: A, range: 1-10, color: "#ffffff" }\n`;
    const scene = parseScene(yaml);
    expect(scene.domains![0].chain).toBe("A");
  });

  it('accepts "*" as a wildcard chain string', () => {
    const yaml = `scene: x\nversion: 1\ndomains:\n  - { name: d, chain: "*", range: 1-10, color: "#ffffff" }\n`;
    const scene = parseScene(yaml);
    expect(scene.domains![0].chain).toBe("*");
  });

  it("accepts a chain list", () => {
    const yaml = `scene: x\nversion: 1\ndomains:\n  - { name: d, chain: [A, B, C], range: 1-10, color: "#ffffff" }\n`;
    const scene = parseScene(yaml);
    expect(scene.domains![0].chain).toEqual(["A", "B", "C"]);
  });

  it("rejects an empty chain list", () => {
    const yaml = `scene: x\nversion: 1\ndomains:\n  - { name: d, chain: [], range: 1-10, color: "#ffffff" }\n`;
    expect(() => parseScene(yaml)).toThrow(/chain list must not be empty/);
  });

  it("rejects a chain list with non-string entries", () => {
    const yaml = `scene: x\nversion: 1\ndomains:\n  - { name: d, chain: [A, 42], range: 1-10, color: "#ffffff" }\n`;
    expect(() => parseScene(yaml)).toThrow(/chain list entries must all be strings/);
  });

  it("rejects a chain field of the wrong type", () => {
    const yaml = `scene: x\nversion: 1\ndomains:\n  - { name: d, chain: 42, range: 1-10, color: "#ffffff" }\n`;
    expect(() => parseScene(yaml)).toThrow(/must be a string or sequence of strings/);
  });

  it("round-trips a wildcard and a list through serialise → parse", () => {
    const yaml = `scene: x\nversion: 1\ndomains:\n  - { name: a, chain: "*",     range: 1-10,  color: "#ffffff" }\n  - { name: b, chain: [A, B],   range: 11-20, color: "#000000" }\n`;
    const first = parseScene(yaml);
    const second = parseScene(serialiseScene(first));
    expect(second).toEqual(first);
  });
});

// --------------------------------------------------------------------------
// superpose: block
// --------------------------------------------------------------------------

describe("Moorhen scene — superpose block", () => {
  it("accepts an SSM entry referencing two files", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: apo,  pdb: 1m17 }
  - { name: holo, pdb: 1m14 }
superpose:
  - { method: ssm, move: holo, onto: apo, movChain: A, refChain: A }
`;
    const scene = parseScene(yaml);
    expect(scene.superpose).toHaveLength(1);
    expect(scene.superpose![0]).toEqual({
      method: "ssm",
      move: "holo",
      onto: "apo",
      movChain: "A",
      refChain: "A",
    });
  });

  it("accepts an LSQ entry with explicit residue matches", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: ref, pdb: 1m17 }
  - { name: mov, pdb: 1m14 }
superpose:
  - method: lsq
    move: mov
    onto: ref
    matches:
      - { refChain: A, refRange: 1-100, movChain: A, movRange: 1-100 }
      - { refChain: B, refRange: 50-150, movChain: B, movRange: 60-160 }
    matchType: ca
`;
    const scene = parseScene(yaml);
    expect(scene.superpose![0]).toEqual({
      method: "lsq",
      move: "mov",
      onto: "ref",
      matchType: "ca",
      matches: [
        { refChain: "A", refRange: "1-100", movChain: "A", movRange: "1-100" },
        { refChain: "B", refRange: "50-150", movChain: "B", movRange: "60-160" },
      ],
    });
  });

  it("rejects an unknown method", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: a, pdb: 1m17 }
  - { name: b, pdb: 1m14 }
superpose:
  - { method: tm-align, move: b, onto: a }
`;
    expect(() => parseScene(yaml)).toThrow(/must be "ssm" or "lsq"/);
  });

  it("rejects move == onto", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: a, pdb: 1m17 }
superpose:
  - { method: ssm, move: a, onto: a, movChain: A, refChain: A }
`;
    expect(() => parseScene(yaml)).toThrow(/cannot superpose a file onto itself/);
  });

  it("rejects an SSM entry missing movChain", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: a, pdb: 1m17 }
  - { name: b, pdb: 1m14 }
superpose:
  - { method: ssm, move: b, onto: a, refChain: A }
`;
    expect(() => parseScene(yaml)).toThrow(/movChain.*required for ssm/);
  });

  it("rejects an LSQ entry missing matches", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: a, pdb: 1m17 }
  - { name: b, pdb: 1m14 }
superpose:
  - { method: lsq, move: b, onto: a }
`;
    expect(() => parseScene(yaml)).toThrow(/matches.*required for lsq/);
  });

  it("rejects an LSQ entry with empty matches", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: a, pdb: 1m17 }
  - { name: b, pdb: 1m14 }
superpose:
  - { method: lsq, move: b, onto: a, matches: [] }
`;
    expect(() => parseScene(yaml)).toThrow(/at least one match entry/);
  });

  it("rejects a malformed residue range in an LSQ match", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: a, pdb: 1m17 }
  - { name: b, pdb: 1m14 }
superpose:
  - method: lsq
    move: b
    onto: a
    matches:
      - { refChain: A, refRange: oops, movChain: A, movRange: 1-100 }
`;
    expect(() => parseScene(yaml)).toThrow(/must be "start-end"/);
  });

  it("rejects superpose entries referencing unknown file names", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: a, pdb: 1m17 }
superpose:
  - { method: ssm, move: ghost, onto: a, movChain: A, refChain: A }
`;
    expect(() => parseScene(yaml)).toThrow(/unknown file "ghost"/);
  });

  it("round-trips an SSM entry through serialise → parse", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: apo,  pdb: 1m17 }
  - { name: holo, pdb: 1m14 }
superpose:
  - { method: ssm, move: holo, onto: apo, movChain: A, refChain: A }
`;
    const first = parseScene(yaml);
    const second = parseScene(serialiseScene(first));
    expect(second.superpose).toEqual(first.superpose);
  });

  it("accepts the LSQ chain+range shorthand", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: a, pdb: 1m17 }
  - { name: b, pdb: 1m14 }
superpose:
  - { method: lsq, move: b, onto: a, chain: A, range: 104-260 }
`;
    const scene = parseScene(yaml);
    expect(scene.superpose![0]).toEqual({
      method: "lsq",
      move: "b",
      onto: "a",
      chain: "A",
      range: "104-260",
    });
  });

  it("rejects mixing `matches` with the chain+range shorthand", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: a, pdb: 1m17 }
  - { name: b, pdb: 1m14 }
superpose:
  - method: lsq
    move: b
    onto: a
    chain: A
    range: 1-100
    matches:
      - { refChain: A, refRange: 1-100, movChain: A, movRange: 1-100 }
`;
    expect(() => parseScene(yaml)).toThrow(
      /use either `matches` or the `chain`\+`range` shorthand/,
    );
  });

  it("rejects a shorthand entry missing range", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: a, pdb: 1m17 }
  - { name: b, pdb: 1m14 }
superpose:
  - { method: lsq, move: b, onto: a, chain: A }
`;
    expect(() => parseScene(yaml)).toThrow(/range.*required when chain is set/);
  });

  it("rejects a shorthand entry with a malformed range", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: a, pdb: 1m17 }
  - { name: b, pdb: 1m14 }
superpose:
  - { method: lsq, move: b, onto: a, chain: A, range: nope }
`;
    expect(() => parseScene(yaml)).toThrow(/must be "start-end"/);
  });

  it("round-trips a shorthand LSQ entry through serialise → parse", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: a, pdb: 1m17 }
  - { name: b, pdb: 1m14 }
superpose:
  - { method: lsq, move: b, onto: a, chain: A, range: 104-260, matchType: main }
`;
    const first = parseScene(yaml);
    const second = parseScene(serialiseScene(first));
    expect(second.superpose).toEqual(first.superpose);
  });
});

// --------------------------------------------------------------------------
// Dictionary scoping (kind: dictionary, globalDictionaries, element.dictionaries)
// --------------------------------------------------------------------------

describe("Moorhen scene — dictionary scoping", () => {
  it("accepts a dictionary file ref with kind: dictionary", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: lig-A, kind: dictionary, url: https://example/dict-A.cif }
`;
    const scene = parseScene(yaml);
    expect(scene.files![0].kind).toBe("dictionary");
    expect(scene.files![0].url).toBe("https://example/dict-A.cif");
  });

  it("rejects an unknown kind value", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: f, kind: weird, url: https://example/f.cif }
`;
    expect(() => parseScene(yaml)).toThrow(/must be "coordinates", "dictionary", or "mtz"/);
  });

  it("rejects pdb: on a dictionary ref (pdb makes no sense for dicts)", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: d, kind: dictionary, pdb: 1abc }
`;
    expect(() => parseScene(yaml)).toThrow(/pdb: is only valid for coordinate refs/);
  });

  it("accepts globalDictionaries referencing a dictionary file", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: coords, pdb: 1abc }
  - { name: cofactor, kind: dictionary, url: https://example/cof.cif }
globalDictionaries: [cofactor]
`;
    const scene = parseScene(yaml);
    expect(scene.globalDictionaries).toEqual(["cofactor"]);
  });

  it("rejects globalDictionaries referencing a non-dictionary file", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: coords, pdb: 1abc }
  - { name: other-coords, pdb: 1xyz }
globalDictionaries: [other-coords]
`;
    expect(() => parseScene(yaml)).toThrow(
      /file "other-coords" is not a dictionary/,
    );
  });

  it("rejects globalDictionaries referencing an unknown file", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: coords, pdb: 1abc }
globalDictionaries: [missing]
`;
    expect(() => parseScene(yaml)).toThrow(/unknown file "missing"/);
  });

  it("accepts an element with dictionaries pointing at dict files", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: coords-A, pdb: 1abc }
  - { name: lig-A, kind: dictionary, url: https://example/A.cif }
elements:
  - file: coords-A
    dictionaries: [lig-A]
    representations:
      - { style: ligands, selection: "//*/LIG", colour: "#2ecc71" }
`;
    const scene = parseScene(yaml);
    expect(scene.elements![0].dictionaries).toEqual(["lig-A"]);
  });

  it("rejects an element.dictionaries referencing a non-dictionary file", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: coords-A, pdb: 1abc }
  - { name: coords-B, pdb: 1xyz }
elements:
  - file: coords-A
    dictionaries: [coords-B]
`;
    expect(() => parseScene(yaml)).toThrow(
      /file "coords-B" is not a dictionary/,
    );
  });

  it("round-trips a scene with dicts through serialise → parse", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: coords-A, pdb: 1abc }
  - { name: lig-A, kind: dictionary, url: https://example/A.cif }
  - { name: cofactor, kind: dictionary, url: https://example/cof.cif }
globalDictionaries: [cofactor]
elements:
  - file: coords-A
    dictionaries: [lig-A]
    representations:
      - { style: ligands, selection: "//*/LIG", colour: "#2ecc71" }
`;
    const first = parseScene(yaml);
    const second = parseScene(serialiseScene(first));
    expect(second).toEqual(first);
  });
});

describe("Moorhen scene — inline dict (cifText)", () => {
  it("accepts a dictionary with cifText as its sole source", () => {
    const yaml = `scene: x
version: 1
files:
  - name: lig
    kind: dictionary
    cifText: |
      data_comp_LIG
      _chem_comp.id LIG
`;
    const scene = parseScene(yaml);
    expect(scene.files![0].kind).toBe("dictionary");
    expect(scene.files![0].cifText).toContain("data_comp_LIG");
  });

  it("rejects cifText on a coordinate ref", () => {
    const yaml = `scene: x
version: 1
files:
  - name: coords
    cifText: some-cif-content
`;
    expect(() => parseScene(yaml)).toThrow(
      /cifText: is only valid for dictionary refs/,
    );
  });

  it("includes cifText in the must-set-one-of message", () => {
    const yaml = `scene: x
version: 1
files:
  - name: orphan
    kind: dictionary
`;
    expect(() => parseScene(yaml)).toThrow(
      /must set one of:.*cifText/,
    );
  });

  it("round-trips an inline-dict scene", () => {
    const yaml = `scene: x
version: 1
files:
  - name: coords-A
    pdb: 1abc
  - name: lig
    kind: dictionary
    cifText: |
      data_comp_LIG
      _chem_comp.id LIG
elements:
  - file: coords-A
    dictionaries: [lig]
`;
    const first = parseScene(yaml);
    const second = parseScene(serialiseScene(first));
    expect(second).toEqual(first);
  });
});

describe("Moorhen scene — bundle refs", () => {
  it("accepts a coord file ref whose only source is bundle:", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: complex, bundle: assets/complex.pdb }
`;
    const scene = parseScene(yaml);
    expect(scene.files![0].bundle).toBe("assets/complex.pdb");
  });

  it("accepts a dictionary ref whose only source is bundle:", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: lig, kind: dictionary, bundle: assets/lig.cif }
`;
    const scene = parseScene(yaml);
    expect(scene.files![0].kind).toBe("dictionary");
    expect(scene.files![0].bundle).toBe("assets/lig.cif");
  });

  it("round-trips a full bundle scene", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: complex, bundle: assets/complex.pdb }
  - { name: lig,     kind: dictionary, bundle: assets/lig.cif }
elements:
  - file: complex
    dictionaries: [lig]
    representations:
      - { style: ligands, selection: "//*/LIG", colour: "#2ecc71" }
`;
    const first = parseScene(yaml);
    const second = parseScene(serialiseScene(first));
    expect(second).toEqual(first);
  });
});

describe("Moorhen scene — bare-int residue ranges", () => {
  it("accepts a bare integer in a domain range and normalises to N-N", () => {
    const yaml = `scene: x
version: 1
domains:
  - { name: catalytic, chain: A, range: 115, color: "#e74c3c" }
`;
    const scene = parseScene(yaml);
    expect(scene.domains![0].range).toBe("115-115");
  });

  it("still accepts the string form for domain range", () => {
    const yaml = `scene: x
version: 1
domains:
  - { name: helix, chain: A, range: 100-130, color: "#e74c3c" }
`;
    const scene = parseScene(yaml);
    expect(scene.domains![0].range).toBe("100-130");
  });

  it("rejects a non-integer number as range", () => {
    // Floats and other numeric oddities shouldn't sneak through —
    // residues are integers.
    const yaml = `scene: x
version: 1
domains:
  - { name: d, chain: A, range: 1.5, color: "#e74c3c" }
`;
    expect(() => parseScene(yaml)).toThrow(/must be a string or integer/);
  });

  it("accepts bare ints in lsq match refRange/movRange", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: a, pdb: 1m17 }
  - { name: b, pdb: 1m14 }
superpose:
  - method: lsq
    move: b
    onto: a
    matches:
      - { refChain: A, refRange: 100, movChain: A, movRange: 100 }
`;
    const scene = parseScene(yaml);
    const sp = scene.superpose![0] as { matches: { refRange: string; movRange: string }[] };
    expect(sp.matches[0].refRange).toBe("100-100");
    expect(sp.matches[0].movRange).toBe("100-100");
  });

  it("accepts a bare int in the lsq chain+range shorthand", () => {
    const yaml = `scene: x
version: 1
files:
  - { name: a, pdb: 1m17 }
  - { name: b, pdb: 1m14 }
superpose:
  - { method: lsq, move: b, onto: a, chain: A, range: 200 }
`;
    const scene = parseScene(yaml);
    const sp = scene.superpose![0] as { range: string };
    expect(sp.range).toBe("200-200");
  });
});

describe("parseScene — maps block", () => {
  it("parses a minimal MTZ map", () => {
    const yaml = `scene: x
version: 1
files:
  - name: mtz1
    kind: mtz
    url: https://example/best.mtz
maps:
  - name: best
    file: mtz1
    columns: { F: FWT, PHI: PHWT }
`;
    const scene = parseScene(yaml);
    expect(scene.maps).toHaveLength(1);
    expect(scene.maps![0]).toMatchObject({
      name: "best",
      file: "mtz1",
      columns: { F: "FWT", PHI: "PHWT" },
    });
  });

  it("captures full render state and difference-map colours", () => {
    const yaml = `scene: x
version: 1
files:
  - name: mtz1
    kind: mtz
    url: https://example/diff.mtz
maps:
  - name: fofc
    file: mtz1
    columns: { F: DELFWT, PHI: PHDELWT }
    isDifference: true
    contourLevel: 3.2
    radius: 17.5
    alpha: 0.85
    style: solid
    positiveColour: "#00ff00"
    negativeColour: "#ff0000aa"
    visible: false
activeMap: fofc
`;
    const scene = parseScene(yaml);
    expect(scene.maps![0]).toMatchObject({
      isDifference: true,
      contourLevel: 3.2,
      radius: 17.5,
      alpha: 0.85,
      style: "solid",
      positiveColour: "#00ff00",
      negativeColour: "#ff0000aa",
      visible: false,
    });
    expect(scene.activeMap).toBe("fofc");
  });

  it("rejects a map referencing a non-existent file", () => {
    const yaml = `scene: x
version: 1
files:
  - name: mtz1
    kind: mtz
    url: https://example/best.mtz
maps:
  - name: best
    file: nope
    columns: { F: FWT, PHI: PHWT }
`;
    expect(() => parseScene(yaml)).toThrow(/unknown file "nope"/);
  });

  it("rejects a map pointing at a coord file", () => {
    const yaml = `scene: x
version: 1
files:
  - name: coords
    url: https://example/model.pdb
maps:
  - name: best
    file: coords
    columns: { F: FWT, PHI: PHWT }
`;
    expect(() => parseScene(yaml)).toThrow(/must be a file with kind: "mtz"/);
  });

  it("rejects a map without F+PHI (and no calcStructFact)", () => {
    const yaml = `scene: x
version: 1
files:
  - name: mtz1
    kind: mtz
    url: https://example/best.mtz
maps:
  - name: best
    file: mtz1
    columns: { F: FWT }
`;
    expect(() => parseScene(yaml)).toThrow(/must set F \+ PHI/);
  });

  it("accepts calcStructFact when Fobs+SigFobs are supplied", () => {
    const yaml = `scene: x
version: 1
files:
  - name: mtz1
    kind: mtz
    url: https://example/data.mtz
maps:
  - name: calc
    file: mtz1
    columns: { Fobs: FP, SigFobs: SIGFP, FreeR: FREE, calcStructFact: true }
`;
    const scene = parseScene(yaml);
    expect(scene.maps![0].columns!.calcStructFact).toBe(true);
  });

  it("parses a real-space map / mask (kind: map, no columns)", () => {
    const yaml = `scene: x
version: 1
files:
  - name: msk
    kind: map
    url: https://example/mask.map
maps:
  - name: domainA
    file: msk
    isMask: true
`;
    const scene = parseScene(yaml);
    expect(scene.maps).toHaveLength(1);
    expect(scene.maps![0]).toMatchObject({ name: "domainA", file: "msk", isMask: true });
    expect(scene.maps![0].columns).toBeUndefined();
  });

  it("rejects columns on a kind: map (real-space) file", () => {
    const yaml = `scene: x
version: 1
files:
  - name: msk
    kind: map
    url: https://example/mask.map
maps:
  - name: domainA
    file: msk
    columns: { F: FWT, PHI: PHWT }
`;
    expect(() => parseScene(yaml)).toThrow(/not allowed for kind: "map"/);
  });

  it("requires kind mtz or map for a maps[] file", () => {
    const yaml = `scene: x
version: 1
files:
  - name: coords
    pdb: 1abc
maps:
  - name: m
    file: coords
`;
    expect(() => parseScene(yaml)).toThrow(/must be a file with kind: "mtz" or "map"/);
  });

  it("rejects an unknown style", () => {
    const yaml = `scene: x
version: 1
files:
  - name: mtz1
    kind: mtz
    url: https://example/best.mtz
maps:
  - name: best
    file: mtz1
    columns: { F: FWT, PHI: PHWT }
    style: chunky
`;
    expect(() => parseScene(yaml)).toThrow(/must be "lines", "solid", or "lit-lines"/);
  });

  it("rejects activeMap referencing a non-existent entry", () => {
    const yaml = `scene: x
version: 1
files:
  - name: mtz1
    kind: mtz
    url: https://example/best.mtz
maps:
  - name: best
    file: mtz1
    columns: { F: FWT, PHI: PHWT }
activeMap: missing
`;
    expect(() => parseScene(yaml)).toThrow(/does not name any entry in maps:/);
  });

  it("round-trips maps through serialise → parse", async () => {
    const { serialiseScene } = await import("../lib/moorhen-scene");
    const scene = parseScene(`scene: round
version: 1
files:
  - name: mtz1
    kind: mtz
    url: https://example/best.mtz
maps:
  - name: best
    file: mtz1
    columns: { F: FWT, PHI: PHWT }
    contourLevel: 1.2
    colour: "#336699"
activeMap: best
`);
    const yaml = serialiseScene(scene);
    const parsed = parseScene(yaml);
    expect(parsed.maps![0]).toMatchObject({
      name: "best",
      file: "mtz1",
      columns: { F: "FWT", PHI: "PHWT" },
      contourLevel: 1.2,
      colour: "#336699",
    });
    expect(parsed.activeMap).toBe("best");
  });
});
