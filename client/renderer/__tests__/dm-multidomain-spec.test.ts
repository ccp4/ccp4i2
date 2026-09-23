/**
 * The dm_multidomain parameter grammar.
 *
 * ASSEMBLY and DOMAINS are stored as strings that the Python wrapper parses
 * (dm_ncs_lib.parse_assembly_rows / parse_segments), and i2run writes by hand.
 * The interface works in instances/roles/segments and serialises through
 * dm-spec, so these tests are what keeps the two ends speaking the same
 * language — particularly the terse form (a homomer must stay bare chain ids
 * and bare ranges) and the role rename, which has to reach both parameters or
 * the cross-reference between them breaks.
 */
import { describe, it, expect } from "vitest";

import {
  IMPLICIT_ROLE,
  Model,
  SegmentParseError,
  addRole,
  formatAssemblyRows,
  formatBodies,
  isTerse,
  overlappingBodies,
  parseAssemblyRows,
  parseBodies,
  parseSegments,
  referenceBounds,
  removeRole,
  renameRole,
  roleLabel,
  setCell,
} from "../components/task/task-interfaces/dm_multidomain/dm-spec";

describe("assembly rows", () => {
  it("reads a homomer as one implicit role", () => {
    const assembly = parseAssemblyRows(["A", "B", "C"]);
    expect(assembly.roles).toEqual([IMPLICIT_ROLE]);
    expect(assembly.instances).toEqual([
      { [IMPLICIT_ROLE]: "A" },
      { [IMPLICIT_ROLE]: "B" },
      { [IMPLICIT_ROLE]: "C" },
    ]);
    expect(isTerse(assembly.roles)).toBe(true);
  });

  it("reads role=chain rows, keeping column order", () => {
    const assembly = parseAssemblyRows(["CDK=A cyclin=B", "CDK=C cyclin=D"]);
    expect(assembly.roles).toEqual(["CDK", "cyclin"]);
    expect(assembly.instances[1]).toEqual({ CDK: "C", cyclin: "D" });
  });

  it("treats a missing role as a partial copy, not an error", () => {
    const assembly = parseAssemblyRows(["A=A B=B", "A=C"]);
    expect(assembly.instances[1]).toEqual({ A: "C" });
  });

  it("accepts commas as separators, as the wrapper does", () => {
    expect(parseAssemblyRows(["CDK=A,cyclin=B"]).instances[0]).toEqual({
      CDK: "A",
      cyclin: "B",
    });
  });

  it("keeps an empty row, so Add copy is not a dead button", () => {
    const assembly = parseAssemblyRows(["A", "", "C"]);
    expect(assembly.instances).toHaveLength(3);
    expect(assembly.instances[1]).toEqual({});
    expect(formatAssemblyRows(assembly)).toEqual(["A", "", "C"]);
  });

  it("round-trips a homomer back to the terse form", () => {
    const rows = ["A", "B", "C"];
    expect(formatAssemblyRows(parseAssemblyRows(rows))).toEqual(rows);
  });

  it("round-trips a hetero-complex", () => {
    const rows = ["CDK=A cyclin=B", "CDK=C cyclin=D"];
    expect(formatAssemblyRows(parseAssemblyRows(rows))).toEqual(rows);
  });

  it("keeps a named single entity explicit, or the name would be lost", () => {
    const assembly = parseAssemblyRows(["CDK=A", "CDK=B"]);
    expect(isTerse(assembly.roles)).toBe(false);
    expect(formatAssemblyRows(assembly)).toEqual(["CDK=A", "CDK=B"]);
  });

  it("labels the unnamed entity after the reference chain", () => {
    const assembly = parseAssemblyRows(["A", "B"]);
    expect(roleLabel(IMPLICIT_ROLE, assembly)).toBe("entity A");
    expect(roleLabel("CDK", assembly)).toBe("CDK");
  });
});

describe("segments", () => {
  it("reads a bare range as the implicit role", () => {
    expect(parseSegments("340-485")).toEqual([
      { role: IMPLICIT_ROLE, lo: 340, hi: 485 },
    ]);
  });

  it("reads a cross-entity body", () => {
    expect(parseSegments("cyclin:10-95,CDK:45-60")).toEqual([
      { role: "cyclin", lo: 10, hi: 95 },
      { role: "CDK", lo: 45, hi: 60 },
    ]);
  });

  it("refuses a spec it cannot read rather than returning an empty body", () => {
    expect(() => parseSegments("CDK:forty-two")).toThrow(SegmentParseError);
  });

  it("keeps an unreadable stored spec instead of silently dropping it", () => {
    const [body] = parseBodies([{ segments: "CDK:forty-two", mode: "average" }]);
    expect(body.error).toBeTruthy();
    expect(formatBodies([body], false)).toEqual([
      { segments: "CDK:forty-two", mode: "average" },
    ]);
  });

  it("falls back to average for an unknown mode", () => {
    expect(parseBodies([{ segments: "1-10", mode: "wat" }])[0].mode).toBe("average");
  });

  it("writes bare ranges for a homomer and roled ranges otherwise", () => {
    const segments = parseSegments("A:1-298");
    expect(formatBodies([{ segments, mode: "average" }], true)).toEqual([
      { segments: "1-298", mode: "average" },
    ]);
    expect(formatBodies([{ segments, mode: "average" }], false)).toEqual([
      { segments: "A:1-298", mode: "average" },
    ]);
  });
});

describe("edits that touch both parameters", () => {
  const model = (): Model => ({
    assembly: parseAssemblyRows(["A=A B=B", "A=C B=D"]),
    bodies: parseBodies([
      { segments: "A:1-298", mode: "average" },
      { segments: "A:45-60,B:175-200", mode: "refine" },
    ]),
  });

  it("renaming a role rewrites the bodies that use it", () => {
    const next = renameRole(model(), "A", "CDK");
    expect(formatAssemblyRows(next.assembly)).toEqual([
      "CDK=A B=B",
      "CDK=C B=D",
    ]);
    expect(formatBodies(next.bodies, false)).toEqual([
      { segments: "CDK:1-298", mode: "average" },
      { segments: "CDK:45-60,B:175-200", mode: "refine" },
    ]);
  });

  it("refuses a rename that would collide with another role", () => {
    const before = model();
    expect(renameRole(before, "A", "B")).toBe(before);
  });

  it("removing a role drops its ranges from every body", () => {
    const next = removeRole(model(), "B");
    expect(next.assembly.roles).toEqual(["A"]);
    expect(formatBodies(next.bodies, false)).toEqual([
      { segments: "A:1-298", mode: "average" },
      { segments: "A:45-60", mode: "refine" },
    ]);
  });

  it("adding a second entity names the implicit one first", () => {
    // A bare token cannot say which of two entities it means, so the terse
    // form has to be given up the moment a second entity appears -- in both
    // parameters at once.
    const terse: Model = {
      assembly: parseAssemblyRows(["A", "B"]),
      bodies: parseBodies([{ segments: "12-485", mode: "average" }]),
    };
    const next = addRole(terse, ["cyclin"]);
    expect(next.assembly.roles).toEqual(["A", "cyclin"]);
    expect(formatAssemblyRows(next.assembly)).toEqual(["A=A", "A=B"]);
    expect(formatBodies(next.bodies, isTerse(next.assembly.roles))).toEqual([
      { segments: "A:12-485", mode: "average" },
    ]);
  });

  it("does not lose the click when the new entity's name is taken", () => {
    // Naming the implicit role consumes a candidate: on a homomer of chains
    // A and B, "add entity" offers A first, which the rename has just used.
    const terse: Model = {
      assembly: parseAssemblyRows(["A", "B"]),
      bodies: parseBodies([{ segments: "12-485", mode: "average" }]),
    };
    const next = addRole(terse, ["A", "B"]);
    expect(next.assembly.roles).toEqual(["A", "B"]);
  });

  it("clearing a cell makes that copy partial rather than dropping the row", () => {
    const next = setCell(model(), 1, "B", "");
    expect(next.assembly.instances[1]).toEqual({ A: "C" });
    expect(formatAssemblyRows(next.assembly)).toEqual(["A=A B=B", "A=C"]);
  });
});

describe("what the strip has to show", () => {
  it("finds bodies claiming the same residues of the same entity", () => {
    const bodies = parseBodies([
      { segments: "A:1-298", mode: "average" },
      { segments: "A:1-100", mode: "average" },
      { segments: "B:1-100", mode: "average" },
    ]);
    expect(overlappingBodies(bodies)).toEqual([[0, 1]]);
  });

  it("bounds a range by the reference chain's own numbering", () => {
    const assembly = parseAssemblyRows(["A=A B=B", "A=C B=D"]);
    const chains = [
      { id: "A", first: 1, last: 298 },
      { id: "B", first: 175, last: 432 },
      { id: "C", first: 1, last: 298 },
    ];
    expect(referenceBounds(assembly, "A", chains)).toEqual({ lo: 1, hi: 298 });
    expect(referenceBounds(assembly, "B", chains)).toEqual({ lo: 175, hi: 432 });
    expect(referenceBounds(assembly, "nope", chains)).toBeNull();
  });
});
