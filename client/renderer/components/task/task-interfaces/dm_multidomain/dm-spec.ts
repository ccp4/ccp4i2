/**
 * The dm_multidomain parameter grammar, as data.
 *
 * ASSEMBLY and DOMAINS are stored as strings in two small ad-hoc languages
 * ("CDK=A cyclin=B" rows, "cyclin:10-95,CDK:45-60" segment specs) that the
 * wrapper parses in `dm_ncs_lib`. i2run and the CLI write them by hand, so the
 * format stays; the interface never asks anyone to type them. This module is
 * the one place that knows the grammar: everything above it works in terms of
 * instances, roles and segments, and serialises through here.
 *
 * Kept in step with server/ccp4i2/wrappers/dm_multidomain/script/dm_ncs_lib.py
 * (parse_segments / parse_assembly_rows / format_assembly_rows). The rules
 * that matter:
 *
 *  - a bare token means the single IMPLICIT role, so a homomer stays terse
 *    ("A", "B", "C" and "340-485") and never grows a vocabulary;
 *  - instance 0 is the reference copy, the one the masks are cut from;
 *  - a role missing from a row is a partial copy (the AmBn case), not an error.
 */

export const IMPLICIT_ROLE = "_";

export type Mode = "average" | "refine" | "exclude";

export const MODES: Mode[] = ["average", "refine", "exclude"];

/** What each mode does, in one line, for the UI to say next to the control. */
export const MODE_HELP: Record<Mode, string> = {
  average: "averaged across the copies with the operators as fitted",
  refine: "averaged, and dm refines the operators as it goes",
  exclude: "left out of averaging altogether",
};

/** One residue range of one entity: the atoms of `role` from `lo` to `hi`. */
export interface Segment {
  role: string;
  lo: number;
  hi: number;
}

/** One NCS copy: which chain plays each role. */
export type Instance = Record<string, string>;

/** One rigid body: the segments that move together, and how to treat them. */
export interface Body {
  segments: Segment[];
  mode: Mode;
}

export interface Assembly {
  /** Copies, in order. `instances[0]` is the reference. */
  instances: Instance[];
  /** Roles in column order (order of first appearance across the rows). */
  roles: string[];
}

export const EMPTY_ASSEMBLY: Assembly = { instances: [], roles: [] };

/**
 * True when the assembly has exactly one, unnamed entity — the case that is
 * written in the terse bare-chain / bare-range form. A single entity the user
 * has named stays explicit, or the name would be lost on the next write.
 */
export const isTerse = (roles: string[]): boolean =>
  roles.length === 1 && roles[0] === IMPLICIT_ROLE;

/** How a role is shown as a column heading. The implicit role has no name of
 *  its own, so it borrows the reference chain's. */
export const roleLabel = (role: string, assembly: Assembly): string => {
  if (role !== IMPLICIT_ROLE) return role;
  const chain = assembly.instances[0]?.[IMPLICIT_ROLE];
  return chain ? `entity ${chain}` : "entity";
};

// ---------------------------------------------------------------------------
// Assembly
// ---------------------------------------------------------------------------

export const parseAssemblyRows = (rows: string[]): Assembly => {
  const instances: Instance[] = [];
  const roles: string[] = [];
  for (const row of rows || []) {
    const instance: Instance = {};
    for (const token of String(row ?? "").replace(/,/g, " ").split(/\s+/)) {
      if (!token) continue;
      const eq = token.indexOf("=");
      const role = eq === -1 ? IMPLICIT_ROLE : token.slice(0, eq).trim() || IMPLICIT_ROLE;
      const chain = (eq === -1 ? token : token.slice(eq + 1)).trim();
      if (!chain) continue;
      instance[role] = chain;
      if (!roles.includes(role)) roles.push(role);
    }
    // An empty row is kept, not dropped: "Add copy" produces one, and a row
    // that vanished the moment it was added would be a dead button. The
    // wrapper skips empty rows when it parses them, so an unfilled copy costs
    // nothing but the space it occupies on screen.
    instances.push(instance);
  }
  return { instances, roles };
};

export const formatAssemblyRows = (assembly: Assembly): string[] => {
  const terse = isTerse(assembly.roles);
  return assembly.instances
    .map((instance) =>
      terse
        ? instance[IMPLICIT_ROLE] ?? ""
        : assembly.roles
            .filter((role) => instance[role])
            .map((role) => `${role}=${instance[role]}`)
            .join(" ")
    );
};

/** Chains already spoken for, so a cell's dropdown can avoid offering them. */
export const usedChains = (assembly: Assembly): Set<string> => {
  const used = new Set<string>();
  for (const instance of assembly.instances) {
    for (const chain of Object.values(instance)) used.add(chain);
  }
  return used;
};

// ---------------------------------------------------------------------------
// Segments and bodies
// ---------------------------------------------------------------------------

export class SegmentParseError extends Error {}

/**
 * "cyclin:10-95,CDK:45-60" -> segments. Throws SegmentParseError on anything
 * it cannot read, so a hand-edited or i2run-written value surfaces as a
 * message rather than as a silently empty body.
 */
export const parseSegments = (spec: string): Segment[] => {
  const out: Segment[] = [];
  for (const raw of String(spec ?? "").split(",")) {
    const token = raw.trim();
    if (!token) continue;
    const colon = token.indexOf(":");
    const role = colon === -1 ? IMPLICIT_ROLE : token.slice(0, colon).trim() || IMPLICIT_ROLE;
    const range = colon === -1 ? token : token.slice(colon + 1);
    const match = /^(-?\d+)\s*-\s*(-?\d+)$/.exec(range.trim());
    if (!match) {
      throw new SegmentParseError(`cannot read the residue range "${token}"`);
    }
    out.push({ role, lo: Number(match[1]), hi: Number(match[2]) });
  }
  return out;
};

export const formatSegments = (segments: Segment[], terse: boolean): string =>
  segments
    .map((s) => (terse || s.role === IMPLICIT_ROLE ? `${s.lo}-${s.hi}` : `${s.role}:${s.lo}-${s.hi}`))
    .join(",");

/** A DOMAINS row as it comes off the container. */
export interface RawBody {
  segments?: string;
  mode?: string;
}

export interface ParsedBody extends Body {
  /** Set when the stored spec could not be parsed; `segments` is then []. */
  error?: string;
  /** The unparseable spec, kept so the user can see what is actually stored. */
  raw?: string;
}

export const parseBodies = (rows: RawBody[]): ParsedBody[] =>
  (rows || []).map((row) => {
    const mode = (MODES as string[]).includes(String(row?.mode))
      ? (row.mode as Mode)
      : "average";
    const spec = String(row?.segments ?? "");
    try {
      return { segments: parseSegments(spec), mode };
    } catch (err) {
      return {
        segments: [],
        mode,
        error: err instanceof Error ? err.message : String(err),
        raw: spec,
      };
    }
  });

export const formatBodies = (bodies: ParsedBody[], terse: boolean): RawBody[] =>
  bodies.map((body) => ({
    segments: body.error ? body.raw ?? "" : formatSegments(body.segments, terse),
    mode: body.mode,
  }));

// ---------------------------------------------------------------------------
// Edits that touch both parameters
//
// Roles are the one thing the two lists share, so any edit to them has to
// rewrite both or the cross-reference breaks — which is the failure the old
// free-text interface made so easy.
// ---------------------------------------------------------------------------

export interface Model {
  assembly: Assembly;
  bodies: ParsedBody[];
}

export const renameRole = (model: Model, from: string, to: string): Model => {
  const name = to.trim();
  if (!name || name === from || model.assembly.roles.includes(name)) return model;
  return {
    assembly: {
      roles: model.assembly.roles.map((role) => (role === from ? name : role)),
      instances: model.assembly.instances.map((instance) => {
        const next: Instance = {};
        for (const [role, chain] of Object.entries(instance)) {
          next[role === from ? name : role] = chain;
        }
        return next;
      }),
    },
    bodies: model.bodies.map((body) => ({
      ...body,
      segments: body.segments.map((s) => (s.role === from ? { ...s, role: name } : s)),
    })),
  };
};

/**
 * Add an entity. The implicit role cannot survive a second entity — a bare
 * token has no way to say which of two things it means — so it is named after
 * the reference chain first, and both lists are rewritten to match.
 *
 * `candidates` are names to try in order (the model's chain ids, so the new
 * column is called after a chain the user can see). The naming of the implicit
 * role can itself consume a candidate, so the free name is chosen after that
 * rename, not before it.
 */
export const addRole = (model: Model, candidates: string[] = []): Model => {
  let next = model;
  if (isTerse(next.assembly.roles)) {
    const reference = next.assembly.instances[0]?.[IMPLICIT_ROLE];
    next = renameRole(next, IMPLICIT_ROLE, reference || "entity1");
  }
  const taken = new Set(next.assembly.roles);
  const name =
    candidates.find((candidate) => candidate && !taken.has(candidate)) ??
    `entity${next.assembly.roles.length + 1}`;
  if (taken.has(name)) return next;
  return { ...next, assembly: { ...next.assembly, roles: [...next.assembly.roles, name] } };
};

export const removeRole = (model: Model, role: string): Model => ({
  assembly: {
    roles: model.assembly.roles.filter((r) => r !== role),
    instances: model.assembly.instances.map((instance) => {
      const next = { ...instance };
      delete next[role];
      return next;
    }),
  },
  bodies: model.bodies.map((body) => ({
    ...body,
    segments: body.segments.filter((s) => s.role !== role),
  })),
});

export const setCell = (
  model: Model,
  instanceIndex: number,
  role: string,
  chain: string
): Model => ({
  ...model,
  assembly: {
    ...model.assembly,
    instances: model.assembly.instances.map((instance, i) => {
      if (i !== instanceIndex) return instance;
      const next = { ...instance };
      if (chain) next[role] = chain;
      else delete next[role];
      return next;
    }),
  },
});

export const addInstance = (model: Model): Model => ({
  ...model,
  assembly: { ...model.assembly, instances: [...model.assembly.instances, {}] },
});

export const removeInstance = (model: Model, index: number): Model => ({
  ...model,
  assembly: {
    ...model.assembly,
    instances: model.assembly.instances.filter((_, i) => i !== index),
  },
});

/**
 * The residues of a role in the reference copy, from the model's own chains.
 * Used to bound the range pickers and to draw the coverage strip, so the
 * interface never offers a residue that does not exist.
 */
export interface ChainInfo {
  id: string;
  first: number | null;
  last: number | null;
  nResidues?: number;
  entity?: number | null;
}

export const referenceBounds = (
  assembly: Assembly,
  role: string,
  chains: ChainInfo[]
): { lo: number; hi: number } | null => {
  const chainId = assembly.instances[0]?.[role];
  if (!chainId) return null;
  const chain = chains.find((c) => c.id === chainId);
  if (!chain || chain.first === null || chain.last === null) return null;
  return { lo: chain.first, hi: chain.last };
};

/** Bodies that claim the same residues of the same role — a typo, almost
 *  always, and the thing the coverage strip is there to make visible. */
export const overlappingBodies = (bodies: ParsedBody[]): Array<[number, number]> => {
  const clashes: Array<[number, number]> = [];
  for (let a = 0; a < bodies.length; a += 1) {
    for (let b = a + 1; b < bodies.length; b += 1) {
      const hit = bodies[a].segments.some((s) =>
        bodies[b].segments.some(
          (t) => s.role === t.role && s.lo <= t.hi && t.lo <= s.hi
        )
      );
      if (hit) clashes.push([a, b]);
    }
  }
  return clashes;
};
