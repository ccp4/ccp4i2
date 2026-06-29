<!-- GENERATED from the Zod schema via buildSceneMarkdown(); do not edit by hand.
     The Grammar section below is authoritative. Edit lib/scene/core.ts (the Zod
     source) and regenerate (UPDATE_SCHEMA=1 vitest run scene-schema). -->
# Moorhen Scene format

A Moorhen **scene** is a portable YAML description of how to render one or more
structures: which files to load, which domains to recognise, what
representations and colours to apply, and where the camera sits. Scenes validate
against the published contract (`moorhen-scene.core.v1.json`, or
`moorhen-scene.ccp4i2.v1.json` for project references) and can be re-applied
across different structures of the same protein.

## Grammar

Fields marked `?` are optional. `A | B` is a choice; `T[]` is a list of `T`.

```yaml
Moorhen scene — YAML grammar (fields, ? = optional):

scene: string  # human-readable scene name
version: number  # must equal 1
authoredIn?: {
  projectId?: string
  projectName?: string
  createdAt?: string
  createdBy?: string
  ccp4i2Version?: string
}
files?: ({
  name: string  # local name referenced by elements/maps
  kind?: "coordinates"|"dictionary"|"mtz"|"map"  # default "coordinates"
  pdb?: string  # PDB id; fetched via proxy on apply
  url?: string  # absolute URL (portable)
  bundle?: string  # asset path inside a .scene.zip
  cifText?: string  # inline CIF (dictionary refs only)
  relativeUrl?: string  # origin-relative URL (/api/…); not portable across deployments
  projectId?: string  # UUID; required with fileId or job+param
  projectName?: string  # advisory
  fileId?: number
  job?: number  # pair with param
  param?: string  # job parameter, e.g. "XYZOUT"
})[]
superpose?: ({
  method: "ssm"
  move: string  # file being transformed
  onto: string  # reference file (unchanged)
  movChain: string
  refChain: string
} | {
  method: "lsq"
  move: string
  onto: string
  matches?: {
    refChain: string
    refRange: string
    movChain: string
    movRange: string
  }[]
  chain?: string  # shorthand chain, both sides
  range?: string  # shorthand range, both sides
  matchType?: "all"|"main"|"ca"  # default "main"
})[]
globalDictionaries?: string[]
domains?: { name: string, selection: string, color: string }[]
elements?: ({
  file: string  # name of a files[] entry
  dictionaries?: string[]
  colour?: string | "by-domain"|"b-factor"|"b-factor-norm"|"af2-plddt"|"secondary-structure"|"jones-rainbow"|"mol-symm" | { selection: string, colour: string }[] | {
    raw: {
      ruleType: string
      args: (string | number)[]
      isMultiColourRule?: boolean
      applyColourToNonCarbonAtoms?: boolean
    }
  }  # molecule-scoped colour: the default for every representation of this file; a representation's own `colour` overrides it
  representations?: ({
    style: "VdwSpheres"|"ligands"|"CAs"|"CBs"|"CDs"|"gaussian"|"allHBonds"|"rama"|"rotamer"|"CRs"|"MolecularSurface"|"DishyBases"|"VdWSurface"|"Calpha"|"unitCell"|"hover"|"environment"|"ligand_environment"|"contact_dots"|"chemical_features"|"ligand_validation"|"glycoBlocks"|"restraints"|"residueSelection"|"MetaBalls"|"adaptativeBonds"|"StickBases"|"residue_environment"|"transformation"  # Moorhen RepresentationStyle, e.g. "CRs", "CBs"
    selection?: string  # CID; default all-atoms
    colour?: string | "by-domain"|"b-factor"|"b-factor-norm"|"af2-plddt"|"secondary-structure"|"jones-rainbow"|"mol-symm" | { selection: string, colour: string }[] | {
      raw: {
        ruleType: string
        args: (string | number)[]
        isMultiColourRule?: boolean
        applyColourToNonCarbonAtoms?: boolean
      }
    }
    alpha?: number  # opacity 0..1 (honoured: governs visibility)
    geometry?: {
      bondRadius?: number  # bond cylinder radius (Å)
      ballRadius?: number  # ball-and-stick atom ball radius (Å)
      vdwScale?: number  # VdW sphere radius multiplier (×element radius)
      probeRadius?: number  # molecular-surface solvent probe radius (Å)
      ribbonCoilThickness?: number  # ribbon coil thickness (Å)
      ribbonHelixWidth?: number  # ribbon helix width (Å)
      ribbonStrandWidth?: number  # ribbon strand width (Å)
      ribbonArrowWidth?: number  # ribbon arrow width (Å)
      ribbonDNARNAWidth?: number  # nucleotide ribbon width (Å)
    }
  })[]
})[]
maps?: ({
  name: string
  file: string  # name of a files[] entry (kind mtz or map)
  columns?: {
    F?: string
    PHI?: string
    Fobs?: string
    SigFobs?: string
    FreeR?: string
    useWeight?: boolean
    calcStructFact?: boolean
  }  # required for mtz, omit for map
  isMask?: boolean
  isDifference?: boolean
  contourLevel?: number  # rmsd-relative
  radius?: number  # contour radius (Å)
  alpha?: number
  style?: "lines"|"solid"|"lit-lines"
  colour?: string  # non-difference maps only
  positiveColour?: string  # hex colour #rrggbb or #rrggbbaa
  negativeColour?: string  # hex colour #rrggbb or #rrggbbaa
  visible?: boolean
})[]
activeMap?: string
view?: {
  origin?: [number, number, number]
  centre?: { file?: string, selection?: string }  # centroid of a selection; beats origin
  quat?: [number, number, number, number]
  zoom?: number
  clipStart?: number
  clipEnd?: number
  fogStart?: number
  fogEnd?: number
  clip?: "auto" | "lock" | { front: number, back: number }
  slab?: { file?: string, selection?: string, pad?: number }  # z-depth window for a selection; beats clip
  background?: string  # hex colour #rrggbb or #rrggbbaa
}
hints?: {
  lighting?: {
    direction?: [number, number, number]  # principal directional light vector → Moorhen lightPosition
    ambient?: string  # ambient light colour
    diffuse?: string  # diffuse light colour
    specular?: string  # specular light colour
    shininess?: number  # specular power → Moorhen specularPower
  }
  effects?: {
    ssao?: { enabled?: boolean, radius?: number, bias?: number }  # screen-space ambient occlusion
    edgeDetect?: {
      enabled?: boolean
      depthThreshold?: number
      normalThreshold?: number
      depthScale?: number
      normalScale?: number
    }
    shadows?: boolean
    depthBlur?: { radius?: number, depth?: number }  # depth-of-field blur
    perspective?: boolean  # perspective vs orthographic
  }
}
resolver?: { onMissingResidues?: "clamp-and-log"|"strict" }
```

## File references

Each `files[]` entry sets **exactly one** source:

| Ref | Meaning | Portable? |
|-----|---------|-----------|
| `pdb` | PDB id, fetched via the PDBe proxy | yes (core) |
| `url` | absolute URL | yes (core) |
| `bundle` | asset inside a `.scene.zip` | yes (core) |
| `cifText` | inline dictionary CIF (`kind: dictionary` only) | yes (core) |
| `relativeUrl` | origin-relative URL (`/api/…`) | within one deployment |
| `fileId` (+`projectId`) | ccp4i2 project file | within one deployment |
| `job`+`param` (+`projectId`) | ccp4i2 job output | within one deployment |

The resolver matches a ref against an already-loaded molecule by its loader URL
(or by `fileId` extracted from it); a `pdb`/`url`/`bundle`/`job+param` ref is
fetched if not present. A molecule's identity is its loader URL, never a bare
path — so `relativeUrl` matches loaded molecules, it does not read local files.

## Selections (Coot CIDs)

Selections everywhere (representation `selection`, `view.centre`/`slab`,
domains) are Coot CIDs:

- `//A` — the whole of chain A
- `//A/703-740` — residues 703–740 of chain A
- `//*/LIG` — every residue named LIG
- `//A/750/CA` — one atom
- join several with `||`: `//A||//B`

A representation draws **its own** `selection` (the whole molecule if omitted);
colour does **not** limit what is drawn — scope the selection to limit it.

## Colours

A `colour` is one of:

1. a hex string — `"#4b8bbe"` (or `#rrggbbaa`);
2. a named scheme — `by-domain`, `b-factor`, `b-factor-norm`, `af2-plddt`,
   `secondary-structure`, `jones-rainbow`, `mol-symm`;
3. a per-selection list — `[{ selection: "//A", colour: "#a08766" }, …]`;
4. a raw Coot rule (escape hatch).

`by-domain` colours by the top-level `domains:` block: define domains
(`name` + `selection` + `color`) and set `colour: by-domain` on the reps.

Colour cascades over two levels (matching Moorhen's molecule-level colour):
`element.colour` is the molecule-wide default inherited by every representation
of that file, and a representation's own `colour` overrides it for that rep.

## Geometry and hints

- `geometry` (per representation) is **honoured** — physical dimensions in
  Ångström a renderer must reproduce (bond/ball radii, ribbon widths, …).
- `hints` (scene-level lighting + effects) are **advisory** — a renderer may
  ignore them and still produce a correct image. `hints.lighting.direction` is
  a direction vector (+x right, +y up, +z toward the viewer); effects are
  scene-authoritative (an `effects` block fully determines effect state).

## Examples

Minimal portable scene:

```yaml
scene: minimal
version: 1
files:
  - { name: prot, pdb: 1ABC }
elements:
  - file: prot
    representations:
      - { style: CRs, selection: "//A", colour: secondary-structure }
```

Domains, honoured geometry, advisory hints (validates under the core schema):

```yaml
scene: gamma-demo
version: 1
files:
  - { name: gamma, pdb: 1B9K }
domains:
  - { name: head,   selection: "//A/703-740",  color: "#4b8bbe" }
  - { name: tail,   selection: "//A/900-1000", color: "#e74c3c" }  # clamps to 938
elements:
  - file: gamma
    representations:
      - { style: CRs, selection: "//A", colour: by-domain }
      - { style: CBs, selection: "//A/750", geometry: { bondRadius: 0.18 } }
hints:
  lighting: { direction: [1, 1, 1], shininess: 24 }
  effects:  { ssao: { enabled: true }, edgeDetect: { enabled: true } }
resolver:
  onMissingResidues: clamp-and-log
```
