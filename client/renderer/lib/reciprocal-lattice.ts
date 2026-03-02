/**
 * Reciprocal lattice math utilities for crystallographic visualization.
 *
 * Provides functions for computing reciprocal lattice vectors from unit cell
 * parameters, projecting Miller indices into Cartesian reciprocal-space
 * coordinates, and indexing reflection data by section plane.
 */

export type Vec3 = [number, number, number];

export type SectionPlane = "hk" | "hl" | "kl";

export interface ReciprocalVectors {
  aStar: Vec3;
  bStar: Vec3;
  cStar: Vec3;
}

export interface SectionPoint {
  /** Cartesian x in reciprocal space (Å⁻¹) */
  u: number;
  /** Cartesian y in reciprocal space (Å⁻¹) */
  v: number;
  /** Miller indices */
  h: number;
  k: number;
  l: number;
  intensity?: number;
  /** True if this reflection is from the asymmetric unit (not symmetry-expanded) */
  isASU?: boolean;
}

/**
 * 2D basis for a reciprocal-space section plane, computed via Gram-Schmidt
 * from the two free reciprocal lattice vectors.
 */
export interface SectionBasis {
  /** Projected coordinates of the u-axis reciprocal vector in 2D */
  uProj: [number, number];
  /** Projected coordinates of the v-axis reciprocal vector in 2D */
  vProj: [number, number];
  /** Orthonormal basis vector 1 (3D unit vector) */
  e1: Vec3;
  /** Orthonormal basis vector 2 (3D unit vector) */
  e2: Vec3;
}

export function degreesToRadians(degrees: number): number {
  return (degrees * Math.PI) / 180;
}

export function crossProduct(v1: Vec3, v2: Vec3): Vec3 {
  return [
    v1[1] * v2[2] - v1[2] * v2[1],
    v1[2] * v2[0] - v1[0] * v2[2],
    v1[0] * v2[1] - v1[1] * v2[0],
  ];
}

export function vectorMagnitude(v: Vec3): number {
  return Math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2);
}

function dot3(a: Vec3, b: Vec3): number {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/**
 * Compute reciprocal lattice vectors from unit cell parameters.
 *
 * @param a,b,c - Cell edge lengths in Angstroms
 * @param alpha,beta,gamma - Cell angles in degrees
 */
export function reciprocalLatticeVectors(
  a: number,
  b: number,
  c: number,
  alpha: number,
  beta: number,
  gamma: number
): ReciprocalVectors {
  const alphaRad = degreesToRadians(alpha);
  const betaRad = degreesToRadians(beta);
  const gammaRad = degreesToRadians(gamma);

  const aVec: Vec3 = [a, 0, 0];
  const bVec: Vec3 = [b * Math.cos(gammaRad), b * Math.sin(gammaRad), 0];
  const cVec: Vec3 = [
    c * Math.cos(betaRad),
    (c * (Math.cos(alphaRad) - Math.cos(betaRad) * Math.cos(gammaRad))) /
      Math.sin(gammaRad),
    (c *
      Math.sqrt(
        1 -
          Math.cos(alphaRad) ** 2 -
          Math.cos(betaRad) ** 2 -
          Math.cos(gammaRad) ** 2 +
          2 *
            Math.cos(alphaRad) *
            Math.cos(betaRad) *
            Math.cos(gammaRad)
      )) /
      Math.sin(gammaRad),
  ];

  const volume =
    aVec[0] * (bVec[1] * cVec[2] - bVec[2] * cVec[1]) -
    aVec[1] * (bVec[0] * cVec[2] - bVec[2] * cVec[0]) +
    aVec[2] * (bVec[0] * cVec[1] - bVec[1] * cVec[0]);

  const aStar: Vec3 = crossProduct(bVec, cVec).map((v) => v / volume) as Vec3;
  const bStar: Vec3 = crossProduct(cVec, aVec).map((v) => v / volume) as Vec3;
  const cStar: Vec3 = crossProduct(aVec, bVec).map((v) => v / volume) as Vec3;

  return { aStar, bStar, cStar };
}

/**
 * Get axis labels for a section plane.
 * Returns [horizontalAxis, verticalAxis, fixedAxis].
 */
export function getSectionAxes(
  plane: SectionPlane
): ["h" | "k" | "l", "h" | "k" | "l", "h" | "k" | "l"] {
  switch (plane) {
    case "hk":
      return ["h", "k", "l"];
    case "hl":
      return ["h", "l", "k"];
    case "kl":
      return ["k", "l", "h"];
  }
}

/**
 * Get the two free reciprocal vectors for a section plane.
 */
function getSectionVectors(
  recipVecs: ReciprocalVectors,
  plane: SectionPlane
): { uVec: Vec3; vVec: Vec3 } {
  switch (plane) {
    case "hk":
      return { uVec: recipVecs.aStar, vVec: recipVecs.bStar };
    case "hl":
      return { uVec: recipVecs.aStar, vVec: recipVecs.cStar };
    case "kl":
      return { uVec: recipVecs.bStar, vVec: recipVecs.cStar };
  }
}

/**
 * Compute a 2D orthonormal basis for a section plane in reciprocal space
 * using Gram-Schmidt on the two free reciprocal lattice vectors.
 *
 * The first basis vector e1 is along the u-axis reciprocal vector.
 * The second e2 is the component of the v-axis vector perpendicular to e1.
 *
 * Returns the basis vectors and the projected coordinates of the reciprocal
 * vectors onto this basis (used for drawing the lattice grid).
 */
export function computeSectionBasis(
  recipVecs: ReciprocalVectors,
  plane: SectionPlane
): SectionBasis {
  const { uVec, vVec } = getSectionVectors(recipVecs, plane);

  // e1 = unit(uVec)
  const e1Mag = vectorMagnitude(uVec);
  const e1: Vec3 = [uVec[0] / e1Mag, uVec[1] / e1Mag, uVec[2] / e1Mag];

  // e2 = unit(vVec - (vVec·e1)e1)
  const vDotE1 = dot3(vVec, e1);
  const e2Raw: Vec3 = [
    vVec[0] - vDotE1 * e1[0],
    vVec[1] - vDotE1 * e1[1],
    vVec[2] - vDotE1 * e1[2],
  ];
  const e2Mag = vectorMagnitude(e2Raw);
  const e2: Vec3 = [e2Raw[0] / e2Mag, e2Raw[1] / e2Mag, e2Raw[2] / e2Mag];

  // Project the reciprocal vectors onto the 2D basis
  const uProj: [number, number] = [dot3(uVec, e1), dot3(uVec, e2)];
  const vProj: [number, number] = [dot3(vVec, e1), dot3(vVec, e2)];

  return { uProj, vProj, e1, e2 };
}

/**
 * Index reflections by the fixed axis value for a given section plane,
 * projecting Miller indices into Cartesian reciprocal-space coordinates
 * using the section basis.
 *
 * The resulting u,v coordinates are in Å⁻¹ and reflect the true geometry
 * of the reciprocal lattice — distance from origin = 1/d.
 */
export function indexReflectionsBySection(
  reflections: { h: number; k: number; l: number; intensity?: number; isASU?: boolean }[],
  plane: SectionPlane,
  basis: SectionBasis
): Map<number, SectionPoint[]> {
  const map = new Map<number, SectionPoint[]>();
  const [uAxis, vAxis, fixedAxis] = getSectionAxes(plane);

  for (const ref of reflections) {
    const fixedVal = ref[fixedAxis];
    const uIdx = ref[uAxis];
    const vIdx = ref[vAxis];

    // Project to Cartesian reciprocal space:
    // pos_2d = uIdx * uProj + vIdx * vProj
    const point: SectionPoint = {
      u: uIdx * basis.uProj[0] + vIdx * basis.vProj[0],
      v: uIdx * basis.uProj[1] + vIdx * basis.vProj[1],
      h: ref.h,
      k: ref.k,
      l: ref.l,
      intensity: ref.intensity,
      isASU: ref.isASU,
    };

    const existing = map.get(fixedVal);
    if (existing) {
      existing.push(point);
    } else {
      map.set(fixedVal, [point]);
    }
  }

  return map;
}

/**
 * Get the range of the fixed axis for a given section plane.
 */
export function getFixedAxisRange(
  hRange: [number, number],
  kRange: [number, number],
  lRange: [number, number],
  plane: SectionPlane
): [number, number] {
  switch (plane) {
    case "hk":
      return lRange;
    case "hl":
      return kRange;
    case "kl":
      return hRange;
  }
}

// ---------------------------------------------------------------------------
// Symmetry expansion
// ---------------------------------------------------------------------------

type Mat3 = [[number, number, number], [number, number, number], [number, number, number]];

/**
 * Parse a crystallographic symmetry operator string (e.g. "-Y,X-Y,Z+2/3")
 * and return the 3×3 rotation matrix **transposed** for reciprocal-space use.
 *
 * In reciprocal space h' = Rᵀ h, so we transpose the real-space rotation.
 * Translations are irrelevant for Miller indices and are ignored.
 */
function parseSymopRotation(symop: string): Mat3 {
  const parts = symop.toUpperCase().split(",");
  if (parts.length !== 3) {
    return [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  }

  // Parse one component string (e.g. "-Y", "X-Y", "1/2+X", "Z+2/3")
  // into coefficients for [X, Y, Z].
  function parseComponent(s: string): [number, number, number] {
    const coeffs: [number, number, number] = [0, 0, 0];
    const axisIdx = { X: 0, Y: 1, Z: 2 } as const;

    // Remove whitespace
    s = s.replace(/\s+/g, "");

    // Walk through the string finding X, Y, Z terms
    // We need to handle patterns like: X, -X, X-Y, -Y+X, 1/2+X, X+1/2
    let i = 0;
    while (i < s.length) {
      const ch = s[i];

      if (ch === "X" || ch === "Y" || ch === "Z") {
        // Found an axis variable — determine its sign from preceding character
        let sign = 1;
        if (i > 0 && s[i - 1] === "-") sign = -1;
        coeffs[axisIdx[ch]] = sign;
        i++;
      } else {
        // Skip non-axis characters (digits, fractions, +, -)
        i++;
      }
    }

    return coeffs;
  }

  // Build the real-space rotation matrix (rows)
  const R: Mat3 = [
    parseComponent(parts[0]),
    parseComponent(parts[1]),
    parseComponent(parts[2]),
  ];

  // Transpose for reciprocal space
  return [
    [R[0][0], R[1][0], R[2][0]],
    [R[0][1], R[1][1], R[2][1]],
    [R[0][2], R[1][2], R[2][2]],
  ];
}

/**
 * Expand reflections from the asymmetric unit to the full reciprocal sphere
 * by applying all space group symmetry operators plus Friedel's law.
 *
 * @param reflections - ASU reflections with Miller indices and optional intensity
 * @param symops - Transformation strings from spacegroups.ts (e.g. ["X,Y,Z", "-Y,X-Y,Z+2/3"])
 * @returns Expanded reflection array with duplicates removed
 */
export function expandReflectionsBySymmetry(
  reflections: { h: number; k: number; l: number; intensity?: number }[],
  symops: string[]
): { h: number; k: number; l: number; intensity?: number; isASU?: boolean }[] {
  // Pre-parse all rotation matrices
  const matrices = symops.map(parseSymopRotation);

  const seen = new Set<string>();
  const result: { h: number; k: number; l: number; intensity?: number; isASU?: boolean }[] = [];

  function addIfNew(h: number, k: number, l: number, intensity: number | undefined, isASU: boolean) {
    const key = `${h},${k},${l}`;
    if (!seen.has(key)) {
      seen.add(key);
      result.push({ h, k, l, intensity, isASU });
    }
  }

  for (const ref of reflections) {
    for (let i = 0; i < matrices.length; i++) {
      const RT = matrices[i];
      const isIdentity = i === 0; // First symop is always identity (X,Y,Z)

      // h' = RT * h
      const hp = Math.round(RT[0][0] * ref.h + RT[0][1] * ref.k + RT[0][2] * ref.l);
      const kp = Math.round(RT[1][0] * ref.h + RT[1][1] * ref.k + RT[1][2] * ref.l);
      const lp = Math.round(RT[2][0] * ref.h + RT[2][1] * ref.k + RT[2][2] * ref.l);

      addIfNew(hp, kp, lp, ref.intensity, isIdentity);
      // Friedel mate
      addIfNew(-hp, -kp, -lp, ref.intensity, false);
    }
  }

  return result;
}
