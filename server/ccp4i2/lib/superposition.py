"""
Rigid-body superposition of one structure onto another, computed here
rather than in the viewer.

A scene builder that overlays several datasets has to put them in one
frame. Coot can do the fit, but then nothing about *which* atoms were
fitted -- the radius, the count, the residues dropped as outliers -- can
be tested or reported. So the fit is done here with gemmi, as a pure
function over two structures and an optional point, and the scene carries
the resulting matrix together with the provenance it was derived from
(``superpose: method: matrix``, see ``docs/campaign-site-scene-design.md``).

The transform is ``x' = mat . x + vec`` about the origin, mapping the
*moving* structure onto the *reference*. That direction is the one easy
thing to get silently wrong: ``gemmi.superpose_positions(pos1, pos2)``
returns the transform taking ``pos2`` onto ``pos1`` (verified on gemmi
0.7.5, and pinned by ``test_round_trip_recovers_the_inverse``), so the
reference goes first.

This module must stay CCP4-free: it runs on the request path of a slim
server that has gemmi from pip and nothing else.
"""
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import gemmi

# A local fit starts at a radius that, in a folded domain, encloses several
# dozen residues -- more than a stable fit needs -- and grows in steps for
# the sparse cases (a shallow surface site, a small protein, a domain edge).
# Past the cap it is no longer a *local* fit in any meaningful sense, so the
# growth stops and the caller is told, rather than a global fit being
# reported as a local one.
FIT_RADIUS_START = 15.0
FIT_RADIUS_STEP = 5.0
FIT_RADIUS_MAX = 30.0

# Three non-collinear points determine a rigid transform in principle; that
# is not a usable floor. Twelve is a floor with margin, not a target.
MIN_FIT_CAS = 12

# Outlier rejection: fit, drop CAs further than this multiple of the RMSD
# from their partner, refit. Two rounds is enough to shed a rim loop that
# genuinely moved without eroding the set on noise; the floor stops a
# near-perfect fit (RMSD ~0) from rejecting everything but exact matches.
OUTLIER_FACTOR = 2.0
OUTLIER_MIN_CUTOFF = 0.5
OUTLIER_ROUNDS = 2

CaKey = Tuple[str, str]  # (chain name, seqid with insertion code)


@dataclass
class FitResult:
    """Outcome of :func:`fit_structures`.

    ``ok`` is False when no trustworthy fit could be made; ``reason`` says
    why in words a ``stats`` block can carry. ``atoms`` is the number of CAs
    in the *final* fit, after outlier rejection, so it is honest about what
    the matrix rests on. ``radius`` is None for a global fit.
    """

    ok: bool
    mat: Optional[List[float]] = None  # 9 numbers, row-major
    vec: Optional[List[float]] = None
    atoms: int = 0
    radius: Optional[float] = None
    rmsd: Optional[float] = None
    reason: Optional[str] = None

    def transform(self) -> gemmi.Transform:
        """The fit as a gemmi transform, so a caller can measure in the
        frame the scene draws the moving structure in."""
        if not self.ok:
            raise ValueError(f"no fit to apply: {self.reason}")
        return gemmi.Transform(
            gemmi.Mat33([self.mat[0:3], self.mat[3:6], self.mat[6:9]]),
            gemmi.Vec3(*self.vec),
        )

    def provenance(self, onto: str) -> dict:
        """The ``fitted`` block of a ``matrix`` superpose entry."""
        out = {"onto": onto, "atoms": self.atoms, "rmsd": self.rmsd}
        if self.radius is not None:
            out["radius"] = self.radius
        return out

    def scene_entry(self, move: str, onto: str) -> dict:
        """A ``superpose[]`` entry applying this fit to file ``move``."""
        if not self.ok:
            raise ValueError(f"no fit to emit: {self.reason}")
        return {
            "method": "matrix",
            "move": move,
            "mat": self.mat,
            "vec": self.vec,
            "fitted": self.provenance(onto),
        }


def ca_positions(structure: gemmi.Structure) -> Dict[CaKey, gemmi.Position]:
    """CA atoms of the first model, keyed by (chain, seqid).

    Only amino-acid residues count, and only a carbon called CA: a calcium
    ion is a residue named CA with an atom named CA, and a fit that quietly
    included ions as backbone would be wrong in a way nothing reports.
    """
    out: Dict[CaKey, gemmi.Position] = {}
    if len(structure) == 0:
        return out
    for chain in structure[0]:
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if info is None or not info.is_amino_acid():
                continue
            atom = res.find_atom("CA", "*", gemmi.Element("C"))
            if atom is None:
                continue
            out[(chain.name, str(res.seqid))] = atom.pos
    return out


def _paired(reference, moving, keys):
    """Positions of ``keys`` from each side, in one order."""
    return [reference[k] for k in keys], [moving[k] for k in keys]


def _fit_with_rejection(ref_pos, mov_pos) -> Tuple[gemmi.SupResult, int]:
    """Superpose, shedding outliers, and return the result with its count.

    A rim loop that moved between datasets would otherwise drag the whole
    frame; dropping the worst-fitting CAs and refitting is what makes a
    local fit robust. The count gate applies to the retained set too: a
    round that would leave too few CAs is not taken, and the previous fit
    stands.
    """
    keep = list(range(len(ref_pos)))
    result = gemmi.superpose_positions(ref_pos, mov_pos)
    for _ in range(OUTLIER_ROUNDS):
        cutoff = max(OUTLIER_FACTOR * result.rmsd, OUTLIER_MIN_CUTOFF)
        retained = [
            i for i in keep
            if gemmi.Position(result.transform.apply(mov_pos[i])).dist(ref_pos[i]) <= cutoff
        ]
        if len(retained) == len(keep) or len(retained) < MIN_FIT_CAS:
            break
        keep = retained
        result = gemmi.superpose_positions(
            [ref_pos[i] for i in keep], [mov_pos[i] for i in keep]
        )
    return result, len(keep)


def fit_structures(
    reference: gemmi.Structure,
    moving: gemmi.Structure,
    centre: Optional[Tuple[float, float, float]] = None,
) -> FitResult:
    """Fit ``moving`` onto ``reference`` on the CAs they share.

    With no ``centre`` the fit is global: every residue that has a CA in
    both structures. With one it is local: the reference's CAs within a
    sphere about that point, grown from ``FIT_RADIUS_START`` until the
    *shared* set reaches ``MIN_FIT_CAS`` or the cap is hit.

    Residues are matched by (chain, seqid), so numbering must correspond.
    When it does not, the intersection collapses and the count gate trips,
    which is the safe failure: a confident fit on coincidentally numbered
    residues would be far worse than none.
    """
    ref_cas = ca_positions(reference)
    mov_cas = ca_positions(moving)
    shared = [k for k in ref_cas if k in mov_cas]

    if centre is None:
        if len(shared) < MIN_FIT_CAS:
            return FitResult(
                ok=False,
                atoms=len(shared),
                reason=f"only {len(shared)} CA atoms in common (need {MIN_FIT_CAS})",
            )
        keys, radius = shared, None
    else:
        origin = gemmi.Position(*centre)
        radius = FIT_RADIUS_START
        while True:
            keys = [k for k in shared if ref_cas[k].dist(origin) <= radius]
            if len(keys) >= MIN_FIT_CAS:
                break
            if radius >= FIT_RADIUS_MAX:
                return FitResult(
                    ok=False,
                    atoms=len(keys),
                    radius=radius,
                    reason=(
                        f"only {len(keys)} shared CA atoms within {radius:g} A "
                        f"of the site (need {MIN_FIT_CAS})"
                    ),
                )
            radius = min(radius + FIT_RADIUS_STEP, FIT_RADIUS_MAX)

    ref_pos, mov_pos = _paired(ref_cas, mov_cas, keys)
    result, retained = _fit_with_rejection(ref_pos, mov_pos)
    mat = result.transform.mat.tolist()
    return FitResult(
        ok=True,
        mat=[float(v) for row in mat for v in row],
        vec=[float(v) for v in result.transform.vec.tolist()],
        atoms=retained,
        radius=radius,
        rmsd=float(result.rmsd),
    )


def fit_files(reference_path, moving_path, centre=None) -> FitResult:
    """:func:`fit_structures` over two coordinate files.

    An unreadable file is a failed fit, not an exception: the scene still
    draws the dataset, untransformed, and says why.
    """
    try:
        reference = gemmi.read_structure(str(reference_path))
        moving = gemmi.read_structure(str(moving_path))
    except Exception as exc:  # noqa: BLE001 - any read failure is "cannot fit"
        return FitResult(ok=False, reason=f"could not read coordinates: {exc}")
    return fit_structures(reference, moving, centre)
