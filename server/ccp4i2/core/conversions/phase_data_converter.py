"""
Converter for phase data format transformations.

Conversions between phase data formats in MTZ files:

Phase Data Files (CPhsDataFile):
- HL (1): Hendrickson-Lattman coefficients (HLA, HLB, HLC, HLD)
- PHIFOM (2): Phase + Figure of Merit (PHI, FOM)

Map Coefficients Files (CMapCoeffsDataFile):
- FPHI (1): Structure factors + Phase (F, PHI)

Conversion Matrix (CPhsDataFile):
                TO
          HL   PHIFOM
FROM HL    ✓      ✓
     PHIFOM ✓      ✓

Implementation:
- Pure gemmi + numpy; no CCP4 binaries (no chltofom) and no clipper. This lets the
  conversion run on a slim, CCP4-free interpreter as well as inside CCP4.
- The numerical routines are a direct port of clipper's Compute_phifom_from_abcd
  and Compute_abcd_from_phifom (clipper/core/hkl_compute.cpp), including its
  centric/acentric handling, 5-degree phase grid and FOM clamp, so the results
  match the CCP4 chltofom tool (which is clipper) to within rounding. Centricity
  and the allowed centric phase are obtained from the space group via gemmi,
  reproducing clipper's HKL_class (clipper/core/coords.cpp).

Numerical agreement with chltofom is checked by
tests/parity/test_phase_conversions_chltofom_parity.py (run under ccp4-python).

References:
- Read, R.J. (1986). Acta Cryst. A42, 140-149. (Phase probability distributions)
- Hendrickson & Lattman (1970). Acta Cryst. B26, 136-143. (HL coefficients)
- clipper (K. Cowtan): core/hkl_compute.cpp, core/coords.cpp, core/clipper_util.cpp
"""

from typing import Optional, Any

import gemmi
import numpy as np

from ccp4i2.core.CCP4ErrorHandling import CException, SEVERITY_ERROR


class PhaseDataConverter:
    """
    Static converter class for phase data format transformations.

    Handles conversions between HL coefficients, PHI/FOM, and F/PHI formats.
    All methods are static and take the file instance as first parameter.
    """

    # Error codes for phase data conversions
    ERROR_CODES = {
        1: {'description': 'Input file does not exist', 'severity': SEVERITY_ERROR},
        2: {'description': 'Cannot determine contentFlag from input file', 'severity': SEVERITY_ERROR},
        3: {'description': 'Unsupported conversion: only HL and PHIFOM formats supported', 'severity': SEVERITY_ERROR},
        8: {'description': 'PHI column not found in input MTZ file', 'severity': SEVERITY_ERROR},
        9: {'description': 'FOM column not found in input MTZ file', 'severity': SEVERITY_ERROR},
        10: {'description': 'Invalid number of HL coefficient columns', 'severity': SEVERITY_ERROR},
    }

    @staticmethod
    def _validate_input_file(phase_file):
        """
        Validate input file before conversion.

        Raises:
            CException: If file doesn't exist or contentFlag cannot be determined
        """
        from pathlib import Path

        input_path = phase_file.getFullPath()
        if not Path(input_path).exists():
            raise CException(PhaseDataConverter, 1, details=f"File: {input_path}")

        content_flag = int(phase_file.contentFlag) if phase_file.contentFlag.isSet() else 0
        if content_flag == 0:
            raise CException(PhaseDataConverter, 2, details=f"File: {input_path}")

    @staticmethod
    def to_hl(phase_file, work_directory: Optional[Any] = None) -> str:
        """
        Convert phase data to HL format (HLA, HLB, HLC, HLD).

        PHIFOM -> HL is computed natively (gemmi/numpy). HL input is returned as-is.

        Raises:
            CException: If validation fails or conversion not supported
        """
        PhaseDataConverter._validate_input_file(phase_file)

        output_path = phase_file._get_conversion_output_path('HL', work_directory=work_directory)
        input_path = phase_file.getFullPath()
        content_flag = int(phase_file.contentFlag)

        if content_flag == 1:  # already HL
            return input_path
        elif content_flag == 2:  # PHIFOM -> HL
            return PhaseDataConverter._phifom_to_hl(input_path, output_path)
        else:
            raise CException(
                PhaseDataConverter, 3,
                details=f"contentFlag={content_flag}, target=HL. Only PHIFOM (2) -> HL supported."
            )

    @staticmethod
    def to_phifom(phase_file, work_directory: Optional[Any] = None) -> str:
        """
        Convert phase data to PHIFOM format (PHI, FOM).

        HL -> PHIFOM is computed natively (gemmi/numpy). PHIFOM input is returned as-is.

        Raises:
            CException: If validation fails or conversion not supported
        """
        PhaseDataConverter._validate_input_file(phase_file)

        output_path = phase_file._get_conversion_output_path('PHIFOM', work_directory=work_directory)
        input_path = phase_file.getFullPath()
        content_flag = int(phase_file.contentFlag)

        if content_flag == 1:  # HL -> PHIFOM
            return PhaseDataConverter._hl_to_phifom(input_path, output_path)
        elif content_flag == 2:  # already PHIFOM
            return input_path
        else:
            raise CException(
                PhaseDataConverter, 3,
                details=f"contentFlag={content_flag}, target=PHIFOM. Only HL (1) -> PHIFOM supported."
            )

    @staticmethod
    def to_fphi(map_coeffs_file, work_directory: Optional[Any] = None) -> str:
        """Return path to FPHI file (CMapCoeffsDataFile is always FPHI; no conversion)."""
        return map_coeffs_file.getFullPath()

    # ========================================================================
    # MTZ I/O (gemmi)
    # ========================================================================

    @staticmethod
    def _new_mtz_like(mtz):
        """Create an empty single-dataset MTZ sharing cell/spacegroup with `mtz`."""
        out = gemmi.Mtz(with_base=False)
        out.spacegroup = mtz.spacegroup
        out.set_cell_for_all(mtz.cell)
        out.add_dataset('HKL_base')
        out.add_column('H', 'H')
        out.add_column('K', 'H')
        out.add_column('L', 'H')
        return out

    @staticmethod
    def _hkl_columns(mtz):
        """Return H, K, L data columns as integer numpy arrays."""
        return (
            np.array(mtz.column_with_label('H'), dtype=np.int32),
            np.array(mtz.column_with_label('K'), dtype=np.int32),
            np.array(mtz.column_with_label('L'), dtype=np.int32),
        )

    @staticmethod
    def _hl_to_phifom(input_path: str, output_path: str) -> str:
        """
        Convert HL coefficients (HLA..HLD) to best phase + FOM (clipper-equivalent).

        Reads the type-'A' HL columns and the space group, computes PHI/FOM via the
        centric/acentric centroid of the HL distribution, and writes PHI (P) + FOM (W).
        """
        mtz = gemmi.read_mtz_file(input_path)

        hl_cols = {}
        for col in mtz.columns:
            if col.type == 'A':
                label = col.label.upper()
                for key in ('HLA', 'HLB', 'HLC', 'HLD'):
                    if key in label and key not in hl_cols:
                        hl_cols[key] = col
                        break
        if len(hl_cols) != 4:
            raise CException(
                PhaseDataConverter, 10,
                details=f"Found {len(hl_cols)} HL columns: {', '.join(sorted(hl_cols.keys()))}"
            )

        hla = np.array(hl_cols['HLA'], dtype=np.float64)
        hlb = np.array(hl_cols['HLB'], dtype=np.float64)
        hlc = np.array(hl_cols['HLC'], dtype=np.float64)
        hld = np.array(hl_cols['HLD'], dtype=np.float64)

        h, k, l_idx = PhaseDataConverter._hkl_columns(mtz)
        centric, allowed = PhaseDataConverter._centric_and_allowed(mtz.spacegroup, h, k, l_idx)
        phi_deg, fom = PhaseDataConverter._phifom_from_abcd(hla, hlb, hlc, hld, centric, allowed)

        out = PhaseDataConverter._new_mtz_like(mtz)
        out.add_column('PHI', 'P')
        out.add_column('FOM', 'W')
        out.set_data(np.column_stack([
            h.astype(np.float32), k.astype(np.float32), l_idx.astype(np.float32),
            phi_deg.astype(np.float32), fom.astype(np.float32)]))
        out.write_to_file(output_path)
        return output_path

    @staticmethod
    def _phifom_to_hl(input_path: str, output_path: str) -> str:
        """
        Convert PHI + FOM to HL coefficients (clipper-equivalent).

        HLA = X.cos(phi), HLB = X.sin(phi), HLC = HLD = 0, where X = atanh(FOM) for
        centric reflections and X = invsim(FOM) for acentric (FOM clamped to 0.9999).
        """
        mtz = gemmi.read_mtz_file(input_path)

        phi_col = next((c for c in mtz.columns if c.type == 'P'), None)
        fom_col = next((c for c in mtz.columns if c.type == 'W'), None)
        if phi_col is None:
            raise CException(PhaseDataConverter, 8, details=f"Input file: {input_path}")
        if fom_col is None:
            raise CException(PhaseDataConverter, 9, details=f"Input file: {input_path}")

        phi_deg = np.array(phi_col, dtype=np.float64)
        fom = np.array(fom_col, dtype=np.float64)

        h, k, l_idx = PhaseDataConverter._hkl_columns(mtz)
        centric, _ = PhaseDataConverter._centric_and_allowed(mtz.spacegroup, h, k, l_idx)
        hla, hlb, hlc, hld = PhaseDataConverter._abcd_from_phifom(phi_deg, fom, centric)

        out = PhaseDataConverter._new_mtz_like(mtz)
        out.add_column('HLA', 'A')
        out.add_column('HLB', 'A')
        out.add_column('HLC', 'A')
        out.add_column('HLD', 'A')
        out.set_data(np.column_stack([
            h.astype(np.float32), k.astype(np.float32), l_idx.astype(np.float32),
            hla.astype(np.float32), hlb.astype(np.float32),
            hlc.astype(np.float32), hld.astype(np.float32)]))
        out.write_to_file(output_path)
        return output_path

    # ========================================================================
    # Numerical core (clipper port)
    # ========================================================================

    @staticmethod
    def _centric_and_allowed(spacegroup, h, k, l):
        """
        Per-reflection centric flag and allowed centric phase (radians).

        Port of clipper HKL_class (clipper/core/coords.cpp): a reflection is centric
        if some symmetry operation maps hkl -> -hkl; the allowed phase is then
        snapped to a 15-degree grid via intr(mod(-0.5*shift, pi) / (pi/12)) * (pi/12),
        where shift = 2*pi*(h . t) for that operation.
        """
        gops = spacegroup.operations()
        ops = list(gops)
        den = float(gemmi.Op.DEN)
        pi = np.pi
        twelfth = pi / 12.0

        n = len(h)
        centric = np.zeros(n, dtype=bool)
        allowed = np.zeros(n, dtype=np.float64)

        for i in range(n):
            hkl = (int(h[i]), int(k[i]), int(l[i]))
            neg = (-hkl[0], -hkl[1], -hkl[2])
            for op in ops:
                rh = op.apply_to_hkl(hkl)
                if (rh[0], rh[1], rh[2]) == neg:
                    t = op.tran
                    hp = (hkl[0] * t[0] + hkl[1] * t[1] + hkl[2] * t[2]) / den
                    shift = 2.0 * pi * hp
                    a_idx = int(np.floor(((-0.5 * shift) % pi) / twelfth + 0.5))
                    allowed[i] = a_idx * twelfth
                    centric[i] = True
                    break
        return centric, allowed

    @staticmethod
    def _phifom_from_abcd(a, b, c, d, centric, allowed):
        """
        PHI (degrees) and FOM from HL coefficients (clipper Compute_phifom_from_abcd).

        Acentric: integrate exp(a cos + b sin + c cos2 + d sin2) over the phase circle
        in 5-degree steps (72 points). Centric: evaluate the two allowed phases only,
        using just the a, b terms. q is truncated to [-700, 700] as in clipper.
        """
        a = np.asarray(a, np.float64)
        b = np.asarray(b, np.float64)
        c = np.asarray(c, np.float64)
        d = np.asarray(d, np.float64)
        centric = np.asarray(centric, bool)

        n = a.shape[0]
        phi = np.zeros(n, np.float64)
        fom = np.zeros(n, np.float64)

        ac = ~centric
        if ac.any():
            ang = np.radians(np.arange(72) * 5.0)
            cosi, sini = np.cos(ang), np.sin(ang)
            cos2, sin2 = np.cos(2 * ang), np.sin(2 * ang)
            q = (a[ac, None] * cosi + b[ac, None] * sini
                 + c[ac, None] * cos2 + d[ac, None] * sin2)
            np.clip(q, -700.0, 700.0, out=q)
            pq = np.exp(q)
            ssum = pq.sum(axis=1)
            scos = (pq * cosi).sum(axis=1)
            ssin = (pq * sini).sum(axis=1)
            phi[ac] = np.arctan2(ssin, scos)
            fom[ac] = np.sqrt(scos * scos + ssin * ssin) / ssum

        cn = centric
        if cn.any():
            ca = np.cos(allowed[cn])
            sa = np.sin(allowed[cn])
            q = np.clip(a[cn] * ca + b[cn] * sa, -700.0, 700.0)
            pq = np.exp(q)
            diff = pq - 1.0 / pq
            ssum = pq + 1.0 / pq
            scos = diff * ca
            ssin = diff * sa
            phi[cn] = np.arctan2(ssin, scos)
            fom[cn] = np.sqrt(scos * scos + ssin * ssin) / ssum

        return np.degrees(phi), np.clip(fom, 0.0, 1.0)

    @staticmethod
    def _abcd_from_phifom(phi_deg, fom, centric):
        """
        HL coefficients from PHI/FOM (clipper Compute_abcd_from_phifom).

        x = min(FOM, 0.9999); centric -> x = atanh(x), acentric -> x = invsim(x);
        HLA = x cos(phi), HLB = x sin(phi), HLC = HLD = 0.
        """
        phi = np.radians(np.asarray(phi_deg, np.float64))
        centric = np.asarray(centric, bool)
        x = np.minimum(np.asarray(fom, np.float64), 0.9999)

        xx = np.where(centric, np.arctanh(np.minimum(x, 0.9999)),
                      PhaseDataConverter._invsim(x))
        hla = xx * np.cos(phi)
        hlb = xx * np.sin(phi)
        zero = np.zeros_like(hla)
        return hla, hlb, zero, zero

    @staticmethod
    def _invsim(x):
        """
        Inverse Sim function: invert FOM = I1(X)/I0(X) for X (acentric).

        Direct port of clipper Util::invsim (clipper/core/clipper_util.cpp), a
        closed-form solution of the relevant cubic. Vectorised over numpy arrays.
        """
        x = np.asarray(x, np.float64)
        x0 = np.abs(x)
        a0 = -7.107935 * x0
        a1 = 3.553967 - 3.524142 * x0
        a2 = 1.639294 - 2.228716 * x0
        a3 = 1.0 - x0
        w = a2 / (3.0 * a3)
        p = a1 / (3.0 * a3) - w * w
        q = -w * w * w + 0.5 * (a1 * w - a0) / a3
        dd = np.sqrt(q * q + p * p * p)
        q1 = q + dd
        q2 = q - dd
        r1 = np.abs(q1) ** (1.0 / 3.0)
        r2 = np.abs(q2) ** (1.0 / 3.0)
        val = (np.where(q1 > 0.0, r1, -r1) + np.where(q2 > 0.0, r2, -r2) - w)
        return np.where(x >= 0.0, val, -val)
