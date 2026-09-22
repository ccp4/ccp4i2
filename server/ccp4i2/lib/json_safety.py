"""Floats that JSON cannot spell.

JSON has no NaN and no infinity. DRF renders with ``allow_nan=False``
(``STRICT_JSON``), so one such float anywhere in a payload raises ValueError
*during rendering* -- after the view has returned, where no view-level
``try/except`` can catch it -- and the endpoint answers 500. Relaxing
``STRICT_JSON`` is not the fix: it emits a bare ``NaN`` token, which
``JSON.parse`` rejects in the browser.

They arrive in real data from more than one direction: gemmi's CifToMtz writes
a NaN dataset wavelength when the structure-factor mmCIF records none, and a
crystallographic program can report a KPI that did not converge. NaN is also
truthy and compares unequal to everything, so neither ``if value:`` nor
``value != 0.0`` keeps it out.

This module holds the one definition. `ccp4i2.lib.kpi_values` applies it to
KPIs, `ccp4i2.lib.utils.files.digest` to file digests, and
`ccp4i2.api.renderers.SafeJSONRenderer` to every response as a backstop. It
imports nothing but the standard library, so it is available to the CCP4-free
request path.
"""

import math
from typing import Any


def is_finite_number(value: Any) -> bool:
    """False only for a float JSON cannot represent: NaN, +inf, -inf.

    Everything else -- ints, strings, None, ordinary floats -- is True. This
    asks whether a value can go on the wire, not whether it is sensible.
    """
    return not (isinstance(value, float) and not math.isfinite(value))


def replace_non_finite(data: Any) -> Any:
    """Copy `data` with every non-finite float replaced by None.

    Dicts and lists are rebuilt, tuples become lists (JSON has no tuple
    either way), and anything else is passed through. The input is not
    mutated.
    """
    if isinstance(data, float):
        return data if math.isfinite(data) else None
    if isinstance(data, dict):
        return {key: replace_non_finite(value) for key, value in data.items()}
    if isinstance(data, (list, tuple)):
        return [replace_non_finite(item) for item in data]
    return data
