"""Non-finite KPI values: what may be stored, and what may go on the wire.

A KPI that came out NaN or infinite is not a measurement. It is a measurement
that failed, and the honest record of that is *no row*, not a row saying "nan".
Three things then make it everyone's problem:

* JSON has no spelling for NaN or infinity, so such a value 500s whatever
  endpoint serves it — see `ccp4i2.lib.json_safety` for why, and for why
  relaxing ``STRICT_JSON`` is not the fix.
* ``JobFloatValue.value`` is ``NOT NULL``. SQLite refuses to store NaN there at
  all (IntegrityError), while PostgreSQL stores it happily. So the same bad KPI
  silently truncates a job's KPIs on the desktop and 500s an endpoint on the
  deployed instances. Infinity, by contrast, stores fine on both.
* The value outlives the job: ``export_project`` writes it to db.xml and the
  importers read it back, so one bad number can propagate between databases.

Hence the same predicate is applied at three depths, each of which would be
sufficient if the others were perfect, and none of which is:

``is_storable_kpi_value``
    The write gate. Gleaning, backfill and the importers drop a non-finite
    value instead of storing it, so it never enters a database.
``kpi_map``
    The read gate. Builds the ``{key: value}`` KPI dicts the API serves,
    omitting anything unservable. Omitting rather than nulling is deliberate:
    the client renders ``kpis.RFactor !== undefined ? ... : '-'``, so a null
    would pass that guard and display a confident ``0.000``.
``ccp4i2.api.renderers.SafeJSONRenderer``
    The backstop, one layer further out: it nulls any non-finite float in any
    response, KPI or not, so no such value can 500 the API again. It uses the
    shared scrub in `ccp4i2.lib.json_safety`, which also serves file digests.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, Iterable

from .json_safety import is_finite_number, replace_non_finite  # noqa: F401

logger = logging.getLogger(__name__)


def is_storable_kpi_value(value: Any) -> bool:
    """Is this KPI value one we are willing to keep?

    False only for floats JSON cannot carry (NaN, +inf, -inf). Everything
    else — ints, strings, ordinary floats — passes; this is a check for
    unrepresentable numbers, not a check for plausible crystallography.
    """
    return is_finite_number(value)


def drop_unstorable(values: Dict[str, Any], context: str = "") -> Dict[str, Any]:
    """Return `values` without its non-finite entries, saying what it dropped.

    `context` is put in the log line so a dropped KPI can be traced back to
    the job or file it came from.
    """
    kept = {}
    for key, value in values.items():
        if is_storable_kpi_value(value):
            kept[key] = value
        else:
            logger.warning(
                "Dropping non-finite KPI %s=%r%s: not a measurement, and not "
                "representable in JSON",
                key, value, f" ({context})" if context else "",
            )
    return kept


def kpi_map(rows: Iterable[Any], context: str = "") -> Dict[str, Any]:
    """Build the ``{key_name: value}`` KPI dict the API serves.

    `rows` are JobFloatValue or JobCharValue instances. ``JobValueKey.name`` is
    that model's primary key, so ``row.key_id`` is already the key name string
    and no join is needed to read it.

    Non-finite values are omitted rather than nulled — see the module docstring.
    They should not be in the database at all, but older rows predate the write
    gate and the importers can carry them in from legacy projects.
    """
    values = {}
    for row in rows:
        if is_storable_kpi_value(row.value):
            values[row.key_id] = row.value
        else:
            logger.warning(
                "Omitting non-finite KPI %s=%r from API payload%s; run "
                "`manage.py prune_nonfinite_kpis --apply` to clear it",
                row.key_id, row.value, f" ({context})" if context else "",
            )
    return values
