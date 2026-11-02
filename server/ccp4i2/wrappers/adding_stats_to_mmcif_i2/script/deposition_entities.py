"""The entity bookkeeping moved to ``ccp4i2.lib.deposition.entities`` when the
``prepare_deposition`` task replaced this one; this module keeps the old import
path working for existing jobs and tests."""
from ccp4i2.lib.deposition.entities import (  # noqa: F401
    entity_summary,
    rebuild_entities,
    rewrite_entity_poly,
)
