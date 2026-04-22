from typing import Iterable, Optional

from compounds.registry.models import Compound, Target


def infer_target_from_compound_ids(compound_ids: Iterable) -> Optional[Target]:
    """Return the single Target shared by all given Compounds, or None.

    Returns None when the compound ids are empty, when no compound has a
    target, or when the compounds disagree on target. Only returns a Target
    when every non-null target_id across the set is the same.
    """
    ids = [cid for cid in compound_ids if cid]
    if not ids:
        return None

    target_ids = set(
        Compound.objects.filter(id__in=ids)
        .exclude(target_id__isnull=True)
        .values_list('target_id', flat=True)
        .distinct()
    )
    if len(target_ids) != 1:
        return None

    (sole_id,) = target_ids
    return Target.objects.filter(id=sole_id).first()
