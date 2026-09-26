"""Coot's headless API, across CCP4 vintages.

CCP4 renamed the molecules-container class: ``molecules_container_py`` exists
in ccp4-20260702 alongside ``molecules_container_t``, and was dropped by
ccp4-20260904, which keeps only ``molecules_container_t``. Four call sites
constructed it by the old name, so moving the container images to the
September CCP4 broke SubstituteLigand's ligand fitting and the three coot_*
wrappers at once -- each with the same AttributeError, raised minutes into a
run rather than at import.

Resolving the name in one place means the next rename is one edit, and means a
build that has neither name says so clearly instead of failing as a missing
attribute inside a pipeline.

The import stays inside the function deliberately: coot_headless_api is an
external Coot API present only in the execution (worker) environment, and
importing it at module scope would break the guarantee that every task plugin
imports with no CCP4 present.
"""

# Newest name first: a build carrying both should use the current one.
_CONTAINER_NAMES = ("molecules_container_t", "molecules_container_py")


def molecules_container(*args, **kwargs):
    """Construct a Coot molecules container, whatever this CCP4 calls it.

    Takes the same arguments as the underlying class (callers pass a single
    boolean for verbosity).
    """
    import coot_headless_api

    for name in _CONTAINER_NAMES:
        factory = getattr(coot_headless_api, name, None)
        if factory is not None:
            return factory(*args, **kwargs)

    raise AttributeError(
        "coot_headless_api exposes none of {}; this CCP4 build cannot be used "
        "for Coot ligand fitting or map morphing. Present attributes: {}".format(
            ", ".join(_CONTAINER_NAMES),
            ", ".join(n for n in dir(coot_headless_api) if not n.startswith("_")),
        )
    )
