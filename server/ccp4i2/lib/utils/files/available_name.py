import pathlib

# Extensions made of more than one dot-separated part. `pathlib.Path.suffix`
# sees only the last part, so "gamma.asu.xml" would otherwise be split into
# the stem "gamma.asu" and the suffix ".xml" -- and slugifying that stem
# destroys the very extension the datatype is identified by, leaving
# "gammaasu.xml". Longest first, so ".asu.xml" wins over a bare ".xml".
COMPOUND_SUFFIXES = (
    ".asu.xml",
    ".scene.xml",
    ".imosflm.xml",
    ".phaser_sol.pkl",
    ".phaser_rlist.pkl",
)


def split_compound_suffix(file_name) -> tuple:
    """Split a file name into (stem, suffix), honouring compound extensions.

    Behaves like ``(Path.stem, Path.suffix)`` except that a name ending in one
    of ``COMPOUND_SUFFIXES`` keeps the whole compound extension as the suffix.

    >>> split_compound_suffix("gamma.asu.xml")
    ('gamma', '.asu.xml')
    >>> split_compound_suffix("rnase.pdb")
    ('rnase', '.pdb')
    """
    name = pathlib.PurePath(str(file_name)).name
    lowered = name.lower()
    for suffix in sorted(COMPOUND_SUFFIXES, key=len, reverse=True):
        if lowered.endswith(suffix) and len(name) > len(suffix):
            return name[: -len(suffix)], name[-len(suffix):]
    path = pathlib.PurePath(name)
    return path.stem, path.suffix


def available_file_name_based_on(file_path: pathlib.Path):
    destination_dir = file_path.parent
    base_stem, base_suffix = split_compound_suffix(file_path)
    # Build by concatenation rather than with_suffix(): with_suffix() would
    # treat the first dot of a compound extension as the start of the suffix
    # and truncate it.
    dest = (destination_dir / f"{base_stem}{base_suffix}").absolute()
    attempt_number = 0
    while dest.exists():
        attempt_number += 1
        dest = (
            destination_dir / f"{base_stem}_{attempt_number}{base_suffix}"
        ).absolute()
    # Check that the filename does lie in the appropriate direction
    return dest
