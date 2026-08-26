"""Gate a test on whether the program it drives can actually be found.

"Installed" has to mean the same thing to the test suite as it does to a job,
or the suite is testing a different application from the one that ships.
CCP4i2 resolves an executable through
:func:`ccp4i2.config.program_discovery.discover_program`, which tries, in
order: a ``{PROG}_EXECUTABLE`` preference, a ``{SUITE}DIR`` preference
(``SHELXDIR`` and friends), the ``exePaths`` list, and finally ``PATH``.

The suite used to probe with ``shutil.which`` --- the fourth of those four, and
the only one that ignores the user's preferences. A developer with ``SHELXDIR``
set correctly had SHELX available to every job while the tests announced
"SHELXE not installed" and skipped. Others did not probe at all and simply
failed, and one was marked ``skip`` unconditionally, which never re-runs however
the machine is configured.

Use ``@requires_program("shelxc")``. A probe self-heals: the test starts running
the moment the program becomes resolvable.
"""
from pytest import mark

from ccp4i2.config.program_discovery import discover_program

#: Programs that gate at least one test here, reported in the run header so a
#: baseline records *how* each was found --- or that it was not. A run where a
#: program resolves through a user preference rather than the CCP4 installation
#: is a different run from one where it resolves through $CBIN, and the
#: difference belongs in the record rather than in someone's memory.
GATED_PROGRAMS = ("shelxc", "shelxd", "shelxe", "phaser", "xds", "xds_par")


def missing_programs(*names):
    """Which of *names* cannot be resolved the way a job would resolve them."""
    return [name for name in names if discover_program(name)["path"] is None]


def requires_program(*names):
    """Skip unless every one of *names* resolves.

    Example::

        @requires_program("shelxc", "shelxd")
        def test_substrdet():
            ...
    """
    missing = missing_programs(*names)
    if not missing:
        return mark.skipif(False, reason="")
    return mark.skipif(
        True,
        reason=("not installed: " + ", ".join(missing) +
                " (looked in the *_EXECUTABLE and *DIR preferences, exePaths "
                "and PATH)"),
    )


def resolution_report(names=GATED_PROGRAMS):
    """One line per program: where it was found, or that it was not.

    Emitted into the pytest header, so a stored baseline says which SHELX ran.
    """
    lines = []
    for name in names:
        found = discover_program(name)
        if found["path"]:
            lines.append(f"  {name:10s} {found['source']:20s} {found['path']}")
        else:
            lines.append(f"  {name:10s} missing")
    return lines
