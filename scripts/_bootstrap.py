"""Import the tree this script lives in, and prove it.

Every script under `scripts/` needs the same three lines to reach `ccp4i2`, and
getting them subtly wrong is not loud --- it is silent and confident.
`snapshot_containers.py` inserted `os.getcwd()` and trusted `import ccp4i2` to
follow. That works from `server/`; run from a worktree root it inserts a
directory with no `ccp4i2` in it, and whatever the interpreter finds next wins.

What it finds next differs per interpreter, and neither announces itself:

  - the dev venv carries an *editable* install pinned to the main checkout, so
    a script run from a worktree silently measured the main tree --- and every
    comparison then reported "byte-identical" because it compared that tree
    with itself
  - `ccp4-python` on a CCP4 build that still bundles the Qt-era i2 finds a
    namespace package at `site-packages/ccp4i2/` which claims the name outright,
    so the new backend is not importable at all

So the fix cannot be "remember to set PYTHONPATH". Locate the tree from
`__file__`, then *check* that the import actually came from there, and refuse to
run if it did not. A tool that cannot say which tree it measured is not
evidence.

Usage:

    from _bootstrap import setup           # scripts/ is on sys.path already
    setup()                                # or setup(django=False)
"""
import os
import sys

#: <repo>/server --- derived from this file, not from the working directory.
SERVER = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "server"
)


def setup(django: bool = True,
          settings: str = "ccp4i2.config.test_settings") -> str:
    """Put this checkout's `server/` first on `sys.path` and verify it won.

    Returns the server directory. Raises SystemExit --- loudly, before any work
    --- if `ccp4i2` resolves anywhere else.
    """
    if SERVER not in sys.path:
        sys.path.insert(0, SERVER)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", settings)

    import ccp4i2
    verify(ccp4i2)

    if django:
        import django as _django
        _django.setup()
    return SERVER


def verify(module) -> None:
    """Refuse to continue if `module` did not come from this checkout."""
    got = getattr(module, "__file__", None)
    if got is None:                       # a namespace package claimed the name
        raise SystemExit(
            f"refusing to run: `import {module.__name__}` produced a namespace "
            f"package with no __file__, searching {list(getattr(module, '__path__', []))}.\n"
            f"Expected the tree at {SERVER}. On a CCP4 build that still bundles "
            "the Qt-era ccp4i2, its site-packages directory claims the name."
        )
    got = os.path.abspath(got)
    want = os.path.join(SERVER, module.__name__) + os.sep
    if not got.startswith(want):
        raise SystemExit(
            f"refusing to run: this script lives in {SERVER}\n"
            f"but `import {module.__name__}` resolved to {got}\n"
            "Any result would describe the wrong tree."
        )
