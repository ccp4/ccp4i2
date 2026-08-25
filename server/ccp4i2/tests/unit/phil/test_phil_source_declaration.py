"""Declaring where a task's PHIL comes from.

Every PHIL wrapper used to implement get_master_phil() with the same few lines
of import-and-return, in one of three shapes. They now declare the shape and
the base class resolves it, so the import and its failure mode are written
once.

Overriding get_master_phil() directly is still supported, for tools none of
the three shapes fit.

Resolution is deliberately generic — import a module, take an attribute — so
most of it can be tested without libtbx, which is CCP4/cctbx-only and absent
from the CCP4-free test run. Only the PHIL_PARAMS_FILE case has to parse PHIL,
and it is skipped when iotbx is missing.
"""

import json

import pytest

from ccp4i2.core.PhilPluginScript import PhilPluginScript


class _Bare(PhilPluginScript):
    TASKNAME = "bare"


class _Scope(PhilPluginScript):
    TASKNAME = "scope"
    # Any importable attribute exercises the resolver; using a stdlib one
    # keeps this test honest about what it covers, and CCP4-free.
    PHIL_SCOPE = "json:dumps"


class _Both(PhilPluginScript):
    TASKNAME = "both"
    PHIL_SCOPE = "json:dumps"
    PHIL_PROGRAM = "somewhere:Program"


class _Missing(PhilPluginScript):
    TASKNAME = "missing"
    PHIL_SCOPE = "a_tool_that_is_not_installed:master_phil"


def resolve(cls):
    """get_master_phil() on an uninstantiated subclass — no plugin needed."""
    return PhilPluginScript.get_master_phil(cls.__new__(cls))


def test_a_declared_scope_resolves():
    """PHIL_SCOPE names a module attribute and returns exactly that object."""
    assert resolve(_Scope) is json.dumps


def test_declaring_nothing_is_an_error():
    with pytest.raises(NotImplementedError, match="PHIL_SCOPE"):
        resolve(_Bare)


def test_declaring_two_sources_is_an_error():
    """Two declarations have no defined precedence, so say so rather than pick."""
    with pytest.raises(ValueError, match="more than one PHIL source"):
        resolve(_Both)


def test_an_uninstalled_tool_yields_none_rather_than_raising():
    """The task must still open so it can say what it needs.

    _merge_phil_parameters() treats None as "no parameters to merge"; a raised
    ImportError here would stop the task loading at all.
    """
    assert resolve(_Missing) is None


def test_params_file_declaration_reads_from_the_package(tmp_path, monkeypatch):
    """PHIL_PARAMS_FILE locates a file shipped inside an installed package."""
    pytest.importorskip("iotbx.phil", reason="needs iotbx (CCP4/cctbx)")

    import sys
    import types

    package = types.ModuleType("fake_tool")
    package.__file__ = str(tmp_path / "__init__.py")
    (tmp_path / "params").mkdir()
    (tmp_path / "params" / "master.params").write_text(
        "threshold = 3\n  .type = int\n"
    )
    monkeypatch.setitem(sys.modules, "fake_tool", package)

    class _File(PhilPluginScript):
        TASKNAME = "file"
        PHIL_PARAMS_FILE = "fake_tool:params/master.params"

    scope = resolve(_File)
    assert [o.path for o in scope.all_definitions()] == ["threshold"]
