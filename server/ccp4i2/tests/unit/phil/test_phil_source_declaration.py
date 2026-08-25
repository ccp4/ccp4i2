"""Declaring where a task's PHIL comes from.

Every PHIL wrapper used to implement get_master_phil() with the same few lines
of import-and-return, in one of three shapes. They now declare the shape and
the base class resolves it, so the import and its failure mode are written
once.

Overriding get_master_phil() directly is still supported, for tools none of
the three shapes fit.
"""

import pytest

from ccp4i2.core.PhilPluginScript import PhilPluginScript


class _Bare(PhilPluginScript):
    TASKNAME = "bare"


class _Scope(PhilPluginScript):
    TASKNAME = "scope"
    PHIL_SCOPE = "libtbx.phil:parse"


class _Both(PhilPluginScript):
    TASKNAME = "both"
    PHIL_SCOPE = "libtbx.phil:parse"
    PHIL_PROGRAM = "somewhere:Program"


class _Missing(PhilPluginScript):
    TASKNAME = "missing"
    PHIL_SCOPE = "a_tool_that_is_not_installed:master_phil"


def test_a_declared_scope_resolves():
    """PHIL_SCOPE names a module attribute and returns exactly that object."""
    from libtbx.phil import parse

    assert PhilPluginScript.get_master_phil(_Scope.__new__(_Scope)) is parse


def test_declaring_nothing_is_an_error():
    with pytest.raises(NotImplementedError, match="PHIL_SCOPE"):
        PhilPluginScript.get_master_phil(_Bare.__new__(_Bare))


def test_declaring_two_sources_is_an_error():
    """Two declarations have no defined precedence, so say so rather than pick."""
    with pytest.raises(ValueError, match="more than one PHIL source"):
        PhilPluginScript.get_master_phil(_Both.__new__(_Both))


def test_an_uninstalled_tool_yields_none_rather_than_raising():
    """The task must still open so it can say what it needs.

    _merge_phil_parameters() treats None as "no parameters to merge"; a raised
    ImportError here would stop the task loading at all.
    """
    assert PhilPluginScript.get_master_phil(_Missing.__new__(_Missing)) is None


def test_params_file_declaration_reads_from_the_package(tmp_path, monkeypatch):
    """PHIL_PARAMS_FILE locates a file shipped inside an installed package."""
    import sys
    import types

    package = types.ModuleType("fake_tool")
    package.__file__ = str(tmp_path / "__init__.py")
    (tmp_path / "params").mkdir()
    (tmp_path / "params" / "master.params").write_text("threshold = 3\n  .type = int\n")
    monkeypatch.setitem(sys.modules, "fake_tool", package)

    class _File(PhilPluginScript):
        TASKNAME = "file"
        PHIL_PARAMS_FILE = "fake_tool:params/master.params"

    scope = PhilPluginScript.get_master_phil(_File.__new__(_File))
    assert [o.path for o in scope.all_definitions()] == ["threshold"]
