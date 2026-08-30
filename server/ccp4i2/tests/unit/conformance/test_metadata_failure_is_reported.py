"""A model that fails to build must say so.

`_apply_metadata_attributes` used to end with `except Exception: pass`,
commented "skip silently to avoid breaking existing code". The failure mode it
produced is not a crash but *silence*: a class whose declared attributes fail to
build yields an object with no parameters, which looks exactly like a class that
declares none. Every subsequent check then passes vacuously.

That is not hypothetical. A one-line NameError introduced while collapsing the
child structures presented as fifty unrelated assertion failures elsewhere, all
of the form "'NoneType' object has no attribute 'set'", and the cause was only
visible by calling the builder directly.

Measured across all 171 registered tasks, that except caught nothing, so the
tolerance protected nothing.

It still does not raise: objects are constructed during introspection --- the
task registry, the container snapshot, the task chooser --- and one malformed
class must not take the registry down. It is logged and recorded instead, and
`CPluginScript.validity()` surfaces it.

Pure Python -- no CCP4 binaries needed.
"""
import pytest

django = pytest.importorskip("django")


@pytest.fixture
def failing_builder(monkeypatch):
    """Make the metadata builder raise, as a declaration bug would."""
    import ccp4i2.core.base_object.class_metadata as cm

    def boom(instance):
        raise NameError("simulated builder bug")

    monkeypatch.setattr(cm, "apply_metadata_to_instance", boom)


def test_nothing_fails_to_build_today():
    """The baseline the change rests on: the tolerance protects nothing."""
    from ccp4i2.core.tasks import TASKS, get_plugin_class

    built, failures = 0, []
    for name in TASKS:
        try:
            plugin = get_plugin_class(name)(parent=None, name=name)
        except Exception:
            continue
        built += 1
        seen = set()

        def walk(obj, depth=0):
            if id(obj) in seen or depth > 8:
                return
            seen.add(id(obj))
            if getattr(obj, "_metadata_failure", None) is not None:
                failures.append(f"{name}: {type(obj).__name__}")
            for child in obj.children():
                walk(child, depth + 1)

        walk(plugin.container)

    assert built > 100, f"only {built} tasks built --- proves little"
    assert not failures, f"objects whose model failed to build: {failures[:6]}"


def test_a_failure_is_recorded_on_the_object(failing_builder):
    from ccp4i2.core.base_object.fundamental_types import CInt

    obj = CInt()
    assert obj._metadata_failure is not None
    assert "NameError" in obj._metadata_failure


def test_construction_still_succeeds(failing_builder):
    """It must not raise: the registry is built by constructing every task."""
    from ccp4i2.core.tasks import get_plugin_class

    plugin = get_plugin_class("freerflag")(parent=None, name="f")
    assert plugin is not None


def test_validity_reports_it(failing_builder):
    """The elevation. Without this the task merely looks empty."""
    from ccp4i2.core.tasks import get_plugin_class

    plugin = get_plugin_class("freerflag")(parent=None, name="f")
    reported = [str(e) for e in plugin.validity().entries()
                if "could not be built" in str(e)]
    assert reported, "a model that failed to build was not reported"
    assert any("NameError" in r for r in reported), \
        "the report does not say what went wrong"


def test_a_missing_metadata_system_is_not_an_error(monkeypatch):
    """ImportError stays distinct: nothing declares attributes, so nothing is
    missing. That is the one case the old code got right."""
    import builtins
    real_import = builtins.__import__

    def no_metadata(name, *args, **kwargs):
        if name.endswith("class_metadata"):
            raise ImportError("simulated")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_metadata)
    from ccp4i2.core.base_object.fundamental_types import CInt

    obj = CInt()
    assert obj._metadata_failure is None
