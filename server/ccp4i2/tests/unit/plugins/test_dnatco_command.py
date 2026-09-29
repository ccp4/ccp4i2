"""How the dnatco wrapper decides what to launch DNATCO with.

The CCP4 launcher (``dnatco.sh``, ``dnatco.bat`` on Windows) does one thing:
``node $CCP4/dnatco/bin/dnatco.js``. Which of the two equivalent routes is
taken matters on Windows, where a job's ``subprocess`` call cannot execute a
batch file at all.
"""
import pytest

from ccp4i2.wrappers.dnatco.script import dnatco


@pytest.fixture
def routes(monkeypatch, tmp_path):
    """Control both routes: which programs resolve, and whether dnatco.js exists."""
    state = {"resolvable": set(), "js": None}

    def fake_resolve(name):
        return f"/fake/bin/{name}" if name in state["resolvable"] else None

    monkeypatch.setattr(dnatco, "resolve_program", fake_resolve)
    monkeypatch.setattr(dnatco, "dnatco_js_path", lambda: state["js"])
    return state


def test_launcher_name_is_per_platform(monkeypatch):
    monkeypatch.setattr(dnatco.os, "name", "posix")
    assert dnatco.dnatco_launcher_name() == "dnatco.sh"
    monkeypatch.setattr(dnatco.os, "name", "nt")
    assert dnatco.dnatco_launcher_name() == "dnatco.bat"


def test_posix_prefers_the_launcher_script(monkeypatch, routes, tmp_path):
    monkeypatch.setattr(dnatco.os, "name", "posix")
    routes["resolvable"] = {"dnatco.sh", "node"}
    routes["js"] = tmp_path / "dnatco.js"
    # Both routes available: the relocatable launcher wins.
    assert dnatco.find_dnatco_command() == "dnatco.sh"


def test_posix_falls_back_to_node(monkeypatch, routes, tmp_path):
    monkeypatch.setattr(dnatco.os, "name", "posix")
    routes["resolvable"] = {"node"}
    routes["js"] = tmp_path / "dnatco.js"
    assert dnatco.find_dnatco_command() == "node"


def test_windows_never_picks_the_batch_file_when_node_can_run_it(monkeypatch, routes, tmp_path):
    # CreateProcess cannot execute a .bat, and a job is spawned with no shell,
    # so on Windows the node route is the one that actually works.
    monkeypatch.setattr(dnatco.os, "name", "nt")
    routes["resolvable"] = {"dnatco.bat", "node"}
    routes["js"] = tmp_path / "dnatco.js"
    assert dnatco.find_dnatco_command() == "node"


def test_missing_install_is_named_after_the_launcher(monkeypatch, routes):
    # Nothing resolves: report "dnatco.sh not found", not "node not found".
    monkeypatch.setattr(dnatco.os, "name", "posix")
    assert dnatco.find_dnatco_command() == "dnatco.sh"
    monkeypatch.setattr(dnatco.os, "name", "nt")
    assert dnatco.find_dnatco_command() == "dnatco.bat"


def test_js_path_needs_ccp4_and_the_file(monkeypatch, tmp_path):
    monkeypatch.delenv("CCP4", raising=False)
    assert dnatco.dnatco_js_path() is None
    monkeypatch.setenv("CCP4", str(tmp_path))
    assert dnatco.dnatco_js_path() is None          # $CCP4 set, DNATCO not shipped
    js = tmp_path / "dnatco" / "bin" / "dnatco.js"
    js.parent.mkdir(parents=True)
    js.write_text("// dnatco\n")
    assert dnatco.dnatco_js_path() == js
