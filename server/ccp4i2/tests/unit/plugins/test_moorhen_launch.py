"""i2run moorhen asks the desktop app to open the session window through
CCP4I2_DESKTOP_LAUNCH; without it, it only logs the route."""

import json
import subprocess

from ccp4i2.wrappers.moorhen.script import moorhen as moorhen_module


def test_launches_the_app_with_open_route(monkeypatch):
    calls = []

    def fake_popen(args, **kwargs):
        calls.append((args, kwargs))
        return object()

    monkeypatch.setattr(subprocess, "Popen", fake_popen)
    environ = {moorhen_module.DESKTOP_LAUNCH_ENV: json.dumps(["/Applications/x.app/Contents/MacOS/x"])}
    assert moorhen_module.launch_session_window("/ccp4i2/moorhen-page/session/7", environ) is not None
    assert calls[0][0] == ["/Applications/x.app/Contents/MacOS/x", "--open-route",
                           "/ccp4i2/moorhen-page/session/7"]
    assert calls[0][1]["start_new_session"] is True


def test_no_launch_command_means_no_spawn(monkeypatch):
    def never(*a, **k):
        raise AssertionError("spawned")

    monkeypatch.setattr(subprocess, "Popen", never)
    assert moorhen_module.launch_session_window("/ccp4i2/x", {}) is None
    assert moorhen_module.launch_session_window("/ccp4i2/x", {moorhen_module.DESKTOP_LAUNCH_ENV: "not json"}) is None
    assert moorhen_module.launch_session_window("/ccp4i2/x", {moorhen_module.DESKTOP_LAUNCH_ENV: "[]"}) is None


def test_spawn_failure_is_logged_not_raised(monkeypatch):
    def failing(*a, **k):
        raise OSError("no such executable")

    monkeypatch.setattr(subprocess, "Popen", failing)
    environ = {moorhen_module.DESKTOP_LAUNCH_ENV: json.dumps(["/missing"])}
    assert moorhen_module.launch_session_window("/ccp4i2/x", environ) is None
