"""
Guard: the Django API/control-plane surface must not import anything from the
CCP4 main distribution.

The server is designed to run on a slim, CCP4-free interpreter (stock CPython +
pip gemmi); only job execution (wrappers/pipelines, which run out-of-process in a
ccp4-python worker/subprocess) may touch the CCP4 native stack. This test imports
the realistic API surface — the request-handling entry points and the transitive
closure of everything under ``ccp4i2.api`` and ``ccp4i2.lib`` — and asserts that
no CCP4-suite-only module gets pulled in.

It catches a regression introduced two ways:
  * on a slim interpreter, a new module-level ``import clipper`` fails with a
    native ModuleNotFoundError during the walk;
  * on ccp4-python (where the native modules ARE importable), the native module
    shows up in ``sys.modules``.

Runs in a fresh subprocess so the result is independent of whatever pytest has
already imported, and is meaningful on both interpreters.

If this test fails, something on an API code path now needs CCP4. Either make it
lazy/gemmi-native, or move it behind the job-execution (worker) boundary.
"""

import os
import subprocess
import sys

import pytest

# CCP4-suite-only modules (no pip equivalent; require the full distribution).
NATIVE = sorted({
    "clipper", "cctbx", "scitbx", "iotbx", "mmtbx", "libtbx", "cbflib_adaptbx",
    "mmdb", "mmdb2", "ccp4mg", "phaser", "phaser_ext", "hklfile", "pyrvapi",
    "ccp4srs", "dials", "dxtbx", "xia2", "simbad", "mrbump", "iris_validation",
    "acedrg", "morda", "arcimboldo",
})

_CHILD = r'''
import sys, importlib, pkgutil
import django
django.setup()

NATIVE = set(%(native)r)

def native_misses(msg):
    return {n for n in NATIVE
            if ("No module named '%%s'" %% n) in msg
            or ("No module named '%%s." %% n) in msg}

errs = set()

# Server startup chain first (asgi -> urls -> viewsets -> digest -> ...).
import ccp4i2.api.urls  # noqa: F401

import ccp4i2.api as _api
import ccp4i2.lib as _lib
for pkg in (_api, _lib):
    for mi in pkgutil.walk_packages(pkg.__path__, pkg.__name__ + "."):
        if ".tests" in mi.name or mi.name.endswith(".tests"):
            continue
        try:
            importlib.import_module(mi.name)
        except Exception as e:               # noqa: BLE001
            errs |= native_misses(str(e))    # only native misses count

present = {m.split(".")[0] for m in sys.modules} & NATIVE
leaked = sorted(present | errs)
if leaked:
    print("LEAKED:" + ",".join(leaked))
    sys.exit(1)
print("CLEAN")
sys.exit(0)
''' % {"native": NATIVE}


def test_api_surface_imports_no_ccp4_native():
    env = dict(os.environ)
    env.setdefault("DJANGO_SETTINGS_MODULE", "ccp4i2.config.settings")
    env.setdefault("CCP4I2_BACKEND", "django")
    # Ensure the server package is importable in the child.
    server_root = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "..", "..", "..")
    )
    env["PYTHONPATH"] = os.pathsep.join(
        [server_root, env.get("PYTHONPATH", "")]
    ).rstrip(os.pathsep)

    result = subprocess.run(
        [sys.executable, "-c", _CHILD],
        env=env, capture_output=True, text=True, timeout=300,
    )
    out = (result.stdout + result.stderr).strip()
    assert result.returncode == 0, (
        "The API surface imported CCP4-native module(s) — it must run on a "
        "slim CCP4-free interpreter. Make the import lazy/gemmi-native or move "
        f"it behind the worker boundary.\n{out}"
    )
    assert "CLEAN" in result.stdout, f"Unexpected guard output:\n{out}"
