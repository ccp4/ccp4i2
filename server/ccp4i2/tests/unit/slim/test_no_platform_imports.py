"""
CCP4i2 imports no cloud platform SDK. Anywhere.

The rule (docs/pandda-campaign-design.md, decision 18): CCP4i2 builds generic
capabilities with hooks; a deployment registers what it has, in its own
package. The moment an ``import azure...`` lands in this tree, the desktop
build has acquired a cloud dependency it cannot use and the next platform
must edit CCP4i2 to exist. So CI reads every module and refuses one.

The check is on SDK packages, not on the name of a deployment's Django app:
``ccp4i2/api/urls.py`` may still mount ``azure_extensions`` URLs when that app
is installed, which is a hook, not an SDK.
"""
import ast
import pathlib

import pytest

PACKAGE = pathlib.Path(__file__).resolve().parents[3]   # .../ccp4i2
FORBIDDEN = ("azure", "boto3", "botocore", "google.cloud", "kubernetes")


def _imports(path):
    tree = ast.parse(path.read_text(encoding="utf-8", errors="replace"), str(path))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                yield alias.name
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            yield node.module


def _forbidden(name):
    return any(name == f or name.startswith(f + ".") for f in FORBIDDEN)


def test_no_module_imports_a_platform_sdk():
    offenders = []
    for path in PACKAGE.rglob("*.py"):
        if "tests" in path.relative_to(PACKAGE).parts:
            continue
        for name in _imports(path):
            if _forbidden(name):
                offenders.append(f"{path.relative_to(PACKAGE)}: import {name}")
    assert not offenders, (
        "platform SDK imported inside ccp4i2 (register a run target in the "
        "deployment's own package instead):\n" + "\n".join(offenders))


def test_the_guard_itself_sees_imports():
    assert _forbidden("azure.servicebus") and _forbidden("boto3")
    assert not _forbidden("azure_extensions.urls")   # a deployment app hook, not an SDK
    assert not _forbidden("ccp4i2.lib.dispatch")
