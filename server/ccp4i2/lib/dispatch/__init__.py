"""Run targets: where a job, or one heavy program, runs.

CCP4i2 ships one target, ``local``, and the hook a deployment uses to add
its own. Nothing here knows any platform (decision 18 of
docs/pandda-campaign-design.md): a deployment registers a target by naming
its class in settings, and CCP4i2 resolves it by name.

    # ccp4i2.config.settings -- what CCP4i2 ships
    CCP4I2_RUN_TARGETS = {"local": "ccp4i2.lib.dispatch.local.LocalTarget"}
    CCP4I2_JOB_TARGET = "local"

    # a deployment's settings overlay -- adds a target CCP4i2 never heard of
    from ccp4i2.config.settings import *
    CCP4I2_RUN_TARGETS = {**CCP4I2_RUN_TARGETS,
                          "azure": "azure_extensions.dispatch.ServiceBusTarget"}
    CCP4I2_JOB_TARGET = "azure"

Two interfaces, in :mod:`.base`: a *job* target moves the whole CCP4i2 job
(``run_job``); a *program* target moves one external program while the job
stays in ``RUNNING_REMOTELY`` (``submit`` / ``poll`` / ``cancel`` / ``logs``).
A class may implement either or both; :func:`available_targets` reports which.

The desktop invariant: with no settings at all, ``job_target_name()`` is
``local`` and ``get_target("local")`` is the shipped subprocess runner. A
missing setting is never an error.
"""
import importlib
import os
from typing import Dict, List

DEFAULT_RUN_TARGETS: Dict[str, str] = {
    "local": "ccp4i2.lib.dispatch.local.LocalTarget",
}
DEFAULT_JOB_TARGET = "local"

_PROGRAM_METHODS = ("submit", "poll", "cancel", "logs")


class UnknownRunTarget(LookupError):
    """A target name that no setting registers."""


class RunTargetError(RuntimeError):
    """A registered target whose class cannot be imported or constructed."""


def _settings():
    """Django settings when configured, else None (registry still works)."""
    try:
        from django.conf import settings
        if settings.configured:
            return settings
    except Exception:  # noqa: BLE001 -- settings unavailable is a normal state
        pass
    return None


def run_target_paths() -> Dict[str, str]:
    """Name -> dotted class path, as the deployment registered them."""
    settings = _settings()
    paths = getattr(settings, "CCP4I2_RUN_TARGETS", None) if settings else None
    if not paths:
        return dict(DEFAULT_RUN_TARGETS)
    return {str(name).lower(): str(path) for name, path in dict(paths).items()}


def job_target_name() -> str:
    """The target that runs a whole job here: settings, then env, then local."""
    settings = _settings()
    name = getattr(settings, "CCP4I2_JOB_TARGET", None) if settings else None
    if not name:
        name = os.environ.get("CCP4I2_JOB_TARGET") or DEFAULT_JOB_TARGET
    return str(name).lower()


def _load(name: str, path: str):
    module_name, _, attr = path.rpartition(".")
    if not module_name or not attr:
        raise RunTargetError(
            f"run target '{name}' is registered as '{path}', which is not a "
            "dotted 'package.module.Class' path")
    try:
        module = importlib.import_module(module_name)
        cls = getattr(module, attr)
    except Exception as err:  # noqa: BLE001 -- report every cause the same way
        raise RunTargetError(
            f"run target '{name}' could not be loaded from '{path}': "
            f"{type(err).__name__}: {err}") from err
    try:
        return cls()
    except Exception as err:  # noqa: BLE001
        raise RunTargetError(
            f"run target '{name}' ({path}) could not be constructed: "
            f"{type(err).__name__}: {err}") from err


def get_target(name: str):
    """The registered target instance for ``name``.

    Raises :class:`UnknownRunTarget` naming what *is* registered, so a typo in
    a deployment's settings is a one-line diagnosis, not a stack trace.
    """
    key = str(name).lower()
    paths = run_target_paths()
    if key not in paths:
        raise UnknownRunTarget(
            f"no run target named '{name}'; registered: "
            + ", ".join(sorted(paths)))
    return _load(key, paths[key])


def runs_jobs(target) -> bool:
    return callable(getattr(target, "run_job", None))


def runs_programs(target) -> bool:
    return all(callable(getattr(target, m, None)) for m in _PROGRAM_METHODS)


def available_targets() -> List[dict]:
    """What this deployment can do, for the API and the UI to report.

    One entry per registered name; a target that fails to load is reported
    with its error rather than hidden, because a mis-registered target is
    exactly what an operator needs to see.
    """
    out = []
    for name, path in sorted(run_target_paths().items()):
        try:
            target = _load(name, path)
        except RunTargetError as err:
            out.append({"name": name, "path": path, "error": str(err),
                        "runs_jobs": False, "runs_programs": False})
            continue
        out.append({"name": name, "path": path,
                    "runs_jobs": runs_jobs(target),
                    "runs_programs": runs_programs(target)})
    return out
