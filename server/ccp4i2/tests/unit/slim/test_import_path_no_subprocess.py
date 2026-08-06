"""
Guard: the file-import / content-detection path must stay gemmi-native and must
never spawn a subprocess.

This is the runtime sibling of test_no_ccp4_native_imports.py. That test catches
a CCP4 coupling introduced by an *import*; this one catches a coupling introduced
by a *subprocess shell-out* to a CCP4 binary, which an import guard cannot see.

The concrete risk this locks down: obs-data import detects its content flag
(IPAIR / IMEAN / FMEAN / FREER / HL ...) by reading MTZ columns with gemmi
(CObsDataFile._introspect_content_flag -> gemmi.read_mtz_file). The *conversion*
between those representations (e.g. French-Wilson I -> F) needs the `servalcat fw`
binary and is deliberately deferred to job run time inside CPluginScript.makeHklin
(see core/conversions/servalcat_converter.py). If anyone wires a binary conversion
onto the detection path -- so that merely importing/digesting a file shells out --
the slim CCP4-free control plane can no longer serve file import, and this test
trips.

Mechanism: block every subprocess entry point, then run the exact detection call
the import path performs (setContentFlag) on real demo MTZ files. gemmi reads MTZ
in-process, so a green run proves detection never left the interpreter.

Pure Python + gemmi; no CCP4 binaries required.
"""
import os
import subprocess
from pathlib import Path

import pytest

from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4XtalData import CObsDataFile, CFreeRDataFile, CPhsDataFile


GAMMA = Path(CCP4Utils.getCCP4I2Dir()) / "demo_data" / "gamma"

# (class, filename, expected contentFlag). The IPAIR case is the load-bearing one:
# detecting intensities is exactly where a careless change might "helpfully" run a
# French-Wilson conversion (servalcat) to also produce F.
DETECTION_CASES = [
    (CObsDataFile, "merged_intensities_native.mtz", 1),   # IPAIR
    (CFreeRDataFile, "freeR.mtz", 1),                     # FREER
    (CPhsDataFile, "initial_phases.mtz", 1),              # HL
]


@pytest.fixture
def no_subprocess(monkeypatch):
    """Make any subprocess/os-level process spawn an immediate failure."""
    def _boom(*args, **kwargs):
        target = args[0] if args else kwargs.get("args") or kwargs.get("command")
        raise AssertionError(
            "the file-import / content-detection path spawned a subprocess "
            f"({target!r}). This breaks the slim CCP4-free invariant: detection "
            "must be gemmi-native. Defer any binary conversion (e.g. servalcat fw) "
            "to CPluginScript.makeHklin at job run time, not the import path."
        )

    for name in ("run", "Popen", "call", "check_call", "check_output", "getoutput", "getstatusoutput"):
        if hasattr(subprocess, name):
            monkeypatch.setattr(subprocess, name, _boom)
    for name in ("system", "popen", "spawnv", "spawnve", "spawnvp", "spawnvpe"):
        if hasattr(os, name):
            monkeypatch.setattr(os, name, _boom)
    yield


def _content_flag(file_obj):
    cf = file_obj.contentFlag
    return cf.value if hasattr(cf, "value") else cf


@pytest.mark.parametrize("klass,filename,expected", DETECTION_CASES)
def test_content_detection_never_shells_out(no_subprocess, klass, filename, expected):
    """setContentFlag() (the import path's detection step) is gemmi-native: it
    detects the right flag AND never spawns a subprocess."""
    mtz = GAMMA / filename
    assert mtz.exists(), f"demo data missing: {mtz}"

    file_obj = klass(file_path=str(mtz))
    file_obj.setContentFlag()  # would raise via no_subprocess if it shelled out

    assert _content_flag(file_obj) == expected, (
        f"{klass.__name__} detected the wrong contentFlag for {filename}"
    )


def test_subprocess_block_is_actually_armed(no_subprocess):
    """Sanity check: the guard fixture really does trip on a subprocess call, so a
    green test above means 'no spawn' rather than 'fixture silently disabled'."""
    with pytest.raises(AssertionError):
        subprocess.run(["true"])
    with pytest.raises(AssertionError):
        os.system("true")
