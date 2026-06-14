"""
Behavioural guard: no task's validity()/runTimeValidity() may require the CCP4
install.

Complements test_no_ccp4_native_imports.py. That guard fences Class A (CCP4-native
*imports*) on the API surface; this one fences Class B (env-driven discovery of
CCP4 reference data) and Class C (CCP4 binary execution) on the validation path,
which the import guard cannot see. The real invariant is "nothing on the
request/validation path requires the CCP4 install" — and validity()/
runTimeValidity() run on the API via the /validation/ and /run_time_validation/
endpoints.

The work is done by _ccp4_free_validation_probe.py, run here in a subprocess whose
CCP4 install is scrubbed (env vars unset + CCP4 bin removed from PATH), so the
check is faithful on a slim interpreter AND under ccp4-python.

If this fails: a task reads CCP4 reference data or runs a CCP4 binary from
validity()/runTimeValidity(). Either defer it (return SUCCEEDED/skip when the
install is absent; the worker re-runs it in process()) or port it to pure
gemmi/numpy. Cf. PR #175 (checkMonomeCoverage defer) and the crossec->gemmi port.
"""

import os
import subprocess
import sys
from pathlib import Path

SERVER_ROOT = Path(__file__).resolve().parents[4]  # .../server
PROBE = Path(__file__).parent / "_ccp4_free_validation_probe.py"

# CCP4-install variables to remove so the probe sees no install.
_CCP4_ENV_VARS = (
    "CCP4", "CLIBD", "CLIBD_MON", "CBIN", "CLIB", "CINCL",
    "CCP4_SCR", "CETC", "CCP4_MASTER", "CCP4I2_TOP",
)


def test_validation_runs_with_no_ccp4_install(tmp_path):
    env = dict(os.environ)
    env.setdefault("DJANGO_SETTINGS_MODULE", "ccp4i2.config.settings")
    env.setdefault("CCP4I2_BACKEND", "django")
    env["CCP4I2_HOME"] = str(tmp_path)  # don't touch the real ~/.ccp4i2
    env["PYTHONPATH"] = os.pathsep.join([str(SERVER_ROOT), env.get("PYTHONPATH", "")]).rstrip(os.pathsep)

    # Scrub the CCP4 install: unset env vars and drop CCP4 bin dirs from PATH so
    # Class C (shutil.which(<binary>)) and Class B (env discovery) both see nothing.
    ccp4 = env.get("CCP4")
    for var in _CCP4_ENV_VARS:
        env.pop(var, None)
    if ccp4:
        env["PATH"] = os.pathsep.join(
            p for p in env.get("PATH", "").split(os.pathsep) if ccp4 not in p
        )

    result = subprocess.run(
        [sys.executable, str(PROBE)], env=env,
        capture_output=True, text=True, timeout=600,
    )
    out = (result.stdout + result.stderr).strip()
    assert result.returncode == 0, (
        "A task's validity()/runTimeValidity() requires the CCP4 install on the "
        "request path — defer it (return SUCCEEDED when absent; worker re-runs in "
        "process()) or port it to pure gemmi/numpy.\n" + out
    )
    assert "CLEAN" in result.stdout, f"Unexpected probe output:\n{out}"
