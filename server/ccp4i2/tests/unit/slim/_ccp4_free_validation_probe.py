"""Probe: run every task's validity()/runTimeValidity() with the CCP4 install
made unavailable, and check that none of them require it.

The slim/CCP4-free API boundary is wider than "no CCP4 imports" (which the
sibling test_no_ccp4_native_imports.py covers). The /validation/ and
/run_time_validation/ endpoints instantiate task plugins and call these methods
on the API, so a method that reads CCP4 reference data (Class B) or shells out to
a CCP4 binary (Class C) breaks on the CCP4-free server — even with no CCP4 import.
That is the checkMonomeCoverage regression (PR #175): it located the monomer
library via $CLIBD_MON/$CCP4 and, with none present, reported every residue as
"no dictionary", blocking refinement at Confirm.

This probe is run by test_validation_ccp4_free.py in a subprocess whose CCP4
install has been scrubbed (env vars unset, CCP4 bin removed from PATH), so it is
meaningful on a slim interpreter AND under ccp4-python.

Exit 0 + "CLEAN" if no task needs the install; exit 1 + offenders otherwise.
"""
import logging
logging.disable(logging.CRITICAL)

import subprocess
import sys
import tempfile

import django
django.setup()

from ccp4i2 import I2_TOP
from ccp4i2.core.tasks import TASKS, get_plugin_class

# The signatures of an install-access failure (missing data file / missing binary).
INSTALL_ERRS = (FileNotFoundError, NotADirectoryError, PermissionError,
                subprocess.CalledProcessError)

failures = []

# --- Part 1: every importable task, empty container -------------------------
# Catches *unconditional* install access in validity()/runTimeValidity().
# Tasks whose plugin module does not import on this interpreter (missing pip dep
# or a CCP4-native import) yield get_plugin_class() -> None and are skipped here;
# they are a separate (Class A) concern for the import guard / dependency set.
for name in sorted(TASKS):
    try:
        cls = get_plugin_class(name)
    except Exception:
        continue
    if cls is None:
        continue
    try:
        plugin = cls(workDirectory=tempfile.mkdtemp())
    except Exception:
        continue
    for meth in ("validity", "runTimeValidity"):
        try:
            getattr(plugin, meth)()
        except INSTALL_ERRS as exc:
            failures.append("%s.%s -> %s: %s" % (name, meth, type(exc).__name__, str(exc)[:140]))
        except Exception:
            pass  # validation-logic errors (returned via CErrorReport) are fine

# --- Part 2: the monomer-coverage path with a real model --------------------
# This is *input-gated*: it only runs when a model is present, so the Part-1
# empty-container sweep would not reach it. Exercises the exact PR #175 scenario
# under a CCP4-free env — it must defer gracefully, not raise FileNotFoundError.
try:
    refmac_cls = get_plugin_class("refmac")
    if refmac_cls is not None:
        plugin = refmac_cls(workDirectory=tempfile.mkdtemp())
        model = str(I2_TOP / "demo_data" / "gamma" / "gamma_model.pdb")
        plugin.checkMonomeCoverage(model)
except INSTALL_ERRS as exc:
    failures.append("checkMonomeCoverage(model) -> %s: %s" % (type(exc).__name__, str(exc)[:140]))
except Exception:
    pass

if failures:
    print("CCP4-FREE VALIDATION FAILURES (need the CCP4 install on the request path):")
    for line in failures:
        print("  " + line)
    sys.exit(1)
print("CLEAN")
sys.exit(0)
