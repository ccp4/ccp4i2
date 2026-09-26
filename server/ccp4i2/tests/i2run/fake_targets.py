"""A program run target for tests that run jobs in child processes.

Registered by dotted path (``ccp4i2.tests.i2run.fake_targets.FakeBatch``)
through ``CCP4I2_RUN_TARGETS`` in the environment; answers ``poll`` with
``CCP4I2_TEST_FAKE_BATCH_STATE`` and ``logs`` with
``CCP4I2_TEST_FAKE_BATCH_LOG``, so a test steers "the run" from outside the
processes that use the target. ``submit`` records what it was given in
``CCP4I2_TEST_FAKE_BATCH_RECORD`` (a JSON file) and never runs anything.
"""
import json
import os


class FakeBatch:
    def submit(self, tree, argv, out_dir, sizing_hint):
        record = os.environ.get("CCP4I2_TEST_FAKE_BATCH_RECORD")
        if record:
            with open(record, "w") as fh:
                json.dump({"tree": str(tree), "argv": list(argv), "out_dir": str(out_dir),
                           "sizing_hint": dict(sizing_hint)}, fh)
        return "fake-batch-1"

    def poll(self, handle):
        return os.environ.get("CCP4I2_TEST_FAKE_BATCH_STATE", "queued")

    def cancel(self, handle):
        return None

    def logs(self, handle):
        return os.environ.get("CCP4I2_TEST_FAKE_BATCH_LOG") or None
