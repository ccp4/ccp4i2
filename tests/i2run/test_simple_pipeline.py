"""
Simple test to debug signal-driven pipeline execution.

This mimics aimless_pipe but with minimal complexity.
"""
from pathlib import Path
from .utils import i2run, demoData
import asyncio


def test_simple_signal_chain():
    """Test a simple single-step pipeline to verify basic signal chain functionality."""

    # Use aimless_pipe as a simple, proven wrapper
    # (pointless has a bug with empty UNMERGEDFILES that triggers UnboundLocalError)
    mtz = demoData("gamma", "gamma_native.mtz")
    args = ["aimless_pipe", "--UNMERGEDFILES", f"file={mtz}"]

    with i2run(args) as job:  # Use random project name to avoid FileExistsError
        # Check that aimless completed and produced expected outputs
        expected_outputs = [
            "FREERFLAG.mtz",
            "HKLOUT_0-observed_data_asIMEAN.mtz",
            "HKLOUT_0-observed_data.mtz",
            "HKLOUT_unmerged.mtz",
        ]
        for output in expected_outputs:
            assert (job / output).exists(), f"aimless did not create {output}"
        print(f"✓ aimless_pipe completed successfully with all expected outputs")
