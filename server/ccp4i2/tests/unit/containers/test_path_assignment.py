"""Assigning a path to a file parameter sets the file, whether the path is a
str or a pathlib.Path.

Only str was coerced: ``out.XYZOUT = Path(...)`` replaced the CPdbDataFile
with a bare PosixPath, silently, so the output could not be annotated and
was never gleaned (SliceNDice, whose outputs were assigned that way).
"""
from pathlib import Path

import pytest

from ccp4i2.core.tasks import get_plugin_class


@pytest.mark.parametrize("value", ["/tmp/model.pdb", Path("/tmp/model.pdb")])
def test_path_sets_the_file(value, tmp_path):
    plugin = get_plugin_class("slicendice")(workDirectory=str(tmp_path), name="t")
    out = plugin.container.outputData
    out.XYZOUT = value
    assert type(out.XYZOUT).__name__ == "CPdbDataFile"
    # An output file is kept by name, relative to the job: the same either way.
    assert out.XYZOUT.getFullPath().endswith("model.pdb")
