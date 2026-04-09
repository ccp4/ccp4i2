"""Test that RESMAX stays NOT_SET and is excluded from XML."""

import tempfile
from pathlib import Path

from ccp4i2.core.tasks import get_plugin_class


def test_resmax_not_set():
    """Verify RESMAX without default stays NOT_SET and is excluded from XML."""
    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir)

        plugin_class = get_plugin_class('freerflag')
        plugin = plugin_class(parent=None, name='freerflag_test')
        plugin.workDirectory = workdir

        resmax = plugin.container.controlParameters.RESMAX
        assert not resmax.isSet(allowDefault=False)

        output_file = workdir / "test_params.xml"
        plugin.saveDataToXml(str(output_file))

        content = output_file.read_text()
        assert '<RESMAX>' not in content, "RESMAX should be excluded from XML when not set"
