"""Re-pointing a plugin's working file attribute.

A pipeline tracks which file the next stage should read in an ordinary
attribute -- ``self.finalCoordinates``, ``self.F_SIGF_TOUSE`` -- and re-points
it after each stage. Plain assignment does not re-point it. A CPluginScript is
itself a CData, so ``self.x = otherFile`` is intercepted by CData.__setattr__,
which copies into the object already there instead of rebinding the name; for
a CDataFile that copy does not carry the path, so the attribute goes on
reading as the first file while the code says otherwise.

Every stage after the first is then discarded in silence. These tests pin both
halves: that the coercion really does behave this way (so the day it changes,
something says so), and that ``_useFile`` defeats it.
"""

import pytest

from ccp4i2.core.CCP4File import CDataFile
from ccp4i2.core.CCP4PluginScript import CPluginScript


@pytest.fixture
def plugin():
    return CPluginScript(name="working_attributes_test")


def coordinates(path):
    dataFile = CDataFile()
    dataFile.setFullPath(path)
    return dataFile


class TestPlainAssignment:
    def test_the_first_assignment_binds(self, plugin):
        """Nothing to coerce into while the attribute is still None."""
        plugin.finalCoordinates = None
        first = coordinates("/tmp/stage1/refined.pdb")

        plugin.finalCoordinates = first

        assert plugin.finalCoordinates is first

    def test_a_reassignment_silently_keeps_the_first_file(self, plugin):
        """The trap, pinned.

        This is not the behaviour anybody writing the assignment intends, and
        it is why `_useFile` exists. If CData's assignment semantics are ever
        changed so that this rebinds, this test should fail and be deleted --
        loudly, rather than leaving `_useFile` as cargo.
        """
        plugin.finalCoordinates = None
        plugin.finalCoordinates = coordinates("/tmp/stage1/refined.pdb")

        plugin.finalCoordinates = coordinates("/tmp/stage2/rerefined.pdb")

        assert str(plugin.finalCoordinates.fullPath).endswith("refined.pdb")
        assert "stage2" not in str(plugin.finalCoordinates.fullPath)

    def test_assigning_an_unset_file_also_does_nothing(self, plugin):
        """A second silent no-op in the same code path.

        CData skips the assignment entirely when the source has no value set,
        so a stage that failed to produce a file leaves the previous stage's
        file in place rather than clearing it.
        """
        plugin.mapToUse = None
        plugin.mapToUse = coordinates("/tmp/stage1/map.mtz")

        plugin.mapToUse = CDataFile()  # a stage that produced nothing

        assert str(plugin.mapToUse.fullPath).endswith("map.mtz")


class TestUseFile:
    def test_repoints_the_attribute(self, plugin):
        plugin.finalCoordinates = None
        plugin._useFile("finalCoordinates", coordinates("/tmp/stage1/refined.pdb"))

        second = coordinates("/tmp/stage2/rerefined.pdb")
        plugin._useFile("finalCoordinates", second)

        assert plugin.finalCoordinates is second
        assert str(plugin.finalCoordinates.fullPath).endswith("rerefined.pdb")

    def test_survives_a_whole_chain_of_stages(self, plugin):
        """What a pipeline actually does: re-point after every stage."""
        plugin.coordinatesForCoot = None
        for stage in ("dimple", "servalcat", "coot"):
            plugin._useFile(
                "coordinatesForCoot", coordinates(f"/tmp/{stage}/{stage}.pdb")
            )

        # The last stage's file, not the first -- which is what the plain
        # assignment in this loop would have left behind.
        assert str(plugin.coordinatesForCoot.fullPath).endswith("coot.pdb")

    def test_the_first_assignment_works_too(self, plugin):
        """So a file can use it uniformly rather than only after the first."""
        plugin.obsToUse = None
        first = coordinates("/tmp/input/observed.mtz")

        plugin._useFile("obsToUse", first)

        assert plugin.obsToUse is first

    def test_available_to_every_plugin(self):
        """It lives on CPluginScript, not on the one pipeline that found it."""
        assert callable(CPluginScript._useFile)
