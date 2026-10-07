"""A Phaser follow-on given an ASU contents file uses it (#678).

The context filled ASUFILE, but COMP_BY stayed at its default (average
solvent content), so the file was ignored and its field hidden.
"""

from types import SimpleNamespace

import pytest


from ccp4i2.wrappers.phaser_phil.script.phaser_shims import AsuCompositionFromContext


class Param:
    def __init__(self, value=None):
        self.value = value

    def isSet(self):
        return self.value is not None

    def set(self, value):
        self.value = value

    def __str__(self):
        return str(self.value)


def plugin(comp_by, asu_path):
    job = AsuCompositionFromContext()
    job.container = SimpleNamespace(inputData=SimpleNamespace(
        COMP_BY=Param(comp_by), ASUFILE=Param(asu_path)))
    return job


def test_default_composition_switches_to_the_asu_file_given():
    job = plugin("DEFAULT", "/p/job_1/ASUCONTENTFILE.asu.xml")
    job.contextApplied()
    assert str(job.container.inputData.COMP_BY) == "ASU"


@pytest.mark.parametrize("chosen", ["SEQUENCES", "MW", "SOLVENT", "ASU"])
def test_a_chosen_composition_is_left_alone(chosen):
    job = plugin(chosen, "/p/job_1/ASUCONTENTFILE.asu.xml")
    job.contextApplied()
    assert str(job.container.inputData.COMP_BY) == chosen


def test_no_asu_file_no_change():
    job = plugin("DEFAULT", None)
    job.contextApplied()
    assert str(job.container.inputData.COMP_BY) == "DEFAULT"
