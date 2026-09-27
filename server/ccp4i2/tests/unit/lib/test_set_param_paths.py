"""Every path form a caller may write names the same parameter (Materia, 2026-09-27)."""
import pytest

from ccp4i2.lib.utils.parameters.set_param import normalize_object_path


@pytest.mark.parametrize("given", [
    "prosmart_refmac.container.inputData.XYZIN",
    "prosmart_refmac.inputData.XYZIN",
    "container.inputData.XYZIN",
    "inputData.XYZIN",
])
def test_all_documented_forms_normalise_to_the_same_path(given):
    assert normalize_object_path(given, "prosmart_refmac") == "prosmart_refmac.inputData.XYZIN"


def test_deeper_paths_and_other_sections():
    assert normalize_object_path("controlParameters.LOCAL_CPUS", "pandda_campaign") == "pandda_campaign.controlParameters.LOCAL_CPUS"
    assert normalize_object_path("inputData.DATASETS[0].DTAG", "pandda_campaign") == "pandda_campaign.inputData.DATASETS[0].DTAG"
    assert normalize_object_path("container.controlParameters.RUN_MODE", "pandda_campaign") == "pandda_campaign.controlParameters.RUN_MODE"


def test_without_a_task_name_the_old_behaviour_holds():
    assert normalize_object_path("prosmart_refmac.container.inputData.XYZIN") == "prosmart_refmac.inputData.XYZIN"
    assert normalize_object_path("inputData.XYZIN") == "inputData.XYZIN"
