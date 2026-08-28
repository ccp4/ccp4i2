"""
Configuring a PHIL-driven task through the API.

Every API test in the tree drives `aimless_pipe`, whose container comes from a
def.xml. PHIL tasks --- phasertng_picard and riker, xia2_dials and xds,
editbfac --- build theirs by a different route, and nothing exercised that
route through the REST layer at all. (`phaser_MR` looks PHIL-ish for its
keywords but stayed on the static def.xml model, so it is not the test here.)

The combination matters because setting a parameter through the API
constructs a plugin every time, then walks the container by path. A change to
construction can hold for one route and fail for the other, and this is the
cheap place to find out: configuration only, no execution, no CCP4 binaries.
"""
import pytest

pytest.importorskip("libtbx.phil", reason="needs libtbx (CCP4/cctbx)")

from ..base import APITestBase  # noqa: E402


class TestPhilTaskConfiguration(APITestBase):
    """phasertng_picard: PHIL in, container out, addressed over REST."""

    task_name = "phasertng_picard"

    def test_a_phil_task_can_be_created_at_all(self):
        self.create_project("phil_config_create")
        job = self.create_job()
        assert job['id']

    def test_its_container_is_reachable_over_the_api(self):
        self.create_project("phil_config_container")
        self.create_job()
        response = self.client.get(f'{self.API_PREFIX}/jobs/{self.job_id}/container/')
        assert response.status_code == 200, response.content
        top = response.json()["data"]["result"]["_value"]
        assert "inputData" in top
        assert "controlParameters" in top

    def test_a_parameter_can_be_set_and_read_back(self):
        """The path that constructs a plugin on every call."""
        self.create_project("phil_config_setparam")
        self.create_job()
        response = self.client.get(f'{self.API_PREFIX}/jobs/{self.job_id}/container/')
        control = response.json()["data"]["result"]["_value"]["controlParameters"]["_value"]

        name = next((n for n, v in control.items()
                     if v.get("_class") in ("CInt", "CFloat", "CString", "CBoolean")), None)
        assert name is not None, (
            "no simple control parameter found; with libtbx present this task "
            f"has PHIL_EXPERT_LEVEL among {sorted(control)[:5]}")

        original = control[name]["_value"]
        value = {"CBoolean": True, "CInt": 3, "CFloat": 1.5}.get(
            control[name]["_class"], "test")
        self.set_param(f"controlParameters.{name}", value)

        after = self.client.get(f'{self.API_PREFIX}/jobs/{self.job_id}/container/')
        now = after.json()["data"]["result"]["_value"]["controlParameters"]["_value"][name]["_value"]
        assert str(now) != str(original) or str(now) == str(value), \
            f"{name} did not take the value: {original!r} -> {now!r}"

    def test_a_file_can_be_uploaded_to_a_phil_task(self, gamma_model_pdb):
        """Upload resolves a path through a PHIL-built container."""
        self.create_project("phil_config_upload")
        self.create_job()
        self.upload_file("inputData.XYZIN", gamma_model_pdb)

        from ccp4i2.db import models
        records = list(models.File.objects.filter(job__id=self.job_id))
        assert records, "the upload produced no file record"
        assert any(r.path.exists() for r in records), \
            "the file record points at nothing on disk"

    def test_validation_runs_against_a_phil_container(self):
        """validity() walks the container; a PHIL tree must survive the walk."""
        self.create_project("phil_config_validation")
        self.create_job()
        result = self.get_validation()
        assert result is not None
