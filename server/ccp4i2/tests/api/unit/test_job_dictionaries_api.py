"""Which dictionaries belong with a job: its own files and its inputs, a
subjob inheriting its pipeline's, and nothing from anywhere else."""

import pytest
from rest_framework.test import APIClient

from ccp4i2.db import models

API = "/api/ccp4i2"
DICT, PDB = "application/refmac-dictionary", "chemical/x-pdb"


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def project(tmp_path):
    directory = tmp_path / "dicts"
    directory.mkdir()
    return models.Project.objects.create(name="dicts", directory=str(directory))


def make_job(project, number, task="servalcat_pipe", parent=None):
    return models.Job.objects.create(project=project, number=number, title=task, task_name=task,
                                     status=models.Job.Status.FINISHED, parent=parent)


def make_file(job, name, type_name, directory=models.File.Directory.JOB_DIR):
    file_type, _ = models.FileType.objects.get_or_create(name=type_name, defaults={"description": type_name})
    return models.File.objects.create(job=job, name=name, directory=directory, type=file_type,
                                      job_param_name=name.split(".")[0])


def use_as_input(job, the_file, param="DICT_LIST[0]"):
    return models.FileUse.objects.create(job=job, file=the_file, role=models.FileUse.Role.IN, job_param_name=param)


def listing(client, url):
    response = client.get(url)
    assert response.status_code == 200, response.content
    return [(f["name"], f["role"]) for f in response.json()["data"]]


def test_a_jobs_own_dictionary_and_its_input_dictionary(client, project):
    acedrg = make_job(project, "1", "LidiaAcedrgNew")
    lig = make_file(acedrg, "LIG.cif", DICT)
    refine = make_job(project, "2")
    use_as_input(refine, lig)
    assert listing(client, f"{API}/jobs/{acedrg.id}/dictionaries/") == [("LIG.cif", "own")]
    assert listing(client, f"{API}/jobs/{refine.id}/dictionaries/") == [("LIG.cif", "input")]


def test_another_jobs_ligand_of_the_same_name_is_not_offered(client, project):
    """Two ligands both called LIG in one project: each refinement sees only its own."""
    first, second = make_job(project, "1", "LidiaAcedrgNew"), make_job(project, "2", "LidiaAcedrgNew")
    lig_a, lig_b = make_file(first, "LIG_a.cif", DICT), make_file(second, "LIG_b.cif", DICT)
    refine_a, refine_b = make_job(project, "3"), make_job(project, "4")
    use_as_input(refine_a, lig_a)
    use_as_input(refine_b, lig_b)
    assert listing(client, f"{API}/jobs/{refine_a.id}/dictionaries/") == [("LIG_a.cif", "input")]
    assert listing(client, f"{API}/jobs/{refine_b.id}/dictionaries/") == [("LIG_b.cif", "input")]


def test_a_job_with_no_dictionary_gets_none(client, project):
    assert listing(client, f"{API}/jobs/{make_job(project, '1').id}/dictionaries/") == []


def test_a_subjob_inherits_its_pipelines_dictionaries(client, project):
    pipeline = make_job(project, "1", "SubstituteLigand")
    make_file(pipeline, "DICTOUT.cif", DICT)
    subjob = make_job(project, "1.2", "servalcat", parent=pipeline)
    assert listing(client, f"{API}/jobs/{subjob.id}/dictionaries/") == [("DICTOUT.cif", "own")]


def test_own_and_input_are_both_listed_without_duplicates(client, project):
    job = make_job(project, "1", "SubstituteLigand")
    own = make_file(job, "DICTOUT.cif", DICT)
    other = make_file(make_job(project, "2", "LidiaAcedrgNew"), "DRG.cif", DICT)
    use_as_input(job, other)
    use_as_input(job, own, param="DICTIN")  # the same file also recorded as an input
    assert listing(client, f"{API}/jobs/{job.id}/dictionaries/") == [("DICTOUT.cif", "own"), ("DRG.cif", "input")]


def test_a_coordinate_files_companions_are_its_jobs_dictionaries(client, project):
    refine = make_job(project, "1")
    lig = make_file(make_job(project, "2", "LidiaAcedrgNew"), "LIG.cif", DICT)
    use_as_input(refine, lig)
    xyzout = make_file(refine, "XYZOUT.cif", PDB)
    assert listing(client, f"{API}/files/{xyzout.id}/companion_dictionaries/") == [("LIG.cif", "input")]
    assert client.get(f"{API}/files/999999/companion_dictionaries/").status_code == 404
