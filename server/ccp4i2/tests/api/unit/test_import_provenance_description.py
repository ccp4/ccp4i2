"""Importing a file may carry a free-text provenance note.

When the client's import-provenance preference is on, the upload POST includes
a `description` form field ("where did this come from?"). It must be stored on
the file's import record (FileImport.description) -- distinct from the
auto-generated File.annotation label -- and must default to blank when the
client sends no note (the preference-off / programmatic case).

Mirrors the upload fixture in test_run_import_keeps_params_clean.py. No CCP4
binaries needed.
"""
import json
import pathlib

import pytest
from django.core.files.uploadedfile import SimpleUploadedFile
from rest_framework.test import APIClient

from ccp4i2.db import models

MTZ = pathlib.Path(__file__).parents[3] / 'demo_data' / 'gamma' / 'gamma_Xe_mosflm.mtz'
TASK = 'aimless_pipe'


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


def _upload(client, tmp_path, extra):
    (tmp_path / 'p').mkdir()
    proj = models.Project.objects.create(name=f'prov_{tmp_path.name}',
                                         directory=str(tmp_path / 'p'))
    r = client.post(f'/api/ccp4i2/projects/{proj.id}/create_task/',
                    data=json.dumps({'task_name': TASK}),
                    content_type='application/json')
    assert r.status_code == 200, r.content
    job = models.Job.objects.select_related('project').get(
        id=r.json()['data']['new_job']['id'])

    r1 = client.post(f'/api/ccp4i2/jobs/{job.id}/set_parameter/',
                     data=json.dumps({
                         'object_path': f'{TASK}.container.inputData.UNMERGEDFILES',
                         'value': [{}]}),
                     content_type='application/json')
    assert r1.status_code == 200, r1.content

    payload = {
        'file': SimpleUploadedFile('gamma_Xe_mosflm.mtz', MTZ.read_bytes()),
        'objectPath': f'{TASK}.container.inputData.UNMERGEDFILES[0].file',
    }
    payload.update(extra)
    r2 = client.post(f'/api/ccp4i2/jobs/{job.id}/upload_file_param/',
                     payload, format='multipart')
    assert r2.status_code == 200 and r2.json().get('success'), r2.content
    return job


def _import_record(job):
    imports = list(models.FileImport.objects.filter(file__job=job))
    assert len(imports) == 1, f'expected one import, got {len(imports)}'
    return imports[0]


def test_description_is_stored_on_the_import_record(client, tmp_path):
    note = 'Collected at DLS i04, processed with xia2/DIALS; from Alice.'
    job = _upload(client, tmp_path, {'description': note})
    fi = _import_record(job)
    assert fi.description == note
    # It must NOT have clobbered the auto-generated File.annotation label.
    assert fi.file.annotation and fi.file.annotation != note


def test_description_defaults_blank_without_a_note(client, tmp_path):
    """The preference-off / programmatic path sends no `description`; the
    import record's description is then simply blank."""
    job = _upload(client, tmp_path, {})
    fi = _import_record(job)
    assert fi.description == ''


def test_blank_description_is_stripped(client, tmp_path):
    """A whitespace-only note is treated as no note, not stored verbatim."""
    job = _upload(client, tmp_path, {'description': '   '})
    fi = _import_record(job)
    assert fi.description == ''
