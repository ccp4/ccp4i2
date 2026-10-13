"""Uploading a file the project already holds uses the copy it has.

The same unmerged file brought into a second job used to be copied into
CCP4_IMPORTED_FILES again, and again for every job after --- the copies people
then delete by hand to get their disk back. An upload identical to an earlier
import (same bytes, type, content flag and subtype) now points the parameter at
that earlier file and keeps no second copy.

Also guards the bug that hid every such match: a slot autopopulated with a
previous job's file went on resolving to that file after the upload, because
getFullPath() prefers dbFileId to the path, so the checksum recorded for the
import was the old file's.

No CCP4 binaries needed.
"""
import hashlib
import json
import pathlib

import pytest
from django.core.files.uploadedfile import SimpleUploadedFile
from rest_framework.test import APIClient

from ccp4i2.db import models

DEMO = pathlib.Path(__file__).parents[3] / 'demo_data' / 'gamma'
MTZ = DEMO / 'gamma_Xe_mosflm.mtz'
OTHER = DEMO / 'HKLOUT_unmerged.mtz'
TASK = 'aimless_pipe'
SLOT = f'{TASK}.container.inputData.UNMERGEDFILES[0].file'


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def project(tmp_path):
    (tmp_path / 'p').mkdir()
    return models.Project.objects.create(name=f'reuse_{tmp_path.name}',
                                         directory=str(tmp_path / 'p'))


def _new_job(client, project):
    r = client.post(f'/api/ccp4i2/projects/{project.id}/create_task/',
                    data=json.dumps({'task_name': TASK}),
                    content_type='application/json')
    assert r.status_code == 200, r.content
    job = models.Job.objects.get(id=r.json()['data']['new_job']['id'])
    r = client.post(f'/api/ccp4i2/jobs/{job.id}/set_parameter/',
                    data=json.dumps({'object_path': f'{TASK}.container.inputData.UNMERGEDFILES',
                                     'value': [{}]}),
                    content_type='application/json')
    assert r.status_code == 200, r.content
    return job


def _upload(client, job, source=MTZ):
    r = client.post(f'/api/ccp4i2/jobs/{job.id}/upload_file_param/',
                    {'file': SimpleUploadedFile(source.name, source.read_bytes()),
                     'objectPath': SLOT},
                    format='multipart')
    assert r.status_code == 200 and r.json().get('success'), r.content
    return r.json()['data']


def _imported_files(project):
    return sorted(p.name for p in (pathlib.Path(project.directory) / 'CCP4_IMPORTED_FILES').iterdir())


def test_a_second_upload_of_the_same_file_uses_the_first_copy(client, project):
    first_job = _new_job(client, project)
    _upload(client, first_job)
    [first] = models.File.objects.filter(job=first_job, fileimport__isnull=False)

    second_job = _new_job(client, project)
    data = _upload(client, second_job)

    assert data['reused']['uuid'] == str(first.uuid)
    assert data['reused']['jobNumber'] == first_job.number
    assert _imported_files(project) == [first.name]
    assert not models.File.objects.filter(job=second_job, fileimport__isnull=False).exists()
    held = data['updated_item']['_value']['dbFileId']['_value']
    assert held.replace('-', '') == str(first.uuid).replace('-', '')


def test_the_checksum_recorded_is_the_uploaded_file_s(client, project):
    """The slot already names another job's file, as autopopulation leaves it."""
    owner = _new_job(client, project)
    owner.directory.mkdir(parents=True, exist_ok=True)
    (owner.directory / OTHER.name).write_bytes(OTHER.read_bytes())
    earlier = models.File.objects.create(
        job=owner, name=OTHER.name, directory=models.File.Directory.JOB_DIR,
        type=models.FileType.objects.get_or_create(name='application/CCP4-unmerged-experimental')[0],
        annotation='an earlier output',
    )
    job = _new_job(client, project)
    r = client.post(f'/api/ccp4i2/jobs/{job.id}/set_parameter/',
                    data=json.dumps({'object_path': SLOT, 'value': {
                        'project': str(project.uuid).replace('-', ''),
                        'baseName': earlier.name,
                        'relPath': owner.rel_path,
                        'dbFileId': str(earlier.uuid).replace('-', ''),
                    }}),
                    content_type='application/json')
    assert r.status_code == 200, r.content

    _upload(client, job)

    [record] = models.FileImport.objects.filter(file__job=job)
    assert record.checksum == hashlib.md5(MTZ.read_bytes()).hexdigest()
