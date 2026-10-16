"""The interface's MetalCoord monomer picker stores its choice as a whole list.

servalcat_pipe runs MetalCoord for each code in
metalCoordPipeline.LIGAND_CODES_SELECTED; an empty list means every metal
site in the model. The picker writes the list through set_parameter, so a
list of codes must arrive intact and an empty one must clear it.

No CCP4 binaries needed.
"""
import json

import pytest
from rest_framework.test import APIClient

from ccp4i2.db import models

TASK = 'servalcat_pipe'
CODES = f'{TASK}.container.metalCoordPipeline.LIGAND_CODES_SELECTED'


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def job(client, tmp_path):
    (tmp_path / 'p').mkdir()
    project = models.Project.objects.create(name=f'metal_{tmp_path.name}',
                                            directory=str(tmp_path / 'p'))
    r = client.post(f'/api/ccp4i2/projects/{project.id}/create_task/',
                    data=json.dumps({'task_name': TASK}),
                    content_type='application/json')
    assert r.status_code == 200, r.content
    return models.Job.objects.get(id=r.json()['data']['new_job']['id'])


def _set(client, job, value):
    r = client.post(f'/api/ccp4i2/jobs/{job.id}/set_parameter/',
                    data=json.dumps({'object_path': CODES, 'value': value}),
                    content_type='application/json')
    assert r.status_code == 200, r.content


def _codes(client, job):
    r = client.get(f'/api/ccp4i2/jobs/{job.id}/container/')
    assert r.status_code == 200, r.content
    node = r.json()['data']['result']['_value']['metalCoordPipeline']['_value']['LIGAND_CODES_SELECTED']
    return [item['_value'] if isinstance(item, dict) else item for item in node['_value']]


def test_a_list_of_codes_is_stored_whole(client, job):
    _set(client, job, ['ZN', 'HEM'])
    assert _codes(client, job) == ['ZN', 'HEM']


def test_an_empty_list_clears_the_choice(client, job):
    _set(client, job, ['ZN'])
    _set(client, job, [])
    assert _codes(client, job) == []
