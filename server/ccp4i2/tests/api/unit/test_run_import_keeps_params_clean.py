"""Starting a job must not rewrite input_params.xml with parameters nobody set.

The last writer standing after PR #348: import_input_files_async runs at the
start of EVERY job (async_run_job calls it before execution) and re-saves
input_params.xml after stamping dbFileId/relPath on the imported inputs. It
passed exclude_unset=False, so the moment a job ran, its parameter file went
from "what the user chose" to a dump of the whole container --- every
untouched integer at zero, every empty file object serialised. Nothing in the
create/set/upload path could catch it, because those were all clean: the file
was corrupted by Run itself.

The sequence here is the New Project autogeneration verbatim (create_task,
set the list to [{}], upload a real MTZ), then the import step exactly as the
runner performs it --- create_plugin_for_job + import_input_files_async, with
a stub in place of the database handler. No CCP4 binaries needed.
"""
import json
import pathlib
import types
import xml.etree.ElementTree as ET

import pytest
from django.core.files.uploadedfile import SimpleUploadedFile
from rest_framework.test import APIClient

from ccp4i2.db import models
from ccp4i2.tests.unit.conformance import harness

MTZ = pathlib.Path(__file__).parents[3] / 'demo_data' / 'gamma' / 'gamma_Xe_mosflm.mtz'
TASK = 'aimless_pipe'


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


class _StubDbHandler:
    """The two async touchpoints the import path reaches for an already-
    uploaded file (it has a dbFileId, so no copy and no registration of a
    new file happens)."""

    async def register_input_file(self, **kwargs):
        return None

    async def find_imported_file_by_checksum(self, **kwargs):
        return None


async def _run_import(job):
    from ccp4i2.lib.async_run_job import create_plugin_for_job
    from ccp4i2.lib.async_import_files import import_input_files_async
    stub = _StubDbHandler()
    plugin = await create_plugin_for_job(job, stub)
    imported = await import_input_files_async(job, plugin, stub)
    return imported


def _configured_job(client, tmp_path):
    (tmp_path / 'p').mkdir()
    # unique per test: the async tests in this file share a database
    proj = models.Project.objects.create(name=f'runclean_{tmp_path.name}',
                                         directory=str(tmp_path / 'p'))
    r = client.post(f'/api/ccp4i2/projects/{proj.id}/create_task/',
                    data=json.dumps({'task_name': TASK}), content_type='application/json')
    assert r.status_code == 200, r.content
    # select_related: the async import path reads job.project attributes,
    # and a lazy ORM fetch is illegal in async context
    job = models.Job.objects.select_related('project').get(id=r.json()['data']['new_job']['id'])

    r1 = client.post(f'/api/ccp4i2/jobs/{job.id}/set_parameter/',
                     data=json.dumps({'object_path': f'{TASK}.container.inputData.UNMERGEDFILES',
                                      'value': [{}]}),
                     content_type='application/json')
    assert r1.status_code == 200, r1.content

    r2 = client.post(f'/api/ccp4i2/jobs/{job.id}/upload_file_param/',
                     {'file': SimpleUploadedFile('gamma_Xe_mosflm.mtz', MTZ.read_bytes()),
                      'objectPath': f'{TASK}.container.inputData.UNMERGEDFILES[0].file'},
                     format='multipart')
    assert r2.status_code == 200 and r2.json().get('success'), r2.content
    return job


def _body(job):
    return ET.parse(str(job.directory / 'input_params.xml')).getroot().find('.//ccp4i2_body')


async def test_starting_the_job_does_not_dump_the_container(client, tmp_path):
    job = await sync_wrap(_configured_job, client, tmp_path)
    before = {e.tag for e in _body(job).iter()}
    assert 'NPROC' not in before, 'file was dirty before the run --- different bug'

    await _run_import(job)

    body = _body(job)
    tags = {e.tag for e in body.iter()}
    dumped = sorted(tags & {'NPROC', 'HKLIN_REF', 'XYZIN_REF', 'FREERFLAG',
                            'SCALES_NTILEX', 'CHOOSE_SOLUTION_NO'})
    assert not dumped, (
        f'starting the job wrote parameters nobody set: {dumped} --- '
        f'an unset NPROC written as 0 blocks the job on any later reload'
    )


async def test_the_imported_file_details_still_reach_the_file(client, tmp_path):
    """The reason the importer saves at all: dbFileId/relPath stamped during
    import must survive the exclude_unset save."""
    job = await sync_wrap(_configured_job, client, tmp_path)
    await _run_import(job)

    body = _body(job)
    item = body.find('inputData/UNMERGEDFILES/CImportUnmerged/file')
    assert item is not None, 'the uploaded file vanished from input_params.xml'
    base = item.find('baseName')
    assert base is not None and base.text.strip() == 'gamma_xe_mosflm.mtz'
    dbid = item.find('dbFileId')
    assert dbid is not None and dbid.text and dbid.text.strip(), \
        'dbFileId did not survive --- the save the importer exists for'


async def sync_wrap(fn, *args):
    from asgiref.sync import sync_to_async
    return await sync_to_async(fn)(*args)
