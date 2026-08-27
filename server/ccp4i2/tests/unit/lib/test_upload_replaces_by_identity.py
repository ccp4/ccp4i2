"""
Adding a file to a list slot must not delete another row's file.

Reconstructed from job 66 of GammaBySAD. The user added unmerged data from
gamma, then from MDM2, deleted a row, and added MDM2 back. The gamma row then
read "File does not exist", its dbFileId matched no record, and
gamma_xe_mosflm_1.mtz was gone from CCP4_IMPORTED_FILES --- while both MDM2
files remained.

    1771  mdm2_unmerged.mtz     job_param_name = UNMERGEDFILES[1].file
    1773  mdm2_unmerged_1.mtz   job_param_name = UNMERGEDFILES[2].file

`job_param_name` is *positional*. The upload path deleted every File record
sharing the slot's name, so the re-add landed at [1], matched gamma's record
by name, and unlinked gamma's file. A position identifies a slot; it does not
identify a file, and any insert or removal renumbers every slot after it.

Two independent guards, either of which prevents the loss:

  - replace what the slot *points at* (its dbFileId), not what shares its name
  - never delete a file another parameter in the same job still refers to

Pure Python -- no CCP4 binaries needed.
"""
import uuid

import pytest

from ccp4i2.lib.utils.files.upload_param import (
    _file_currently_held_by, _is_referenced_elsewhere,
)


class _FakeFile:
    """Just enough of a CDataFile for the identity helpers."""
    def __init__(self, file_id=None):
        self.dbFileId = file_id


class _FakeContainer:
    def __init__(self, files):
        self._files = files

    def find_all_files(self):
        return self._files


@pytest.fixture(autouse=True)
def _value_dict(monkeypatch):
    """value_dict_for_object, without needing a real container."""
    monkeypatch.setattr(
        'ccp4i2.lib.utils.parameters.value_dict.value_dict_for_object',
        lambda obj: {'dbFileId': getattr(obj, 'dbFileId', None)},
    )


# --- the guard that would have prevented the loss on its own ----------------

def test_a_file_referenced_by_another_row_is_not_deleted():
    """The gamma row still pointed at it when the MDM2 re-add tried to."""
    shared = str(uuid.uuid4())
    gamma = _FakeFile(shared)
    slot_being_written = _FakeFile(shared)
    container = _FakeContainer([gamma, slot_being_written])

    assert _is_referenced_elsewhere(container, shared, slot_being_written)


def test_a_file_referenced_by_nobody_else_may_go():
    only = _FakeFile(str(uuid.uuid4()))
    other = _FakeFile(str(uuid.uuid4()))
    container = _FakeContainer([only, other])

    assert not _is_referenced_elsewhere(container, only.dbFileId, only)


def test_the_slot_itself_does_not_count_as_another_reference():
    """Otherwise nothing could ever be replaced."""
    target = _FakeFile(str(uuid.uuid4()))
    container = _FakeContainer([target])
    assert not _is_referenced_elsewhere(container, target.dbFileId, target)


def test_an_unanswerable_question_does_not_authorise_deletion():
    """False means "go ahead and unlink", so it is the wrong way to fail."""
    class Hostile:
        def find_all_files(self):
            raise RuntimeError('no')

    assert _is_referenced_elsewhere(Hostile(), 'x', object()) is True


def test_an_unreadable_sibling_does_not_authorise_deletion():
    """It might be the reference that matters."""
    class Opaque:
        @property
        def dbFileId(self):
            raise RuntimeError('no')

    container = _FakeContainer([Opaque()])
    assert _is_referenced_elsewhere(container, 'x', object()) is True


# --- replacing by identity rather than by slot name -------------------------

def test_an_empty_slot_holds_nothing_to_delete():
    assert list(_file_currently_held_by(_FakeFile(None), job=object())) == []


def test_a_malformed_file_id_is_not_treated_as_a_match():
    """A stale or garbled id must not select some other record."""
    result = _file_currently_held_by(_FakeFile('not-a-uuid'), job=object())
    assert list(result) == []


def test_a_dict_valued_file_id_is_not_treated_as_a_match():
    """The serialiser sometimes yields a dict for an unset file."""
    result = _file_currently_held_by(_FakeFile({'unset': True}), job=object())
    assert list(result) == []
