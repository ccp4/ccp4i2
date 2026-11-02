"""Pointing a file object at another file stops it naming the old record.

getFullPath() prefers dbFileId to the path, so a file object that held a
database file -- autopopulated from an earlier job, or copied with copyData --
kept resolving to the old file after setFullPath named a new one. molrep_pipe's
second MOLREP read the unreindexed data that way (#785); an upload recorded the
old file's checksum (#817).
"""
import uuid

from ccp4i2.core.CCP4ModelData import CPdbDataFile

OLD_ID = str(uuid.uuid4())


def _held(path):
    f = CPdbDataFile()
    f.setFullPath(str(path))
    f.dbFileId.set(OLD_ID)
    return f


def test_another_file_drops_the_old_id(tmp_path):
    old = tmp_path / "old.pdb"
    new = tmp_path / "new.pdb"
    old.write_text("END\n")
    new.write_text("END\n")
    f = _held(old)
    f.setFullPath(str(new))
    assert not f.dbFileId.isSet()
    assert f.getFullPath() == str(new)


def test_the_same_file_keeps_its_id(tmp_path):
    path = tmp_path / "model.pdb"
    path.write_text("END\n")
    f = _held(path)
    f.setFullPath(str(path))
    assert str(f.dbFileId) .replace("-", "") == OLD_ID.replace("-", "")


def test_a_file_without_an_id_is_unaffected(tmp_path):
    path = tmp_path / "model.pdb"
    path.write_text("END\n")
    f = CPdbDataFile()
    f.setFullPath(str(path))
    assert not f.dbFileId.isSet()
    assert f.getFullPath() == str(path)
