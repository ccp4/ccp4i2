"""Imported files that have gone from the project, and getting them back.

Every file a job imports is copied into the project's CCP4_IMPORTED_FILES, so
that the project carries its own inputs and nothing outside it can change them.
People delete those copies anyway --- most often a second copy of a large
unmerged file, to save space --- and the job that names one then fails in
whatever program first tries to read it, with an error that says nothing about
why.

Nothing here can stop that, and the original outside the project is no safer:
whoever deletes the project's copy can delete the original too. What it can do
is notice before the job runs, and

  - put the copy back from where it was imported from, when that original is
    still there and still holds exactly the bytes that were imported (its
    checksum is FileImport.checksum), and say that it did;
  - otherwise say plainly which file has gone, where it came from, and what to
    do about it, rather than letting the program fail.

FileImport.name holds the source as the importer knew it: an absolute path for
a desktop import by path, for i2run, and for projects from Qt-i2 (whose
ImportFiles.SourceFilename was the full path); only a basename for a browser
upload, whose original location the server never learns. Only an absolute
path can be restored from.

A derived import --- an MTZ split down to the columns one parameter asked for
--- is not byte-identical to its source, so it is reported, never rebuilt here:
rebuilding it would mean guessing the columns.
"""

import enum
import hashlib
import logging
import os
import tempfile
import shutil
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.db import models

logger = logging.getLogger(f"ccp4i2:{__name__}")

#: Error codes reported on the parameter that names the missing file.
CODE_WILL_BE_RESTORED = 110
CODE_CANNOT_BE_RESTORED = 111


class Missing(enum.Enum):
    #: The original is where it was imported from, byte for byte.
    RESTORABLE = "restorable"
    #: The original is there, but the project held something derived from it.
    DERIVED = "derived"
    #: Something is at the original path, but not what was imported.
    CHANGED = "changed"
    #: Nothing is at the original path.
    GONE = "gone"
    #: Uploaded through a browser, or fetched: no path to look in.
    UNKNOWN = "unknown"


@dataclass
class MissingImport:
    file: models.File
    source: Optional[Path]
    state: Missing

    @property
    def restorable(self) -> bool:
        return self.state is Missing.RESTORABLE

    def describe(self) -> str:
        gone = (
            f"{self.file.name}, imported by job {self.file.job.number}, has been "
            f"deleted from the project's CCP4_IMPORTED_FILES."
        )
        source = self.source
        if self.state is Missing.RESTORABLE:
            return (f"{gone} The file it was imported from, {source}, is unchanged, "
                    f"so it will be copied back when the job runs.")
        if self.state is Missing.DERIVED:
            return (f"{gone} It held data taken from {source}, which is still there, "
                    f"but the selection cannot be rebuilt automatically: import "
                    f"{source.name} again.")
        if self.state is Missing.CHANGED:
            return (f"{gone} It was imported from {source}, which has changed since, "
                    f"so it cannot be put back: import the file again.")
        if self.state is Missing.GONE:
            return (f"{gone} It was imported from {source}, which is no longer there "
                    f"either: import the file again.")
        return (f"{gone} It came from {self.file.fileimport.name or 'an upload'}, "
                f"and the project does not know where that is now: import the file "
                f"again.")


# Checksumming a large unmerged file is the only slow step here, and validation
# is polled while a task is edited. Key on what would change if the bytes did.
_md5_cache: Dict[Tuple[str, int, int], str] = {}


def _md5(path: Path) -> str:
    try:
        stat = path.stat()
        key = (str(path), stat.st_size, stat.st_mtime_ns)
        if key not in _md5_cache:
            digest = hashlib.md5()
            with open(path, "rb") as handle:
                for chunk in iter(lambda: handle.read(1 << 20), b""):
                    digest.update(chunk)
            _md5_cache[key] = digest.hexdigest()
        return _md5_cache[key]
    except OSError as err:
        logger.warning("Could not checksum %s: %s", path, err)
        return ""


def _source_of(file_import: models.FileImport) -> Optional[Path]:
    if not file_import.name:
        return None
    source = Path(file_import.name)
    return source if source.is_absolute() else None


def classify(file_record: models.File) -> Optional[MissingImport]:
    """How an imported file stands, or None if it is where it should be (or
    is not an import at all)."""
    if file_record.directory != models.File.Directory.IMPORT_DIR:
        return None
    try:
        file_import = file_record.fileimport
    except models.FileImport.DoesNotExist:
        return None
    if file_record.path.exists():
        return None

    source = _source_of(file_import)
    if source is None:
        state = Missing.UNKNOWN
    elif not source.is_file():
        state = Missing.GONE
    else:
        digest = _md5(source)
        if digest and digest == file_import.checksum:
            state = Missing.RESTORABLE
        elif digest and digest == file_import.source_checksum:
            state = Missing.DERIVED
        else:
            state = Missing.CHANGED
    return MissingImport(file=file_record, source=source, state=state)


def missing_imports(file_ids: Iterable[uuid.UUID]) -> List[MissingImport]:
    """Of these files, the imported ones no longer on disk."""
    records = (
        models.File.objects.filter(
            uuid__in=list(file_ids),
            directory=models.File.Directory.IMPORT_DIR,
            fileimport__isnull=False,
        )
        .select_related("fileimport", "job", "job__project")
    )
    found = []
    for record in records:
        state = classify(record)
        if state is not None:
            found.append(state)
    return found


def restore(missing: MissingImport) -> bool:
    """Copy a restorable import back into place. Verifies what it wrote and
    leaves nothing half-written behind; never raises."""
    if not missing.restorable:
        return False
    target = missing.file.path
    expected = missing.file.fileimport.checksum
    temporary = None
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        handle, temporary = tempfile.mkstemp(
            dir=target.parent, prefix=f".{target.name}.", suffix=".restoring"
        )
        os.close(handle)
        shutil.copyfile(missing.source, temporary)
        # mkstemp makes it owner-only; the original's mode is the better guess.
        shutil.copymode(missing.source, temporary)
        digest = hashlib.md5()
        with open(temporary, "rb") as written:
            for chunk in iter(lambda: written.read(1 << 20), b""):
                digest.update(chunk)
        if digest.hexdigest() != expected:
            logger.warning("Copy of %s did not match what was imported; not restoring",
                           missing.source)
            return False
        os.replace(temporary, target)
        temporary = None
        logger.warning("Restored imported file %s from %s", target, missing.source)
        return True
    except OSError as err:
        logger.warning("Could not restore %s from %s: %s", target, missing.source, err)
        return False
    finally:
        if temporary is not None:
            try:
                os.unlink(temporary)
            except OSError:
                pass


# ---------------------------------------------------------------------------
# A job's container


def _input_files(container) -> list:
    """Every file the job reads: all sections but outputData (the same walk
    async_import_files makes)."""
    if hasattr(container, "inputData"):
        found = []
        for section in container.children():
            if section.objectName() == "outputData":
                continue
            if hasattr(section, "find_all_files"):
                found.extend(section.find_all_files())
        return found
    return container.find_all_files()


def _db_file_id(file_obj) -> Optional[uuid.UUID]:
    db_file_id = getattr(file_obj, "dbFileId", None)
    if db_file_id is None or not db_file_id.isSet():
        return None
    try:
        return uuid.UUID(str(db_file_id).strip())
    except ValueError:
        return None


def missing_imports_in(container) -> List[Tuple[object, MissingImport]]:
    """(file object, what has happened to it) for each input of this
    container whose imported copy has gone."""
    by_id: Dict[uuid.UUID, list] = {}
    try:
        for file_obj in _input_files(container):
            file_id = _db_file_id(file_obj)
            if file_id is not None:
                by_id.setdefault(file_id, []).append(file_obj)
        if not by_id:
            return []
        return [
            (file_obj, missing)
            for missing in missing_imports(by_id)
            for file_obj in by_id[missing.file.uuid]
        ]
    except Exception:
        # Advisory: never let this check break validation or a run.
        logger.exception("Could not check for deleted imported files")
        return []


def report_missing_imports(container, error, task_name: str) -> None:
    """Add to ``error`` (a CErrorReport) one entry per input whose imported
    copy has gone: a warning if it will be put back, an error if not. Replaces
    the bare "File does not exist" a mustExist file would otherwise get, which
    says nothing about why or what to do."""
    for file_obj, missing in missing_imports_in(container):
        name = file_obj.object_path()
        error._errors = [
            each for each in error._errors
            if not (each.get("name") == name and each.get("code") == 102)
        ]
        error.append(
            klass=task_name,
            code=CODE_WILL_BE_RESTORED if missing.restorable else CODE_CANNOT_BE_RESTORED,
            details=missing.describe(),
            name=name,
            severity=(CCP4ErrorHandling.SEVERITY_WARNING if missing.restorable
                      else CCP4ErrorHandling.SEVERITY_ERROR),
        )


def restore_missing_imports(container) -> List[MissingImport]:
    """Put back every restorable import this container's inputs name. Returns
    those that are still missing afterwards, for the caller to report."""
    still_missing = []
    seen = set()
    for _file_obj, missing in missing_imports_in(container):
        if missing.file.uuid in seen:
            continue
        seen.add(missing.file.uuid)
        if not restore(missing):
            logger.warning("Imported file is missing: %s", missing.describe())
            still_missing.append(missing)
    return still_missing
