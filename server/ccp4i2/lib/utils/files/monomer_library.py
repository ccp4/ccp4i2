"""Find a CCP4 monomer library entry, locally or from the published library.

The task interfaces draw library monomers (Make Link's "Monomer library"
source) from the entry's CIF, parsed on the server. Desktop installs read it
from the CCP4 suite. A served deployment whose API server carries no CCP4 --
the slim server, where CCP4 lives only in the worker that runs jobs -- has no
library on disk, and every residue came back as "not in the library".

The library is published at https://github.com/MonomerLibrary/monomers with
the same layout and a tag per CCP4 release, so when no local copy exists the
entry is fetched from there at a pinned tag and kept in a small cache.

Order of lookup:
  1. ``$CLIBD_MON`` (the CCP4 variable for the library directory)
  2. ``$CCP4/lib/data/monomers``
  3. the published library at ``CCP4I2_MONOMER_LIBRARY_REF``, cached

Stdlib only (``urllib``): this runs on the slim server too.
"""

import logging
import os
import re
import tempfile
import urllib.error
import urllib.request
from pathlib import Path
from typing import Optional

logger = logging.getLogger(f"ccp4i2:{__name__}")

#: The tag the remote fallback reads. ccp4-9.0.015 is the library that ships
#: with the validated CCP4 build (ccp4-20260904), so what the interface draws
#: matches what a job on that build will use.
DEFAULT_REF = "ccp4-9.0.015"
REMOTE_BASE = "https://raw.githubusercontent.com/MonomerLibrary/monomers"
USER_AGENT = "ccp4i2 (monomer library)"
TIMEOUT = 15

#: Monomer codes: letters, digits, and the few punctuation characters the
#: PDB chemical component dictionary uses. Anything else is never a code and
#: must not reach a path or a URL.
_CODE_RE = re.compile(r"^[A-Z0-9_-]{1,8}$")

#: Codes that are reserved file names on Windows are stored as CODE_CODE.cif
#: (the library's windows_reserved_words.txt).
_WINDOWS_RESERVED = frozenset({"CON", "PRN", "AUX", "NUL"})


class MonomerLibraryError(Exception):
    """A lookup that cannot be answered, with the HTTP status it earns."""

    def __init__(self, status, message):
        super().__init__(message)
        self.status = status


def normalise_code(code) -> str:
    """Upper-cased code, or a 400 for anything that is not a monomer code."""
    code = str(code or "").strip().upper()
    if not _CODE_RE.match(code):
        raise MonomerLibraryError(400, f"Not a monomer code: {code!r}")
    return code


def relative_path(code: str) -> str:
    """Where ``code`` lives inside the library, e.g. ``c/CYS.cif``."""
    file_name = f"{code}_{code}.cif" if code in _WINDOWS_RESERVED else f"{code}.cif"
    return f"{code[0].lower()}/{file_name}"


def local_library_dir() -> Optional[Path]:
    """The library directory of a local CCP4 install, if there is one."""
    clibd_mon = os.environ.get("CLIBD_MON", "")
    if clibd_mon and Path(clibd_mon).is_dir():
        return Path(clibd_mon)
    ccp4 = os.environ.get("CCP4", "")
    if ccp4 and (Path(ccp4) / "lib" / "data" / "monomers").is_dir():
        return Path(ccp4) / "lib" / "data" / "monomers"
    return None


def library_ref() -> str:
    return os.environ.get("CCP4I2_MONOMER_LIBRARY_REF", "") or DEFAULT_REF


def cache_dir() -> Path:
    return Path(tempfile.gettempdir()) / "ccp4i2-monomer-library" / library_ref()


def find_monomer_cif(code, opener=None) -> Path:
    """Path to the entry's CIF, fetching it if no local library exists.

    Raises MonomerLibraryError: 400 for a malformed code, 404 when the
    library has no such entry, 502 when the published library is unreachable.
    """
    code = normalise_code(code)
    rel = relative_path(code)

    local = local_library_dir()
    if local is not None:
        path = local / rel
        if not path.is_file():
            raise MonomerLibraryError(404, f"Monomer '{code}' not found in library")
        return path

    cached = cache_dir() / rel
    if cached.is_file():
        return cached

    url = f"{REMOTE_BASE}/{library_ref()}/{rel}"
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    opener = opener or urllib.request.urlopen
    try:
        with opener(request, timeout=TIMEOUT) as response:
            content = response.read()
    except urllib.error.HTTPError as err:
        if err.code == 404:
            raise MonomerLibraryError(404, f"Monomer '{code}' not found in library")
        raise MonomerLibraryError(502, f"Monomer library answered {err.code} for {code}")
    except (urllib.error.URLError, OSError) as err:
        raise MonomerLibraryError(502, f"Could not reach the monomer library: {err}")

    # Write then rename, so a concurrent reader never sees half a file.
    cached.parent.mkdir(parents=True, exist_ok=True)
    partial = cached.with_name(f"{cached.name}.{os.getpid()}.part")
    partial.write_bytes(content)
    os.replace(partial, cached)
    logger.info("Fetched monomer %s from %s", code, url)
    return cached
