"""
Finding a monomer library entry: locally, or from the published library.

Written after the Make Link interface on a served deployment said every
library residue ("CYS", "LYS", ...) was "not in the CCP4 monomer library".
The API server there carries no CCP4, so the endpoint, which only read
$CCP4/lib/data/monomers, had nothing to read. It now falls back to the
published library (github.com/MonomerLibrary/monomers) at a pinned tag.

Pure Python -- no network (the opener is faked), no CCP4.
"""
import io
import urllib.error

import pytest

from ccp4i2.lib.utils.files import monomer_library as ml
from ccp4i2.lib.utils.files.monomer_library import MonomerLibraryError, find_monomer_cif


@pytest.fixture(autouse=True)
def no_local_library(tmp_path, monkeypatch):
    """No CCP4, no CLIBD_MON, and a private cache: every test starts remote."""
    monkeypatch.delenv("CCP4", raising=False)
    monkeypatch.delenv("CLIBD_MON", raising=False)
    monkeypatch.delenv("CCP4I2_MONOMER_LIBRARY_REF", raising=False)
    monkeypatch.setattr(ml.tempfile, "gettempdir", lambda: str(tmp_path / "tmp"))


class FakeOpener:
    """Records the URLs asked for; answers from a dict of URL -> bytes."""

    def __init__(self, files=None, error=None):
        self.files = files or {}
        self.error = error
        self.urls = []

    def __call__(self, request, timeout=None):
        self.urls.append(request.full_url)
        if self.error is not None:
            raise self.error
        if request.full_url not in self.files:
            raise urllib.error.HTTPError(request.full_url, 404, "Not Found", {}, None)
        return io.BytesIO(self.files[request.full_url])


def remote_url(rel, ref=ml.DEFAULT_REF):
    return f"{ml.REMOTE_BASE}/{ref}/{rel}"


@pytest.mark.parametrize("code, rel", [
    ("CYS", "c/CYS.cif"),
    ("cys", "c/CYS.cif"),
    ("0G6", "0/0G6.cif"),
    # Reserved file names on Windows are stored doubled.
    ("CON", "c/CON_CON.cif"),
    ("NUL", "n/NUL_NUL.cif"),
])
def test_relative_path(code, rel):
    assert ml.relative_path(ml.normalise_code(code)) == rel


@pytest.mark.parametrize("code", ["", "   ", "../etc", "C/S", "CYS.cif", "A" * 9])
def test_malformed_code_is_400_and_never_fetched(code):
    opener = FakeOpener()
    with pytest.raises(MonomerLibraryError) as err:
        find_monomer_cif(code, opener=opener)
    assert err.value.status == 400
    assert opener.urls == []


def test_local_ccp4_library_wins(tmp_path, monkeypatch):
    lib = tmp_path / "ccp4" / "lib" / "data" / "monomers"
    (lib / "c").mkdir(parents=True)
    (lib / "c" / "CYS.cif").write_text("local")
    monkeypatch.setenv("CCP4", str(tmp_path / "ccp4"))
    opener = FakeOpener()

    assert find_monomer_cif("CYS", opener=opener) == lib / "c" / "CYS.cif"
    assert opener.urls == []


def test_clibd_mon_wins_over_ccp4(tmp_path, monkeypatch):
    mon = tmp_path / "mon"
    (mon / "c").mkdir(parents=True)
    (mon / "c" / "CYS.cif").write_text("clibd")
    monkeypatch.setenv("CLIBD_MON", str(mon))
    monkeypatch.setenv("CCP4", str(tmp_path / "elsewhere"))

    assert find_monomer_cif("CYS").read_text() == "clibd"


def test_missing_from_local_library_is_404_without_fetching(tmp_path, monkeypatch):
    """A local library is authoritative: no network call behind its back."""
    (tmp_path / "mon").mkdir()
    monkeypatch.setenv("CLIBD_MON", str(tmp_path / "mon"))
    opener = FakeOpener()
    with pytest.raises(MonomerLibraryError) as err:
        find_monomer_cif("ZZZ", opener=opener)
    assert err.value.status == 404
    assert opener.urls == []


def test_no_local_library_fetches_once_then_caches():
    opener = FakeOpener({remote_url("c/CYS.cif"): b"remote cys"})

    first = find_monomer_cif("CYS", opener=opener)
    second = find_monomer_cif("cys", opener=opener)

    assert first == second
    assert first.read_bytes() == b"remote cys"
    assert opener.urls == [remote_url("c/CYS.cif")]
    assert not list(first.parent.glob("*.part"))


def test_ref_is_configurable(monkeypatch):
    monkeypatch.setenv("CCP4I2_MONOMER_LIBRARY_REF", "ccp4-9.0.012")
    opener = FakeOpener({remote_url("c/CYS.cif", "ccp4-9.0.012"): b"older"})
    assert find_monomer_cif("CYS", opener=opener).read_bytes() == b"older"


def test_remote_404_is_not_found():
    with pytest.raises(MonomerLibraryError) as err:
        find_monomer_cif("ZZZ", opener=FakeOpener())
    assert err.value.status == 404


@pytest.mark.parametrize("error", [
    urllib.error.URLError("no route to host"),
    TimeoutError("timed out"),
    urllib.error.HTTPError("u", 503, "Unavailable", {}, None),
])
def test_unreachable_library_is_502_not_404(error):
    """Unreachable must not read as "no such monomer"."""
    with pytest.raises(MonomerLibraryError) as err:
        find_monomer_cif("CYS", opener=FakeOpener(error=error))
    assert err.value.status == 502
