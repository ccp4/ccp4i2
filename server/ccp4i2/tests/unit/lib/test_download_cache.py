"""
The i2run suite's download cache.

Written after Zenodo was unreachable on 2026-08-26 and two xia2 tests errored
in *setup* with a read timeout. In a baseline that reads as two failures, and
distinguishing "the code broke" from "a server was down" cost a round of
curl-ing the host by hand.

Caching also stops the suite re-fetching the same fixtures on every run. Files
above a size cap are still fetched fresh, because the xia2 archive is 335 MB
and uses twenty frames out of it.

Pure Python -- no network, no CCP4.
"""
import io
from pathlib import Path

import pytest

from ccp4i2.tests.i2run import utils


@pytest.fixture(name="cache")
def cache_fixture(tmp_path, monkeypatch):
    """Point the cache at a temporary directory."""
    monkeypatch.setattr(utils, "_download_cache_dir", lambda: tmp_path)
    monkeypatch.delenv("CCP4I2_TEST_REFETCH", raising=False)
    monkeypatch.delenv("CCP4I2_TEST_CACHE_MAX_MB", raising=False)
    return tmp_path


class _Headers:
    def __init__(self, filename=None):
        self._filename = filename

    def get_filename(self):
        return self._filename


class _Response(io.BytesIO):
    """Just enough of urlopen()'s result: a body and a headers object."""
    def __init__(self, payload, filename=None):
        super().__init__(payload)
        self.headers = _Headers(filename)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


def _serve(monkeypatch, payload, counter):
    def fake_urlopen(url, timeout=None):
        counter.append(url)
        return _Response(payload)
    monkeypatch.setattr(utils, "urlopen", fake_urlopen)


def test_a_file_is_fetched_once_and_reused(cache, monkeypatch):
    fetches = []
    _serve(monkeypatch, b"reflection data", fetches)

    with utils.download("https://example.org/db/1abc/1abc.mtz") as first:
        assert Path(first).read_bytes() == b"reflection data"
    with utils.download("https://example.org/db/1abc/1abc.mtz") as second:
        assert Path(second).read_bytes() == b"reflection data"

    assert len(fetches) == 1, "the second call should not have gone to the network"
    assert first == second


def test_the_cached_file_survives_the_context_manager(cache, monkeypatch):
    """The old helper deleted its temp file on exit; a cache that did the same
    would be no cache at all."""
    _serve(monkeypatch, b"payload", [])
    with utils.download("https://example.org/a.cif") as path:
        pass
    assert Path(path).is_file()


def test_two_urls_with_the_same_basename_do_not_collide(cache, monkeypatch):
    _serve(monkeypatch, b"one", [])
    with utils.download("https://a.example/db/1abc/final.mtz") as first:
        Path(first).write_bytes(b"one")
    _serve(monkeypatch, b"two", [])
    with utils.download("https://b.example/db/9xyz/final.mtz") as second:
        assert first != second
        assert Path(second).read_bytes() == b"two"
    assert Path(first).read_bytes() == b"one"


def test_a_large_file_is_not_cached(cache, monkeypatch):
    """335 MB to use twenty frames is not a trade worth making by default."""
    monkeypatch.setenv("CCP4I2_TEST_CACHE_MAX_MB", "0.001")   # 1 KB
    fetches = []
    _serve(monkeypatch, b"x" * 5000, fetches)

    with utils.download("https://example.org/big.tar.bz2") as path:
        assert Path(path).is_file()
    assert not Path(path).exists(), "a file too big to keep must be cleaned up"
    assert list(cache.iterdir()) == []

    with utils.download("https://example.org/big.tar.bz2"):
        pass
    assert len(fetches) == 2, "an uncached file is fetched again"


def test_no_limit_when_the_cap_is_zero(cache, monkeypatch):
    monkeypatch.setenv("CCP4I2_TEST_CACHE_MAX_MB", "0")
    _serve(monkeypatch, b"y" * 5000, [])
    with utils.download("https://example.org/big.tar.bz2") as path:
        pass
    assert Path(path).is_file()


def test_refetch_ignores_what_is_cached(cache, monkeypatch):
    fetches = []
    _serve(monkeypatch, b"first", fetches)
    with utils.download("https://example.org/a.cif"):
        pass

    monkeypatch.setenv("CCP4I2_TEST_REFETCH", "1")
    _serve(monkeypatch, b"second", fetches)
    with utils.download("https://example.org/a.cif") as path:
        assert Path(path).read_bytes() == b"second"
    assert len(fetches) == 2


def test_a_truncated_cache_entry_is_not_reused(cache, monkeypatch):
    """An interrupted run must not leave half a file to be trusted later."""
    fetches = []
    _serve(monkeypatch, b"whole file", fetches)
    with utils.download("https://example.org/a.cif") as path:
        pass
    Path(path).write_bytes(b"")          # as an interrupted download would leave it

    with utils.download("https://example.org/a.cif") as again:
        assert Path(again).read_bytes() == b"whole file"
    assert len(fetches) == 2


def test_the_servers_name_is_what_the_cache_keeps(cache, monkeypatch):
    """Robetta serves a PDB from a URL ending ".php".

    Naming the cached copy after the URL cost it the extension, and every
    reader in CCP4i2 tells a format by extension: the first cached run failed
    with "Unknown format of ...models_download.php" on a file that was
    perfectly good. The readable half of the name is not knowable before the
    request, so lookup is by key and the name comes from the response.
    """
    def fake_urlopen(url, timeout=None):
        return _Response(b"ATOM      1  N   MET A   1\n", filename="model.pdb")
    monkeypatch.setattr(utils, "urlopen", fake_urlopen)

    url = "https://robetta.bakerlab.org/models_download.php?id=14697"
    with utils.download(url) as path:
        first = Path(path)
    assert first.name.endswith("_model.pdb"), first.name

    # And the copy is found again, despite the name being unguessable from
    # the URL alone.
    def refuse(url, timeout=None):
        raise AssertionError("should have been served from the cache")
    monkeypatch.setattr(utils, "urlopen", refuse)
    with utils.download(url) as path:
        assert Path(path) == first


def test_a_url_with_no_useful_name_still_caches(cache, monkeypatch):
    def fake_urlopen(url, timeout=None):
        return _Response(b"data")
    monkeypatch.setattr(utils, "urlopen", fake_urlopen)

    with utils.download("https://example.org/") as path:
        assert Path(path).is_file()
