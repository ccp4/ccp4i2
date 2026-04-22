"""
Unit tests for the HGNC REST client.

Uses a monkey-patched _fetch_json so no network calls are made.
"""

from __future__ import annotations

import json
from io import BytesIO
from unittest import mock
from urllib.error import HTTPError, URLError

import pytest

from compounds.registry import hgnc


# ----- fixtures / helpers -------------------------------------------------- #

EGFR_PAYLOAD = {
    "responseHeader": {"status": 0, "QTime": 1},
    "response": {
        "numFound": 1,
        "start": 0,
        "docs": [
            {
                "symbol": "EGFR",
                "hgnc_id": "HGNC:3236",
                "name": "epidermal growth factor receptor",
                "alias_symbol": ["ERBB", "ERBB1"],
                "prev_symbol": ["HER1"],
                "uniprot_ids": ["P00533"],
                "ensembl_gene_id": "ENSG00000146648",
            }
        ],
    },
}

EMPTY_PAYLOAD = {
    "responseHeader": {"status": 0},
    "response": {"numFound": 0, "start": 0, "docs": []},
}


def _patched_urlopen(payload: dict):
    """Return a mock of urllib.request.urlopen that yields `payload`."""
    body = json.dumps(payload).encode("utf-8")

    class _Resp:
        def __init__(self, data):
            self._data = data

        def read(self):
            return self._data

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    def _urlopen(req, timeout=None):
        return _Resp(body)

    return _urlopen


# ----- fetch_by_symbol ----------------------------------------------------- #


def test_fetch_by_symbol_returns_parsed_record():
    with mock.patch.object(hgnc, "urlopen", _patched_urlopen(EGFR_PAYLOAD)):
        rec = hgnc.fetch_by_symbol("EGFR")

    assert rec is not None
    assert rec.symbol == "EGFR"
    assert rec.hgnc_id == "HGNC:3236"
    assert rec.name == "epidermal growth factor receptor"
    assert rec.aliases == ["ERBB", "ERBB1", "HER1"]  # union, order preserved
    assert rec.uniprot_ids == ["P00533"]
    assert rec.ensembl_gene_id == "ENSG00000146648"


def test_fetch_by_symbol_returns_none_on_empty():
    with mock.patch.object(hgnc, "urlopen", _patched_urlopen(EMPTY_PAYLOAD)):
        assert hgnc.fetch_by_symbol("NOPE") is None


def test_fetch_by_symbol_normalizes_case():
    """Lowercased input should be uppercased before the HGNC call."""
    calls = []

    def _urlopen(req, timeout=None):
        calls.append(req.full_url)
        body = json.dumps(EMPTY_PAYLOAD).encode("utf-8")

        class _R:
            def read(self):
                return body

            def __enter__(self):
                return self

            def __exit__(self, *a):
                return False

        return _R()

    with mock.patch.object(hgnc, "urlopen", _urlopen):
        hgnc.fetch_by_symbol(" egfr ")

    assert len(calls) == 1
    assert calls[0].endswith("/fetch/symbol/EGFR")


def test_fetch_by_symbol_empty_input_returns_none_without_call():
    called = []

    def _urlopen(req, timeout=None):
        called.append(req)
        raise AssertionError("urlopen should not be called for empty input")

    with mock.patch.object(hgnc, "urlopen", _urlopen):
        assert hgnc.fetch_by_symbol("") is None
        assert hgnc.fetch_by_symbol("   ") is None

    assert called == []


# ----- search --------------------------------------------------------------- #


def test_search_returns_list_of_records():
    payload = {
        "response": {
            "numFound": 2,
            "docs": [
                {"symbol": "EGFR", "hgnc_id": "HGNC:3236"},
                {"symbol": "AREG", "hgnc_id": "HGNC:651"},
            ],
        }
    }
    with mock.patch.object(hgnc, "urlopen", _patched_urlopen(payload)):
        recs = hgnc.search("epidermal")

    assert [r.symbol for r in recs] == ["EGFR", "AREG"]
    assert recs[0].hgnc_id == "HGNC:3236"


def test_search_empty_query_returns_empty_list():
    assert hgnc.search("") == []
    assert hgnc.search("   ") == []


# ----- error surface -------------------------------------------------------- #


def test_http_error_wraps_as_hgnc_error():
    def _urlopen(req, timeout=None):
        raise HTTPError(req.full_url, 503, "Service Unavailable", {}, BytesIO(b""))

    with mock.patch.object(hgnc, "urlopen", _urlopen):
        with pytest.raises(hgnc.HGNCError) as exc:
            hgnc.fetch_by_symbol("EGFR")

    assert "HTTP 503" in str(exc.value)


def test_url_error_wraps_as_hgnc_error():
    def _urlopen(req, timeout=None):
        raise URLError("name resolution failed")

    with mock.patch.object(hgnc, "urlopen", _urlopen):
        with pytest.raises(hgnc.HGNCError) as exc:
            hgnc.fetch_by_symbol("EGFR")

    assert "network error" in str(exc.value)


def test_non_json_body_wraps_as_hgnc_error():
    class _R:
        def read(self):
            return b"<html>nope</html>"

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    with mock.patch.object(hgnc, "urlopen", lambda req, timeout=None: _R()):
        with pytest.raises(hgnc.HGNCError) as exc:
            hgnc.fetch_by_symbol("EGFR")

    assert "non-JSON" in str(exc.value)


# ----- parser edge cases --------------------------------------------------- #


def test_parse_doc_handles_missing_alias_fields():
    payload = {
        "response": {
            "numFound": 1,
            "docs": [
                {"symbol": "BARE", "hgnc_id": "HGNC:1"},  # no aliases / uniprot
            ],
        }
    }
    with mock.patch.object(hgnc, "urlopen", _patched_urlopen(payload)):
        rec = hgnc.fetch_by_symbol("BARE")

    assert rec.aliases == []
    assert rec.uniprot_ids == []
    assert rec.ensembl_gene_id == ""


def test_parse_doc_dedupes_aliases():
    payload = {
        "response": {
            "numFound": 1,
            "docs": [
                {
                    "symbol": "DUP",
                    "alias_symbol": ["A", "B", "A"],
                    "prev_symbol": ["B", "C"],
                }
            ],
        }
    }
    with mock.patch.object(hgnc, "urlopen", _patched_urlopen(payload)):
        rec = hgnc.fetch_by_symbol("DUP")

    assert rec.aliases == ["A", "B", "C"]  # deduped, order preserved
