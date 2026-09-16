"""Autopopulation gives sibling slots distinct files.

The two half-map inputs used to both grab file_id_list[0] -- the same file --
because each slot was populated independently. A shared "claimed" set makes each
slot prefer a candidate no sibling has taken. Safe across all slots at once: by
the mini-MTZ data model two slots only ever share a candidate pool when they
have the same (mimeType, subType, contentFlag).
"""

import types

import ccp4i2.lib.utils.parameters.set_input_by_context as mod


class _Dobj:
    def __init__(self, name="MAPIN", subtype=5):
        self._name = name
        self._subtype = subtype

    def qualifiers(self, key):
        return {
            "requiredSubType": self._subtype,
            "requiredContentFlag": None,
            "mimeTypeName": "application/CCP4-map",
        }.get(key)

    def objectPath(self):
        return f"x.{self._name}"


def _job():
    return types.SimpleNamespace(project=types.SimpleNamespace(uuid="proj-uuid"))


def _patch(monkeypatch, candidates):
    recorded = []
    monkeypatch.setattr(mod, "get_file_by_job_context", lambda **kw: list(candidates))
    monkeypatch.setattr(mod, "_set_file_from_db",
                        lambda dobj, fid, job: recorded.append((dobj._name, fid)))
    return recorded


def test_two_slots_get_distinct_files(monkeypatch):
    recorded = _patch(monkeypatch, ["halfmap-A", "halfmap-B"])
    claimed = set()
    mod._populate_file_from_context(_Dobj("MAPIN1"), "ctx", _job(), claimed)
    mod._populate_file_from_context(_Dobj("MAPIN2"), "ctx", _job(), claimed)
    assert [fid for _, fid in recorded] == ["halfmap-A", "halfmap-B"]


def test_falls_back_to_first_when_only_one_candidate(monkeypatch):
    # One file, two slots: never populate less than before -> both get it.
    recorded = _patch(monkeypatch, ["only-map"])
    claimed = set()
    mod._populate_file_from_context(_Dobj("MAPIN1"), "ctx", _job(), claimed)
    mod._populate_file_from_context(_Dobj("MAPIN2"), "ctx", _job(), claimed)
    assert [fid for _, fid in recorded] == ["only-map", "only-map"]


def test_without_claimed_set_unchanged(monkeypatch):
    # Backwards-compatible: no claimed set -> first candidate, as before.
    recorded = _patch(monkeypatch, ["a", "b"])
    mod._populate_file_from_context(_Dobj(), "ctx", _job(), None)
    assert recorded == [("MAPIN", "a")]
