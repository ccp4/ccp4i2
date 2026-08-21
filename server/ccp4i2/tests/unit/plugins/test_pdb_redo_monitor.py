"""Unit tests for the PDB-REDO polling loop (no network).

The behaviour under test is the difference between the three ways waiting can
go wrong, because conflating them is what made the original bare ``except:``
unhelpful:

* the remote run failed        -> PDBRedoJobStopped (salvage the logs)
* our credentials stopped working -> PDBRedoAuthError (user must act)
* we lost contact              -> PDBRedoUnreachable (remote run is probably fine)

Plus the property that matters most in practice: a transient blip during a
multi-hour run must NOT fail the job.
"""

import pytest

from ccp4i2.wrappers.pdb_redo_api.script import test_api


class FakeResponse:
    def __init__(self, status_code=200, payload=None):
        self.status_code = status_code
        self._payload = payload

    def json(self):
        if self._payload is None:
            raise ValueError("no json")
        return self._payload

    def raise_for_status(self):
        if self.status_code >= 400:
            import requests

            raise requests.HTTPError(f"HTTP {self.status_code}")


@pytest.fixture
def no_sleep(monkeypatch):
    """Make the loop run instantly, and record the delays it asked for."""
    slept = []
    monkeypatch.setattr(test_api.time, "sleep", lambda d: slept.append(d))
    return slept


def _responses(monkeypatch, sequence):
    """Drive requests.get from a list of responses/exceptions."""
    items = list(sequence)
    calls = {"n": 0}

    def fake_get(*args, **kwargs):
        calls["n"] += 1
        item = items.pop(0) if items else items_last[0]
        items_last[0] = item
        if isinstance(item, Exception):
            raise item
        return item

    items_last = [items[-1] if items else FakeResponse(200, {"status": "ended"})]
    monkeypatch.setattr(test_api.requests, "get", fake_get)
    return calls


def test_returns_when_the_run_ends(monkeypatch, no_sleep):
    _responses(
        monkeypatch,
        [
            FakeResponse(200, {"status": "starting"}),
            FakeResponse(200, {"status": "running"}),
            FakeResponse(200, {"status": "ended"}),
        ],
    )
    assert test_api.monitor(42, "id", "secret") == "ended"


def test_remote_failure_is_distinguishable(monkeypatch, no_sleep):
    _responses(monkeypatch, [FakeResponse(200, {"status": "stopped"})])
    with pytest.raises(test_api.PDBRedoJobStopped) as excinfo:
        test_api.monitor(42, "id", "secret")
    # The run id must travel with the error: it is the one thing the user needs.
    assert excinfo.value.run_id == 42


def test_revoked_token_fails_immediately(monkeypatch, no_sleep):
    """A 401 will not fix itself, so it must not be retried for ten minutes."""
    calls = _responses(
        monkeypatch,
        [FakeResponse(401), FakeResponse(200, {"status": "ended"})],
    )
    with pytest.raises(test_api.PDBRedoAuthError):
        test_api.monitor(42, "id", "secret")
    assert calls["n"] == 1


def test_transient_failure_does_not_fail_the_job(monkeypatch, no_sleep):
    """The whole point: a blip mid-run must be ridden out, not fatal."""
    import requests

    _responses(
        monkeypatch,
        [
            FakeResponse(200, {"status": "running"}),
            requests.ConnectionError("dropped"),
            FakeResponse(502),
            FakeResponse(200, None),  # malformed body
            FakeResponse(200, {"status": "ended"}),
        ],
    )
    assert test_api.monitor(42, "id", "secret") == "ended"


def test_gives_up_after_enough_consecutive_failures(monkeypatch, no_sleep):
    import requests

    _responses(monkeypatch, [requests.ConnectionError("down")])
    with pytest.raises(test_api.PDBRedoUnreachable) as excinfo:
        test_api.monitor(42, "id", "secret", max_consecutive_failures=4)
    assert excinfo.value.run_id == 42
    # "may still be running" is the actionable half of the message.
    assert "still be running" in str(excinfo.value)


def test_a_success_resets_the_failure_count(monkeypatch, no_sleep):
    """Failures must be *consecutive*; a run that blips repeatedly over hours
    but keeps recovering should still be waited on."""
    import requests

    err = requests.ConnectionError("blip")
    _responses(
        monkeypatch,
        [
            err, err, err,
            FakeResponse(200, {"status": "running"}),
            err, err, err,
            FakeResponse(200, {"status": "ended"}),
        ],
    )
    assert test_api.monitor(42, "id", "secret", max_consecutive_failures=4) == "ended"


def test_polling_backs_off_and_is_capped(monkeypatch, no_sleep):
    """5s flat for a multi-hour run is ~720 requests/hour against pdb-redo.eu."""
    _responses(
        monkeypatch,
        [FakeResponse(200, {"status": "running"})] * 30
        + [FakeResponse(200, {"status": "ended"})],
    )
    test_api.monitor(42, "id", "secret", poll_initial=5, poll_max=60)
    assert no_sleep[0] == 5
    assert max(no_sleep) <= 60
    assert no_sleep[-1] == 60  # reached the cap
    # Over the same wall-clock, far fewer requests than a flat 5s poll.
    assert sum(no_sleep) / len(no_sleep) > 20


def test_status_change_restores_responsiveness(monkeypatch, no_sleep):
    """After backing off, a transition should shorten the interval again."""
    _responses(
        monkeypatch,
        [FakeResponse(200, {"status": "running"})] * 10
        + [FakeResponse(200, {"status": "finishing"})]
        + [FakeResponse(200, {"status": "ended"})],
    )
    test_api.monitor(42, "id", "secret", poll_initial=5, poll_max=60)
    assert no_sleep[-1] == 5
