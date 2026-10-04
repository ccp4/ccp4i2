"""The MCP facade served by the app at /mcp (ccp4i2/agent/http.py), and the
request-state codec that lets the MCP SDK run on CCP4's cryptography."""
import asyncio
import os

import pytest

pytest.importorskip("mcp", reason="needs the mcp package")

from ccp4i2.agent import request_state  # noqa: E402
from ccp4i2.agent.http import MCP_PATH, with_mcp  # noqa: E402


# --- the codec --------------------------------------------------------------

def test_key_derivation_is_cryptographys_hkdf():
    # The point of the module: the same key as cryptography's HKDF, computed
    # without the call CCP4's cryptography 2.8 cannot make.
    from cryptography.hazmat.primitives.hashes import SHA256
    from cryptography.hazmat.primitives.kdf.hkdf import HKDF
    try:
        modern = HKDF(algorithm=SHA256(), length=32, salt=None, info=b"x").derive
    except TypeError:
        pytest.skip("this cryptography needs a backend; the comparison needs a modern one")
    for secret in (b"\0" * 32, os.urandom(32), os.urandom(64)):
        expected = HKDF(algorithm=SHA256(), length=32, salt=None,
                        info=request_state._KDF_INFO).derive(secret)
        assert request_state.hkdf_sha256(secret, request_state._KDF_INFO) == expected
    assert modern  # (built once above only to detect the API)


def test_tokens_round_trip_and_tampering_fails():
    from mcp.server.request_state import InvalidRequestState
    codec = request_state.AESGCMCodec()
    token = codec.seal(b'{"state": 1}')
    assert token.startswith("v1.") and codec.unseal(token) == b'{"state": 1}'
    tampered = token[:-2] + ("A" if token[-2] != "A" else "B") + token[-1]
    with pytest.raises(InvalidRequestState):
        codec.unseal(tampered)
    with pytest.raises(InvalidRequestState):
        request_state.AESGCMCodec().unseal(token)  # another process's key
    with pytest.raises(InvalidRequestState):
        codec.unseal("v2.abc")


# --- the mount ----------------------------------------------------------------

class _Recorder:
    def __init__(self):
        self.scopes = []

    async def __call__(self, scope, receive, send):
        self.scopes.append(scope)
        if scope["type"] == "http":
            await send({"type": "http.response.start", "status": 200, "headers": []})
            await send({"type": "http.response.body", "body": b"django"})


def _call(app, scope):
    sent = []

    async def receive():
        return {"type": "http.request", "body": b"", "more_body": False}

    async def send(message):
        sent.append(message)

    asyncio.run(app(scope, receive, send))
    return sent


def _http(path, authorization=None):
    headers = [(b"authorization", authorization.encode())] if authorization else []
    return {"type": "http", "path": path, "method": "POST", "headers": headers,
            "server": ("127.0.0.1", 8000)}


def test_everything_but_mcp_goes_to_django():
    django = _Recorder()
    app = with_mcp(django)
    sent = _call(app, _http("/api/ccp4i2/projects/"))
    assert sent[-1]["body"] == b"django" and len(django.scopes) == 1
    _call(app, _http("/mcpanything"))  # not the MCP path
    assert len(django.scopes) == 2


def test_on_the_desktop_mcp_needs_the_session_token(monkeypatch):
    monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", "secret-token")
    django = _Recorder()
    app = with_mcp(django)
    for auth in (None, "Bearer wrong", "secret-token"):
        sent = _call(app, _http(MCP_PATH, auth))
        assert sent[0]["status"] == 401
    assert not django.scopes  # never reached Django either
