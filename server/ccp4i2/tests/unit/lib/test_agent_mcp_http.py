"""The MCP facade served by the app at /mcp/ccp4i2 (ccp4i2/agent/http.py), and the
request-state codec that lets the MCP SDK run on CCP4's cryptography."""
import asyncio
import os

import pytest

pytest.importorskip("mcp", reason="needs the mcp package")

from ccp4i2.agent import request_state  # noqa: E402
from ccp4i2.agent import http, mcp_server  # noqa: E402
from ccp4i2.agent.http import MCP_PATH, mcp_path, with_mcp  # noqa: E402


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


def _http(path, authorization=None, forwarded=None):
    headers = [(b"authorization", authorization.encode())] if authorization else []
    if forwarded:
        headers.append((b"x-forwarded-for", forwarded.encode()))
    return {"type": "http", "path": path, "method": "POST", "headers": headers,
            "server": ("127.0.0.1", 8000), "client": ("10.1.2.3", 51515)}


@pytest.fixture
def deployment(monkeypatch):
    """A deployment that asked for MCP: no desktop token, CCP4I2_MCP=1."""
    monkeypatch.delenv("CCP4I2_LOCAL_SESSION_TOKEN", raising=False)
    monkeypatch.setenv("CCP4I2_MCP", "1")


def test_the_path_is_scoped():
    # Beside /api/ccp4i2, so an app serving CCP4i2 can serve its own MCP too
    assert MCP_PATH == "/mcp/ccp4i2"


def test_everything_but_mcp_goes_to_django(deployment):
    django = _Recorder()
    app = with_mcp(django)
    sent = _call(app, _http("/api/ccp4i2/projects/"))
    assert sent[-1]["body"] == b"django" and len(django.scopes) == 1
    _call(app, _http("/mcp/ccp4i2anything"))  # not the MCP path
    _call(app, _http("/mcp", "Bearer x"))  # the old, unscoped path
    assert len(django.scopes) == 3


def test_a_deployment_serves_mcp_only_when_it_asks(monkeypatch):
    django = _Recorder()
    monkeypatch.delenv("CCP4I2_LOCAL_SESSION_TOKEN", raising=False)
    monkeypatch.delenv("CCP4I2_MCP", raising=False)
    assert with_mcp(django) is django  # updating CCP4i2 adds no endpoint
    monkeypatch.setenv("CCP4I2_MCP", "1")
    assert with_mcp(django) is not django


def test_mcp_can_be_turned_off_on_the_desktop_too(monkeypatch):
    django = _Recorder()
    monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", "secret-token")
    assert with_mcp(django) is not django  # on by default there
    monkeypatch.setenv("CCP4I2_MCP", "0")
    assert with_mcp(django) is django


def test_a_deployment_refuses_a_request_without_credentials(deployment):
    django = _Recorder()
    sent = _call(with_mcp(django), _http(MCP_PATH))
    assert sent[0]["status"] == 401 and not django.scopes


def test_tool_calls_go_as_the_caller_with_their_address(deployment, monkeypatch):
    seen = {}

    async def fake_mcp_app(scope, receive, send):
        seen.update(mcp_server._caller.get())
        await send({"type": "http.response.start", "status": 200, "headers": []})
        await send({"type": "http.response.body", "body": b"mcp"})

    monkeypatch.setattr(mcp_server.server, "streamable_http_app", lambda **kw: fake_mcp_app)
    _call(with_mcp(_Recorder()), _http(MCP_PATH, "Bearer user-token", forwarded="203.0.113.7"))
    assert seen["authorization"] == "Bearer user-token"
    assert seen["forwarded_for"] == "203.0.113.7, 10.1.2.3"  # the chain, then the peer
    assert seen["base"] == "http://127.0.0.1:8000/api/ccp4i2"


def test_an_http_caller_never_borrows_the_server_token(monkeypatch):
    monkeypatch.setenv("CCP4I2_TOKEN", "server-token")
    token = mcp_server._caller.set({"authorization": None, "forwarded_for": "10.1.2.3"})
    try:
        assert mcp_server._auth_headers() == {"X-Forwarded-For": "10.1.2.3"}
    finally:
        mcp_server._caller.reset(token)
    assert mcp_server._auth_headers() == {"Authorization": "Bearer server-token"}  # stdio


def test_on_the_desktop_mcp_needs_the_session_token(monkeypatch):
    monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", "secret-token")
    django = _Recorder()
    app = with_mcp(django)
    for auth in (None, "Bearer wrong", "secret-token"):
        sent = _call(app, _http(MCP_PATH, auth))
        assert sent[0]["status"] == 401
    assert not django.scopes  # never reached Django either


@pytest.mark.parametrize("setting, expected", [
    (None, "/mcp/ccp4i2"),
    ("", "/mcp/ccp4i2"),
    ("/", "/mcp/ccp4i2"),
    ("  ", "/mcp/ccp4i2"),
    ("/tools/ccp4i2/", "/tools/ccp4i2"),
    ("tools/ccp4i2", "/tools/ccp4i2"),
])
def test_a_host_can_move_the_mount_but_never_onto_its_root(setting, expected):
    assert mcp_path(setting) == expected
