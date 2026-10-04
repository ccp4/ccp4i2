"""The MCP SDK's request-state codec, with a key derivation CCP4's cryptography can run.

The SDK seals multi-round request state with AES-256-GCM under a key derived
by HKDF-SHA256. Its derivation calls cryptography's HKDF without a
``backend``, which the cryptography CCP4 ships (2.8, compiled, and not ours to
replace: requirements-runtime.txt excludes it) refuses, so constructing the
server failed in ccp4-python. HKDF is a few lines of HMAC; this derives the
same key with the standard library, and keeps the SDK's token format and
AES-GCM (in cryptography since 2.0). Passed to MCPServer through the SDK's own
hook, ``RequestStateSecurity(codec=...)``.
"""
import base64
import hashlib
import hmac
import os

from cryptography.exceptions import InvalidTag
from cryptography.hazmat.primitives.ciphers.aead import AESGCM
from mcp.server.request_state import InvalidRequestState

# As the SDK's (mcp/server/request_state.py), so tokens are the same format
_KDF_INFO = b"mcp/request-state/v1/aes-256-gcm"
_KID_INFO = b"mcp/request-state/v1/kid:"
_TOKEN_PREFIX = "v1."
_KID_LEN = 4
_NONCE_LEN = 12


def hkdf_sha256(secret: bytes, info: bytes, length: int = 32) -> bytes:
    """HKDF-SHA256 (RFC 5869) with no salt, as cryptography's HKDF(salt=None)."""
    prk = hmac.new(b"\0" * hashlib.sha256().digest_size, secret, hashlib.sha256).digest()
    okm, block, counter = b"", b"", 1
    while len(okm) < length:
        block = hmac.new(prk, block + info + bytes([counter]), hashlib.sha256).digest()
        okm += block
        counter += 1
    return okm[:length]


def _b64u(data: bytes) -> str:
    return base64.urlsafe_b64encode(data).decode().rstrip("=")


def _b64u_decode(text: str) -> bytes:
    raw = base64.urlsafe_b64decode(text + "=" * (-len(text) % 4))
    if _b64u(raw) != text:
        raise ValueError("non-canonical base64url")
    return raw


class AESGCMCodec:
    """AES-256-GCM request-state tokens under one ephemeral key (per process)."""

    def __init__(self, secret: bytes = None):
        key = hkdf_sha256(secret if secret is not None else os.urandom(32), _KDF_INFO)
        self._kid = hashlib.sha256(_KID_INFO + key).digest()[:_KID_LEN]
        self._aead = AESGCM(key)

    def seal(self, payload: bytes) -> str:
        nonce = os.urandom(_NONCE_LEN)
        sealed = self._aead.encrypt(nonce, payload, _TOKEN_PREFIX.encode() + self._kid)
        return _TOKEN_PREFIX + _b64u(self._kid + nonce + sealed)

    def unseal(self, token: str) -> bytes:
        if not token.startswith(_TOKEN_PREFIX):
            raise InvalidRequestState("malformed")
        try:
            raw = _b64u_decode(token[len(_TOKEN_PREFIX):])
        except ValueError as exc:
            raise InvalidRequestState("malformed") from exc
        if len(raw) < _KID_LEN + _NONCE_LEN + 16:
            raise InvalidRequestState("malformed")
        kid, nonce, sealed = raw[:_KID_LEN], raw[_KID_LEN:_KID_LEN + _NONCE_LEN], raw[_KID_LEN + _NONCE_LEN:]
        if kid != self._kid:
            raise InvalidRequestState("unknown key")
        try:
            return self._aead.decrypt(nonce, sealed, _TOKEN_PREFIX.encode() + kid)
        except InvalidTag:
            raise InvalidRequestState("seal") from None


def request_state_security():
    from mcp.server.request_state import RequestStateSecurity
    return RequestStateSecurity(codec=AESGCMCodec())
