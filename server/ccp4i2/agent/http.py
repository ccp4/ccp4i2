"""The MCP facade served by the app itself, at /mcp/ccp4i2 beside the REST API.

An agent connects with a URL and the same Bearer token as the REST API
(the desktop app shows both in Help > About; they change at every launch):

    claude mcp add --transport http ccp4i2 http://127.0.0.1:<port>/mcp/ccp4i2 \\
        --header "Authorization: Bearer <token>"

The path is scoped (/mcp/ccp4i2, not /mcp), as the REST API is
(/api/ccp4i2), so an application that serves CCP4i2 can serve MCP servers of
its own beside it.

Served on the desktop (CCP4I2_LOCAL_SESSION_TOKEN set), and elsewhere only
when the deployment asks for it with CCP4I2_MCP=1; CCP4I2_MCP=0 turns it off
anywhere. A deployment opts in to a new endpoint; it does not acquire one by
updating CCP4i2.

Stateless HTTP, because the desktop app runs uvicorn with two workers and
consecutive requests may reach either. These requests do not pass through
Django's middleware, but every tool call goes back to this server's REST API
over loopback with the caller's own Authorization, so the API's
authentication and permissions decide everything, here as for the app. The
caller's address goes with it as X-Forwarded-For, so the API's logs show who
called rather than 127.0.0.1. At /mcp/ccp4i2 itself, a request with no
Authorization is refused, and on the desktop the session token is checked,
before a session can even list tools.
"""
import hmac
import logging
import os

logger = logging.getLogger(f"ccp4i2:{__name__}")

DEFAULT_MCP_PATH = "/mcp/ccp4i2"


def mcp_path(setting):
    """The mount point: CCP4I2_MCP_PATH if given, else the default.

    This is an ASGI prefix split, not a Django route, so a host that composes
    CCP4i2 into a larger application (as Materia does) cannot re-route it the
    way it can a URLconf; this setting is its lever. A trailing slash is
    dropped, and an empty or bare "/" value falls back to the default, so a
    misconfigured setting cannot mount the facade over the host's root.
    """
    path = (setting or "").strip().rstrip("/")
    if not path:
        return DEFAULT_MCP_PATH
    return path if path.startswith("/") else "/" + path


MCP_PATH = mcp_path(os.environ.get("CCP4I2_MCP_PATH"))

_ON = ("1", "true", "yes", "on")
_OFF = ("0", "false", "no", "off")


def _header(scope, name):
    name = name.lower().encode()
    for key, value in scope.get("headers") or []:
        if key.lower() == name:
            return value.decode("latin-1")
    return None


async def _refuse(send, status, text):
    await send({"type": "http.response.start", "status": status,
                "headers": [(b"content-type", b"text/plain; charset=utf-8")]})
    await send({"type": "http.response.body", "body": text.encode()})


def mcp_enabled():
    """Whether this process serves MCP: CCP4I2_MCP decides when set, else
    only the desktop (which sets CCP4I2_LOCAL_SESSION_TOKEN) does."""
    setting = os.environ.get("CCP4I2_MCP", "").strip().lower()
    if setting in _ON:
        return True
    if setting in _OFF:
        return False
    return bool(os.environ.get("CCP4I2_LOCAL_SESSION_TOKEN"))


def _forwarded_for(scope):
    """X-Forwarded-For for the loopback calls: any chain the request carried,
    then the address it came from (as a proxy appends its own peer)."""
    chain = _header(scope, "x-forwarded-for")
    peer = (scope.get("client") or (None,))[0]
    parts = [p for p in (chain, peer) if p]
    return ", ".join(parts) or None


def with_mcp(django_app):
    """The ASGI application: MCP_PATH to the MCP facade, the rest to Django.

    Django alone when MCP is not enabled here (mcp_enabled), or when the MCP
    server cannot be built, which is said in the log.
    """
    if not mcp_enabled():
        return django_app
    # Never let the agent route stop the app starting (the desktop's readiness
    # probe is this module's import): an old Python, a missing platform
    # package or an SDK change means MCP is off, said in the log, not an app
    # that will not start.
    try:
        from .mcp_server import _caller, server
        mcp_app = server.streamable_http_app(streamable_http_path=MCP_PATH, stateless_http=True,
                                             json_response=True)
    except Exception as err:  # noqa: BLE001 - deliberately broad, see above
        logger.warning("MCP over HTTP is off: %s: %s", type(err).__name__, err)
        return django_app

    async def application(scope, receive, send):
        if scope["type"] == "lifespan":
            # Django takes no lifespan; the MCP app starts its session manager on it
            return await mcp_app(scope, receive, send)
        path = scope.get("path", "")
        if scope["type"] != "http" or not (path == MCP_PATH or path.startswith(MCP_PATH + "/")):
            return await django_app(scope, receive, send)
        authorization = _header(scope, "authorization")
        expected = os.environ.get("CCP4I2_LOCAL_SESSION_TOKEN")
        if expected is not None:
            provided = (authorization or "")[len("Bearer "):] if (authorization or "").startswith(
                "Bearer ") else ""
            if not hmac.compare_digest(provided, expected):
                return await _refuse(send, 401, "This CCP4i2 needs its session token "
                                                "(Help > About): Authorization: Bearer <token>")
        elif not authorization:
            # The API would refuse every tool call anyway; refuse the session too
            return await _refuse(send, 401, "Authorization required: the same credentials "
                                            "as the CCP4i2 REST API")
        host, port = scope.get("server") or ("127.0.0.1", 80)
        if host in ("0.0.0.0", "::"):
            host = "127.0.0.1"
        token = _caller.set({"base": f"http://{host}:{port}/api/ccp4i2",
                             "authorization": authorization,
                             "forwarded_for": _forwarded_for(scope)})
        try:
            await mcp_app(scope, receive, send)
        finally:
            _caller.reset(token)

    return application
