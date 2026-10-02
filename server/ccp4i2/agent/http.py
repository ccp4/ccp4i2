"""The MCP facade served by the app itself, at /mcp beside the REST API.

An agent connects with a URL and the same Bearer token as the REST API
(the desktop app shows both in Help > About; they change at every launch):

    claude mcp add --transport http ccp4i2 http://127.0.0.1:<port>/mcp \\
        --header "Authorization: Bearer <token>"

Stateless HTTP, because the desktop app runs uvicorn with two workers and
consecutive requests may reach either. Every tool call goes back to this
server's REST API over loopback with the caller's own Authorization, so the
API's authentication and validation decide everything, here as for the app.
On the desktop (CCP4I2_LOCAL_SESSION_TOKEN set) /mcp also checks the token
itself, as the REST middleware does, before a session can even list tools.
"""
import hmac
import logging
import os

logger = logging.getLogger(f"ccp4i2:{__name__}")

MCP_PATH = "/mcp"


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


def with_mcp(django_app):
    """The ASGI application: /mcp to the MCP facade, the rest to Django.

    If the MCP server cannot be built here, the app serves Django alone and
    says why in its log.
    """
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
        host, port = scope.get("server") or ("127.0.0.1", 80)
        if host in ("0.0.0.0", "::"):
            host = "127.0.0.1"
        token = _caller.set({"base": f"http://{host}:{port}/api/ccp4i2",
                             "authorization": authorization})
        try:
            await mcp_app(scope, receive, send)
        finally:
            _caller.reset(token)

    return application
