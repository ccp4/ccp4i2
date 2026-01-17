"""
Cross-Origin Resource Policy (CORP) middleware.

Adds the Cross-Origin-Resource-Policy header to responses to allow resources
to be loaded from pages that have Cross-Origin-Embedder-Policy set.

This is required for Moorhen WebAssembly pages which need COEP: require-corp
to enable SharedArrayBuffer for multi-threading.
"""


class CORPMiddleware:
    """
    Middleware that adds Cross-Origin-Resource-Policy header to all responses.

    Uses 'cross-origin' to allow resources to be loaded from any origin,
    which is needed when frontend and backend are on different ports/domains.
    """

    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        response = self.get_response(request)
        # Add CORP header to allow cross-origin loading under COEP
        response["Cross-Origin-Resource-Policy"] = "cross-origin"
        return response
