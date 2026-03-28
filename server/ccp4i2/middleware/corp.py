# Copyright (C) 2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
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
