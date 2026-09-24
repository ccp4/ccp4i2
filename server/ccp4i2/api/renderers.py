"""JSON rendering for the API.

`SafeJSONRenderer` is the last line of defence against a non-finite float
reaching the renderer. See `ccp4i2.lib.json_safety` for why one can, and why
neither relaxing ``STRICT_JSON`` nor catching the error in the view works.
"""

import logging

from rest_framework.renderers import JSONRenderer

from ..lib.json_safety import replace_non_finite

logger = logging.getLogger(__name__)


class SafeJSONRenderer(JSONRenderer):
    """A strict JSON renderer that will not 500 on NaN or infinity.

    The common path is untouched: the response is rendered strictly, exactly as
    ``JSONRenderer`` would, and pays nothing for this class. Only when that
    raises — which for strict JSON means a float that has no JSON spelling — is
    the payload walked, the offending values replaced by null, and the render
    retried. A ``ValueError`` from any other cause survives the retry and
    propagates as before.

    Null, not omission, is the right repair here: this sees arbitrary payloads
    and cannot know whether a caller distinguishes "absent" from "present but
    unrepresentable". Code that does know — the KPI maps — drops the key
    upstream, and never reaches this path.
    """

    def render(self, data, accepted_media_type=None, renderer_context=None):
        try:
            return super().render(data, accepted_media_type, renderer_context)
        except ValueError:
            request = (renderer_context or {}).get("request")
            logger.warning(
                "Response for %s contains a float JSON cannot represent "
                "(NaN or infinity); serving null in its place. This is a data "
                "defect upstream, not a rendering choice.",
                getattr(request, "path", "<unknown path>"),
                exc_info=True,
            )
            return super().render(
                replace_non_finite(data), accepted_media_type, renderer_context
            )
