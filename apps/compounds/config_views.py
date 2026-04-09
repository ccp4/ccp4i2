"""
Configuration endpoint for compounds app.

Exposes deployment-specific settings to the frontend, allowing the UI
to adapt to different compound ID prefix configurations.
"""

from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import AllowAny
from rest_framework.response import Response

from compounds.formatting import (
    get_compound_prefix,
    get_compound_digits,
    format_compound_id,
)


@api_view(['GET'])
@permission_classes([AllowAny])
def compound_config(request):
    """
    Get deployment configuration for compounds app.

    This endpoint is public (no authentication required) so the frontend
    can fetch configuration before the user is authenticated.

    Returns:
        {
            "compound_id_prefix": "NCL",
            "compound_id_digits": 8,
            "compound_id_example": "NCL-00026042"
        }
    """
    prefix = get_compound_prefix()
    digits = get_compound_digits()
    example_number = 26042

    return Response({
        'compound_id_prefix': prefix,
        'compound_id_digits': digits,
        'compound_id_example': format_compound_id(example_number),
    })
