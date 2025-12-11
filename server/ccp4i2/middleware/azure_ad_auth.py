# Azure AD Group-Based Authorization Middleware
# Add this to your Django settings and middleware

import os
import logging
from django.http import JsonResponse
from django.conf import settings

logger = logging.getLogger(__name__)

# settings.py additions
AZURE_AD_CONFIG = {
    "CLIENT_ID": os.getenv("AZURE_CLIENT_ID"),
    "CLIENT_SECRET": os.getenv("AZURE_CLIENT_SECRET"),
    "TENANT_ID": os.getenv("AZURE_TENANT_ID"),
    "ADMIN_GROUP_ID": os.getenv("AZURE_ADMIN_GROUP_ID"),
    "USER_GROUP_ID": os.getenv("AZURE_USER_GROUP_ID"),
    "REDIRECT_URI": os.getenv(
        "AZURE_REDIRECT_URI",
        "https://your-app.azurecontainerapps.io/.auth/login/aad/callback",
    ),
}

# Required packages to add to requirements.txt:
# msal==1.26.0
# PyJWT==2.8.0


class AzureADGroupMiddleware:
    """
    Middleware to enforce Azure AD group-based access control
    """

    def __init__(self, get_response):
        self.get_response = get_response
        self.public_paths = [
            "/health/",
            "/.auth/",
            "/static/",
            "/media/",
            "/api/health/",
        ]
        self.admin_paths = [
            "/admin/",
            "/api/admin/",
        ]

    def __call__(self, request):
        # Skip middleware for public paths
        if any(request.path.startswith(path) for path in self.public_paths):
            return self.get_response(request)

        # Check if user is authenticated via Azure AD
        if not self._is_authenticated(request):
            return self._redirect_to_login(request)

        # Get user groups
        user_groups = self._get_user_groups(request)

        # Check admin paths
        if any(request.path.startswith(path) for path in self.admin_paths):
            if not self._is_admin(user_groups):
                return self._access_denied("Admin access required")

        # Check if user has any required group membership
        if not self._has_required_groups(user_groups):
            return self._access_denied(
                "Access denied. Please contact your administrator."
            )

        return self.get_response(request)

    def _is_authenticated(self, request):
        """Check if user is authenticated via Azure AD"""
        # Check for Azure AD authentication headers
        auth_header = request.META.get("HTTP_AUTHORIZATION", "")
        if not auth_header.startswith("Bearer "):
            return False

        token = auth_header[7:]  # Remove 'Bearer ' prefix

        # Validate token (simplified - you might want to use azure-identity library)
        try:
            # Decode and validate JWT token
            import jwt
            from jwt import PyJWKClient

            # Get Microsoft's public keys
            jwks_client = PyJWKClient(
                f"https://login.microsoftonline.com/{settings.AZURE_AD_CONFIG['TENANT_ID']}/discovery/v2.0/keys"
            )

            # Decode token
            header = jwt.get_unverified_header(token)
            key = jwks_client.get_signing_key(header["kid"])
            decoded_token = jwt.decode(token, key.key, algorithms=["RS256"])

            # Store user info in request
            request.azure_user = {
                "oid": decoded_token.get("oid"),
                "groups": decoded_token.get("groups", []),
                "name": decoded_token.get("name"),
                "email": decoded_token.get("preferred_username"),
            }

            return True

        except Exception as e:
            logger.error(f"Token validation failed: {e}")
            return False

    def _get_user_groups(self, request):
        """Get user's Azure AD groups"""
        if hasattr(request, "azure_user"):
            return request.azure_user.get("groups", [])

        # Fallback: query Microsoft Graph API for groups
        try:
            access_token = self._get_access_token()
            user_id = request.azure_user["oid"]

            import requests

            headers = {
                "Authorization": f"Bearer {access_token}",
                "Content-Type": "application/json",
            }

            # Get user's group memberships
            response = requests.get(
                f"https://graph.microsoft.com/v1.0/users/{user_id}/memberOf",
                headers=headers,
            )

            if response.status_code == 200:
                groups_data = response.json()
                return [
                    group["id"]
                    for group in groups_data.get("value", [])
                    if group.get("@odata.type") == "#microsoft.graph.group"
                ]

        except Exception as e:
            logger.error(f"Failed to get user groups: {e}")

        return []

    def _is_admin(self, user_groups):
        """Check if user is in admin group"""
        admin_group_id = settings.AZURE_AD_CONFIG.get("ADMIN_GROUP_ID")
        return admin_group_id in user_groups

    def _has_required_groups(self, user_groups):
        """Check if user has any required group membership"""
        required_groups = [
            settings.AZURE_AD_CONFIG.get("ADMIN_GROUP_ID"),
            settings.AZURE_AD_CONFIG.get("USER_GROUP_ID"),
        ]
        return any(group_id in user_groups for group_id in required_groups if group_id)

    def _get_access_token(self):
        """Get access token for Microsoft Graph API"""
        try:
            from msal import ConfidentialClientApplication

            app = ConfidentialClientApplication(
                client_id=settings.AZURE_AD_CONFIG["CLIENT_ID"],
                client_credential=settings.AZURE_AD_CONFIG["CLIENT_SECRET"],
                authority=f"https://login.microsoftonline.com/{settings.AZURE_AD_CONFIG['TENANT_ID']}",
            )

            result = app.acquire_token_for_client(
                scopes=["https://graph.microsoft.com/.default"]
            )

            if "access_token" in result:
                return result["access_token"]
            else:
                raise Exception(f"Failed to acquire token: {result}")

        except Exception as e:
            logger.error(f"Failed to get access token: {e}")
            raise

    def _redirect_to_login(self, request):
        """Redirect to Azure AD login"""
        login_url = f"https://login.microsoftonline.com/{settings.AZURE_AD_CONFIG['TENANT_ID']}/oauth2/v2.0/authorize"
        params = {
            "client_id": settings.AZURE_AD_CONFIG["CLIENT_ID"],
            "response_type": "code",
            "redirect_uri": settings.AZURE_AD_CONFIG["REDIRECT_URI"],
            "scope": "openid profile email https://graph.microsoft.com/User.Read https://graph.microsoft.com/GroupMember.Read.All",
            "state": request.path,
        }

        from urllib.parse import urlencode

        auth_url = f"{login_url}?{urlencode(params)}"

        return JsonResponse(
            {"error": "Authentication required", "login_url": auth_url}, status=401
        )

    def _access_denied(self, message):
        """Return access denied response"""
        return JsonResponse({"error": "Access denied", "message": message}, status=403)
