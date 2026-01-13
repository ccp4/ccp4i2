"""
Azure AD JWT validation middleware for CCP4i2.

This middleware validates Azure AD JWT tokens on incoming requests when
CCP4I2_REQUIRE_AUTH=true is set. It supports both:

1. Authorization header: "Bearer <token>"
2. X-MS-TOKEN-AAD-ACCESS-TOKEN header (set by Azure Container Apps Easy Auth)

Configuration environment variables:
- CCP4I2_REQUIRE_AUTH: Set to "true" to enable authentication (default: false)
- AZURE_AD_TENANT_ID: Your Azure AD tenant ID
- AZURE_AD_CLIENT_ID: Your Azure AD app registration client ID

The middleware validates:
- Token signature (using Azure AD's public keys)
- Token expiration
- Audience (must match client ID)
- Issuer (must match Azure AD tenant)
"""

import os
import json
import time
import logging
import ssl
import certifi
from functools import lru_cache
from typing import Optional, Tuple
from urllib.request import urlopen
from urllib.error import URLError

from django.contrib.auth import get_user_model
from django.http import JsonResponse, HttpRequest

logger = logging.getLogger(__name__)


class AzureADTokenValidator:
    """Validates Azure AD JWT tokens."""

    def __init__(self, tenant_id: str, client_id: str):
        self.tenant_id = tenant_id
        self.client_id = client_id
        self.issuer = f"https://login.microsoftonline.com/{tenant_id}/v2.0"
        self.jwks_uri = f"https://login.microsoftonline.com/{tenant_id}/discovery/v2.0/keys"
        self._keys_cache: Optional[dict] = None
        self._keys_cache_time: float = 0
        self._keys_cache_ttl: float = 3600  # 1 hour

    def _get_signing_keys(self) -> dict:
        """Fetch Azure AD's public signing keys (JWKS)."""
        now = time.time()

        # Return cached keys if still valid
        if self._keys_cache and (now - self._keys_cache_time) < self._keys_cache_ttl:
            return self._keys_cache

        try:
            # Use certifi's certificate bundle for SSL verification
            ssl_context = ssl.create_default_context(cafile=certifi.where())
            with urlopen(self.jwks_uri, timeout=10, context=ssl_context) as response:
                jwks = json.loads(response.read().decode("utf-8"))
                self._keys_cache = {key["kid"]: key for key in jwks.get("keys", [])}
                self._keys_cache_time = now
                logger.debug(f"Fetched {len(self._keys_cache)} signing keys from Azure AD")
                return self._keys_cache
        except URLError as e:
            logger.error(f"Failed to fetch Azure AD signing keys: {e}")
            # Return stale cache if available, otherwise raise
            if self._keys_cache:
                logger.warning("Using stale signing keys cache")
                return self._keys_cache
            raise

    def validate_token(self, token: str) -> Tuple[bool, Optional[dict], Optional[str]]:
        """
        Validate a JWT token.

        Returns:
            Tuple of (is_valid, claims_dict, error_message)
        """
        try:
            import jwt
            from jwt import PyJWK
        except ImportError:
            logger.error("PyJWT not installed. Run: pip install PyJWT[crypto]")
            return False, None, "Server configuration error: PyJWT not installed"

        try:
            # Decode header to get key ID
            unverified_header = jwt.get_unverified_header(token)
            kid = unverified_header.get("kid")

            if not kid:
                return False, None, "Token missing key ID (kid)"

            # Get the signing key
            keys = self._get_signing_keys()
            if kid not in keys:
                # Key not found - try refreshing the cache
                self._keys_cache = None
                keys = self._get_signing_keys()
                if kid not in keys:
                    return False, None, f"Unknown signing key: {kid}"

            key_data = keys[kid]

            # Construct signing key directly from cached JWK data
            # (avoids PyJWKClient making its own HTTPS request)
            jwk = PyJWK.from_dict(key_data)
            signing_key = jwk.key

            # Decode and validate the token
            claims = jwt.decode(
                token,
                signing_key,
                algorithms=["RS256"],
                audience=self.client_id,
                issuer=self.issuer,
                options={
                    "verify_signature": True,
                    "verify_exp": True,
                    "verify_aud": True,
                    "verify_iss": True,
                },
            )

            logger.debug(f"Token validated for subject: {claims.get('sub', 'unknown')}")
            return True, claims, None

        except jwt.ExpiredSignatureError:
            return False, None, "Token has expired"
        except jwt.InvalidAudienceError:
            return False, None, "Invalid token audience"
        except jwt.InvalidIssuerError:
            return False, None, "Invalid token issuer"
        except jwt.InvalidSignatureError:
            return False, None, "Invalid token signature"
        except jwt.DecodeError as e:
            return False, None, f"Token decode error: {e}"
        except Exception as e:
            logger.exception("Unexpected error validating token")
            return False, None, f"Token validation error: {e}"


# Singleton validator instance
_validator: Optional[AzureADTokenValidator] = None


def get_validator() -> Optional[AzureADTokenValidator]:
    """Get or create the token validator instance."""
    global _validator

    if _validator is not None:
        return _validator

    tenant_id = os.environ.get("AZURE_AD_TENANT_ID")
    client_id = os.environ.get("AZURE_AD_CLIENT_ID")

    if not tenant_id or not client_id:
        logger.warning(
            "AZURE_AD_TENANT_ID and AZURE_AD_CLIENT_ID must be set for authentication. "
            "Authentication is disabled."
        )
        return None

    _validator = AzureADTokenValidator(tenant_id, client_id)
    return _validator


def is_auth_required() -> bool:
    """Check if authentication is required."""
    return os.environ.get("CCP4I2_REQUIRE_AUTH", "").lower() in ("true", "1", "yes")


def extract_token(request: HttpRequest) -> Optional[str]:
    """
    Extract JWT token from request.

    Checks in order:
    1. Authorization header (Bearer token)
    2. X-MS-TOKEN-AAD-ACCESS-TOKEN header (Azure Easy Auth)
    """
    # Check Authorization header
    auth_header = request.headers.get("Authorization", "")
    if auth_header.startswith("Bearer "):
        return auth_header[7:]

    # Check Azure Easy Auth header
    easy_auth_token = request.headers.get("X-MS-TOKEN-AAD-ACCESS-TOKEN")
    if easy_auth_token:
        return easy_auth_token

    return None


class AzureADAuthMiddleware:
    """
    Django middleware for Azure AD JWT validation.

    Validates JWT tokens on incoming requests when CCP4I2_REQUIRE_AUTH=true.
    Exempt paths (like health checks) are configurable.
    """

    # Paths that don't require authentication
    EXEMPT_PATHS = [
        "/health",
        "/healthz",
        "/ready",
        "/api/health",
        "/api/ccp4i2/health",
    ]

    def __init__(self, get_response):
        self.get_response = get_response
        self.auth_required = is_auth_required()

        if self.auth_required:
            logger.info("Azure AD authentication is ENABLED")
            validator = get_validator()
            if not validator:
                logger.error(
                    "Authentication required but AZURE_AD_TENANT_ID/AZURE_AD_CLIENT_ID not set!"
                )
        else:
            logger.info("Azure AD authentication is DISABLED")

    def __call__(self, request: HttpRequest):
        # Skip auth if not required
        if not self.auth_required:
            return self.get_response(request)

        # Skip exempt paths
        path = request.path.rstrip("/")
        if any(path == exempt.rstrip("/") for exempt in self.EXEMPT_PATHS):
            return self.get_response(request)

        # Extract token
        token = extract_token(request)
        if not token:
            return JsonResponse(
                {
                    "success": False,
                    "error": "Authentication required. Provide Authorization: Bearer <token>",
                },
                status=401,
            )

        # Validate token
        validator = get_validator()
        if not validator:
            return JsonResponse(
                {
                    "success": False,
                    "error": "Server authentication not configured",
                },
                status=500,
            )

        is_valid, claims, error = validator.validate_token(token)
        if not is_valid:
            return JsonResponse(
                {"success": False, "error": f"Authentication failed: {error}"},
                status=401,
            )

        # Attach claims to request for downstream use
        request.azure_ad_claims = claims
        azure_ad_sub = claims.get("sub")
        request.azure_ad_user_id = azure_ad_sub

        # Try multiple claim fields where email might be found
        email = (
            claims.get("email") or
            claims.get("preferred_username") or
            claims.get("upn") or  # User Principal Name
            claims.get("unique_name")  # Legacy claim
        )

        # Fallback: check X-User-Email header (sent by frontend from MSAL account info)
        # This is secure because:
        # 1. The JWT is already validated - we know WHO this is via cryptographic 'sub'
        # 2. We use 'sub' as the primary key, not email
        # 3. Email is just for display/lookup convenience
        if not email:
            header_email = request.headers.get("X-User-Email")
            if header_email:
                logger.info(f"Using X-User-Email header for sub={azure_ad_sub[:8]}...")
                email = header_email

        if not email:
            # Last resort: generate placeholder from verified sub
            logger.warning(f"No email found for sub={azure_ad_sub}. Claims: {list(claims.keys())}")
            email = f"user_{azure_ad_sub[:8]}@azuread.local"

        request.azure_ad_email = email

        # Get or create Django user
        # Use 'sub' as the unique identifier (cryptographically verified)
        # This prevents email header spoofing from allowing impersonation
        User = get_user_model()
        username = f"aad_{azure_ad_sub[:32]}"  # Stable username from sub
        user, created = User.objects.get_or_create(
            username=username,
            defaults={
                "email": email.lower(),
                "first_name": claims.get("given_name", ""),
                "last_name": claims.get("family_name", ""),
            }
        )

        if created:
            logger.info(f"Created user from Azure AD: {email} (sub={azure_ad_sub[:8]}...)")
        else:
            # Update email and name if changed
            updated = False
            if user.email != email.lower():
                user.email = email.lower()
                updated = True
            if claims.get("given_name") and user.first_name != claims.get("given_name"):
                user.first_name = claims.get("given_name")
                updated = True
            if claims.get("family_name") and user.last_name != claims.get("family_name"):
                user.last_name = claims.get("family_name")
                updated = True
            if updated:
                user.save(update_fields=["email", "first_name", "last_name"])

        request.user = user
        return self.get_response(request)
