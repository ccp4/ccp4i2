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

import json
import logging
import os
import ssl
import time
from typing import Optional, Tuple
from urllib.error import URLError
from urllib.request import urlopen

import certifi
from django.contrib.auth import get_user_model
from django.http import HttpRequest, HttpResponse

from ..exceptions import AuthenticationFailed, AuthorizationFailed
from .base import BaseAuthMiddleware

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
    3. Query parameter access_token (for file downloads/anchor links)
    """
    # Check Authorization header
    auth_header = request.headers.get("Authorization", "")
    if auth_header.startswith("Bearer "):
        return auth_header[7:]

    # Check Azure Easy Auth header
    easy_auth_token = request.headers.get("X-MS-TOKEN-AAD-ACCESS-TOKEN")
    if easy_auth_token:
        return easy_auth_token

    # Check query parameter (for file downloads - anchor links don't send headers)
    query_token = request.GET.get("access_token")
    if query_token:
        return query_token

    return None


class AzureADAuthMiddleware(BaseAuthMiddleware):
    """
    Django middleware for Azure AD JWT validation. Cloud auth path.

    Activates when ``CCP4I2_REQUIRE_AUTH=true`` (otherwise no-op, deferring
    to whatever middleware comes next in the chain). When active, every
    non-exempt request must carry a valid Azure AD JWT bearer token; the
    middleware validates it, optionally enforces group-membership rules,
    and gets-or-creates a Django user keyed on the cryptographic ``sub``.

    The dev_admin auto-login path was deliberately split out into a
    separate ``DevAdminMiddleware`` so a misconfigured cloud deploy
    (REQUIRE_AUTH unset) cannot fall through to creating a superuser
    automatically. See ``ccp4i2_auth.middleware.dev_admin`` for the dev
    path; the CCP4i2 settings module is responsible for picking exactly
    one auth middleware based on the deployment shape.
    """

    # Paths that don't require authentication even when this middleware is active.
    EXEMPT_PATHS = [
        "/health",
        "/healthz",
        "/ready",
        "/api/health",
        "/api/ccp4i2/health",
        "/api/ccp4i2/version",
    ]

    def __init__(self, get_response):
        super().__init__(get_response)
        if self.is_active():
            logger.info("Azure AD authentication is ENABLED")
            validator = get_validator()
            if not validator:
                logger.error(
                    "Authentication required but AZURE_AD_TENANT_ID/AZURE_AD_CLIENT_ID not set!"
                )
        else:
            logger.info("Azure AD authentication is DISABLED (CCP4I2_REQUIRE_AUTH not set)")

    def is_active(self) -> bool:
        return is_auth_required()

    def __call__(self, request: HttpRequest) -> HttpResponse:
        # Filter exempt paths *before* delegating to BaseAuthMiddleware so
        # health checks and version probes always bypass authentication.
        if self.is_active():
            path = request.path.rstrip("/")
            if any(path == exempt.rstrip("/") for exempt in self.EXEMPT_PATHS):
                return self.get_response(request)
        return super().__call__(request)

    def authenticate(self, request: HttpRequest):
        token = extract_token(request)
        if not token:
            raise AuthenticationFailed(
                "Authentication required. Provide Authorization: Bearer <token>"
            )

        validator = get_validator()
        if not validator:
            # Server misconfiguration — AZURE_AD_TENANT_ID/CLIENT_ID missing.
            # We surface this as a 401 because to the *caller* the result is
            # the same as a token rejection: no auth happened, retry doesn't
            # help. The error message and the operator-facing log line above
            # tell the operator what to fix.
            raise AuthenticationFailed("Server authentication not configured")

        is_valid, claims, error = validator.validate_token(token)
        if not is_valid:
            raise AuthenticationFailed(error)

        # Attach claims to request for downstream use.
        request.azure_ad_claims = claims
        azure_ad_sub = claims.get("sub")
        request.azure_ad_user_id = azure_ad_sub

        # Groups / Teams membership authorization.
        self._enforce_group_membership(claims, azure_ad_sub)

        email = self._extract_email(claims, request, azure_ad_sub)
        request.azure_ad_email = email

        return self._get_or_create_user(claims, azure_ad_sub, email)

    # --- helper extractions kept private ------------------------------------

    @staticmethod
    def _enforce_group_membership(claims: dict, azure_ad_sub: str) -> None:
        """Raise AuthorizationFailed if ALLOWED_AZURE_AD_GROUPS is set and
        the user is not a member of any allowed group.

        Requires:
        1. Azure AD app configured to emit 'groups' claim (Token Configuration).
        2. ALLOWED_AZURE_AD_GROUPS env var with comma-separated group IDs.
        3. Azure AD Premium P1/P2 (for groups claim in tokens).
        """
        allowed_groups_str = os.environ.get("ALLOWED_AZURE_AD_GROUPS", "")
        allowed_groups = [g.strip() for g in allowed_groups_str.split(",") if g.strip()]
        if not allowed_groups:
            logger.debug("Groups authorization not configured (ALLOWED_AZURE_AD_GROUPS not set)")
            return

        logger.debug(f"Groups authorization enabled. Allowed groups: {allowed_groups}")

        # Group claims overage — user is in >200 groups, Azure AD substitutes
        # _claim_names/_claim_sources for the full list. We can't validate
        # locally; surface a friendly 403 so an admin can move the user to a
        # dedicated app-access group with fewer members.
        if "_claim_names" in claims or "_claim_sources" in claims:
            logger.warning(
                f"Group claims overage detected for user {azure_ad_sub[:8]}. "
                "User is in >200 groups - cannot validate Teams membership from token. "
                "Consider using a dedicated app access group with fewer members."
            )
            raise AuthorizationFailed(
                "Your account has too many group memberships to verify "
                "automatically. Please contact your administrator to be "
                "added to a dedicated application access group."
            )

        user_groups = claims.get("groups", [])
        logger.debug(f"User {azure_ad_sub[:8]}... has groups: {user_groups}")

        if not any(group_id in allowed_groups for group_id in user_groups):
            logger.warning(
                f"Access denied for user {azure_ad_sub[:8]}... - not in authorized groups. "
                f"User groups: {user_groups[:5]}{'...' if len(user_groups) > 5 else ''}, "
                f"Required: {allowed_groups}"
            )
            raise AuthorizationFailed(
                "You are not a member of an authorized group. This application "
                "requires membership in the Newcastle Drug Discovery Unit team. "
                "Please contact your administrator to request access."
            )

        logger.info(f"User {azure_ad_sub[:8]}... authorized via Teams/Groups membership")

    @staticmethod
    def _extract_email(claims: dict, request: HttpRequest, azure_ad_sub: str) -> str:
        # Try multiple claim fields where email might be found.
        email = (
            claims.get("email")
            or claims.get("preferred_username")
            or claims.get("upn")  # User Principal Name
            or claims.get("unique_name")  # Legacy claim
        )

        # Fallback: X-User-Email header (sent by frontend from MSAL account info).
        # Secure because (1) JWT is already validated so we know WHO via 'sub';
        # (2) 'sub' is the primary key, not email; (3) email is just for display.
        if not email:
            header_email = request.headers.get("X-User-Email")
            if header_email:
                logger.info(f"Using X-User-Email header for sub={azure_ad_sub[:8]}...")
                email = header_email

        if not email:
            logger.warning(f"No email found for sub={azure_ad_sub}. Claims: {list(claims.keys())}")
            email = f"user_{azure_ad_sub[:8]}@azuread.local"

        return email

    @staticmethod
    def _get_or_create_user(claims: dict, azure_ad_sub: str, email: str):
        # Extract name from claims (priority: given_name/family_name → 'name' → empty).
        first_name = claims.get("given_name", "")
        last_name = claims.get("family_name", "")
        if not first_name and not last_name:
            full_name = claims.get("name", "")
            if full_name:
                name_parts = full_name.strip().split()
                if len(name_parts) >= 2:
                    first_name = name_parts[0]
                    last_name = " ".join(name_parts[1:])
                elif len(name_parts) == 1:
                    first_name = name_parts[0]

        # Use 'sub' as the unique identifier (cryptographically verified) —
        # prevents email-header spoofing from enabling impersonation.
        User = get_user_model()
        username = f"aad_{azure_ad_sub[:32]}"  # Stable username from sub
        user, created = User.objects.get_or_create(
            username=username,
            defaults={
                "email": email.lower(),
                "first_name": first_name,
                "last_name": last_name,
            },
        )
        if created:
            logger.info(f"Created user from Azure AD: {email} (sub={azure_ad_sub[:8]}...)")

        # Update email and name if they changed or were initially missing.
        updated = False
        if user.email != email.lower():
            user.email = email.lower()
            updated = True
        if first_name and user.first_name != first_name:
            user.first_name = first_name
            updated = True
        if last_name and user.last_name != last_name:
            user.last_name = last_name
            updated = True
        if updated:
            user.save(update_fields=["email", "first_name", "last_name"])

        return user
