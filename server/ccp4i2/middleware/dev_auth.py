"""
Development authentication middleware.

Auto-authenticates requests in DEBUG mode for local development.
This middleware should ONLY be used in development environments.
"""

from django.conf import settings
from django.contrib.auth import get_user_model


class DevAuthMiddleware:
    """
    Middleware that auto-authenticates requests in development mode.

    If DEBUG=True and the user is not authenticated, this middleware
    will automatically log them in as the dev user (configured via
    DEV_USER_EMAIL setting).
    """

    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        # Only run in DEBUG mode and if user is not already authenticated
        # Check hasattr because this may run before AuthenticationMiddleware
        if settings.DEBUG:
            user = getattr(request, 'user', None)
            if user is None or not user.is_authenticated:
                User = get_user_model()
                dev_email = getattr(settings, 'DEV_USER_EMAIL', 'dev@localhost')

                # Get or create the dev user
                user, created = User.objects.get_or_create(
                    email=dev_email,
                    defaults={
                        'username': dev_email.split('@')[0],
                        'is_staff': True,
                        'is_superuser': True,
                    }
                )

                if created:
                    print(f"Created dev user: {dev_email}")

                request.user = user

        return self.get_response(request)
