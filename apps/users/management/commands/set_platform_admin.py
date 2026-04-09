"""
Management command to grant or revoke platform admin status.

Usage:
    python manage.py set_platform_admin user@example.com --grant
    python manage.py set_platform_admin user@example.com --revoke
    python manage.py set_platform_admin --list
"""

from django.contrib.auth import get_user_model
from django.core.management.base import BaseCommand, CommandError

from users.models import UserProfile

User = get_user_model()


class Command(BaseCommand):
    """Grant or revoke platform admin status."""

    help = "Grant or revoke platform admin status for users"

    def add_arguments(self, parser):
        parser.add_argument(
            'email',
            type=str,
            nargs='?',
            help='User email address',
        )
        parser.add_argument(
            '--grant',
            action='store_true',
            help='Grant platform admin status',
        )
        parser.add_argument(
            '--revoke',
            action='store_true',
            help='Revoke platform admin status',
        )
        parser.add_argument(
            '--list',
            action='store_true',
            help='List all platform admins',
        )

    def handle(self, *args, **options):
        email = options.get('email')
        grant = options['grant']
        revoke = options['revoke']
        list_admins = options['list']

        if list_admins:
            self._list_admins()
            return

        if not email:
            raise CommandError("Email address required (or use --list)")

        if grant and revoke:
            raise CommandError("Cannot use both --grant and --revoke")

        if not grant and not revoke:
            raise CommandError("Must specify --grant or --revoke")

        # Find user by email
        user = User.objects.filter(email__iexact=email).first()
        if not user:
            # Try by username
            user = User.objects.filter(username__iexact=email).first()

        if not user:
            raise CommandError(f"User not found: {email}")

        # Get or create profile
        profile, _ = UserProfile.objects.get_or_create(user=user)

        if grant:
            profile.is_platform_admin = True
            profile.save()
            self.stdout.write(self.style.SUCCESS(
                f"Granted platform admin to: {user.email or user.username}"
            ))
        else:
            profile.is_platform_admin = False
            profile.save()
            self.stdout.write(self.style.WARNING(
                f"Revoked platform admin from: {user.email or user.username}"
            ))

    def _list_admins(self):
        """List all platform admins."""
        admins = User.objects.filter(
            profile__is_platform_admin=True
        ).select_related('profile').order_by('email')

        self.stdout.write("\nPlatform Admins:")
        self.stdout.write("-" * 60)

        if not admins:
            self.stdout.write("  (none)")
        else:
            for user in admins:
                name = f"{user.first_name} {user.last_name}".strip() or user.username
                self.stdout.write(f"  {user.email or user.username} - {name}")

        self.stdout.write("-" * 60)
        self.stdout.write(f"Total: {admins.count()}")
