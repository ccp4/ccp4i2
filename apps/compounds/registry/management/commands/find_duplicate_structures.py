"""
Management command to find compounds with duplicate canonical SMILES.

Groups compounds by rdkit_smiles and reports any groups with more than
one compound, so they can be merged manually.

Usage:
    python manage.py find_duplicate_structures
"""

from django.apps import apps
from django.db.models import Count
from django.core.management.base import BaseCommand


class Command(BaseCommand):
    help = 'Find compounds that share the same canonical SMILES (duplicate structures)'

    def handle(self, *args, **options):
        Compound = apps.get_model('registry', 'Compound')

        # First check how many compounds are missing rdkit_smiles
        missing = Compound.objects.filter(
            smiles__isnull=False
        ).exclude(smiles='').filter(rdkit_smiles__isnull=True).count()

        if missing:
            self.stderr.write(self.style.WARNING(
                f'{missing} compounds have SMILES but no rdkit_smiles. '
                f'Run "python manage.py populate_inchi" first to compute canonical SMILES.'
            ))

        # Find rdkit_smiles values that appear more than once
        duplicates = (
            Compound.objects
            .filter(rdkit_smiles__isnull=False)
            .values('rdkit_smiles')
            .annotate(count=Count('id'))
            .filter(count__gt=1)
            .order_by('-count')
        )

        total_groups = duplicates.count()
        if total_groups == 0:
            self.stdout.write(self.style.SUCCESS('No duplicate structures found.'))
            return

        self.stdout.write(self.style.WARNING(
            f'Found {total_groups} structure(s) registered more than once:\n'
        ))

        total_duplicates = 0
        for dup in duplicates:
            compounds = (
                Compound.objects
                .filter(rdkit_smiles=dup['rdkit_smiles'])
                .select_related('target')
                .order_by('reg_number')
            )
            total_duplicates += compounds.count() - 1  # -1 because one is the "original"

            self.stdout.write(
                f'  Canonical SMILES: {dup["rdkit_smiles"]}'
            )
            self.stdout.write(f'  {dup["count"]} registrations:')
            for c in compounds:
                target_name = c.target.name if c.target else '(no target)'
                batch_count = c.batches.count()
                self.stdout.write(
                    f'    {c.formatted_id} (reg #{c.reg_number}) '
                    f'- target: {target_name}, '
                    f'batches: {batch_count}, '
                    f'registered: {c.registered_at:%Y-%m-%d}'
                )
            self.stdout.write('')

        self.stdout.write(self.style.WARNING(
            f'Summary: {total_groups} duplicate groups, '
            f'{total_duplicates} extra registration(s) to review.'
        ))
