"""
Management command to populate rdkit_smiles and molecular_weight
for compounds that have SMILES but are missing computed fields.

Usage:
    python manage.py populate_inchi
    python manage.py populate_inchi --dry-run
    python manage.py populate_inchi --force   # recalculate ALL, even if already set
"""

from django.apps import apps
from django.core.management.base import BaseCommand


class Command(BaseCommand):
    help = 'Populate rdkit_smiles and molecular_weight for compounds with SMILES'

    def add_arguments(self, parser):
        parser.add_argument(
            '--dry-run', action='store_true',
            help='Show what would be updated without making changes',
        )
        parser.add_argument(
            '--force', action='store_true',
            help='Recalculate all fields even if already populated',
        )

    def handle(self, *args, **options):
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors
        except ImportError:
            self.stderr.write(self.style.ERROR('RDKit is not available'))
            return

        # Check if InChI support is available
        MolToInchi = None
        try:
            from rdkit.Chem import inchi as _inchi_mod
            if getattr(_inchi_mod, 'INCHI_AVAILABLE', False):
                MolToInchi = _inchi_mod.MolToInchi
        except (ImportError, AttributeError):
            pass

        if MolToInchi:
            self.stdout.write('InChI support: available')
        else:
            self.stdout.write('InChI support: NOT available (will populate rdkit_smiles only)')

        dry_run = options['dry_run']
        force = options['force']

        Compound = apps.get_model('registry', 'Compound')

        # Find compounds needing update
        qs = Compound.objects.filter(smiles__isnull=False).exclude(smiles='')
        if not force:
            qs = qs.filter(rdkit_smiles__isnull=True)

        total = qs.count()
        if total == 0:
            self.stdout.write(self.style.SUCCESS(
                'All compounds already have rdkit_smiles populated.'
            ))
            return

        self.stdout.write(f'{"[DRY RUN] " if dry_run else ""}Processing {total} compounds...')

        updated = 0
        errors = 0

        for compound in qs.iterator():
            try:
                mol = Chem.MolFromSmiles(compound.smiles)
                if mol is None:
                    self.stderr.write(
                        f'  {compound.formatted_id}: invalid SMILES "{compound.smiles}"'
                    )
                    errors += 1
                    continue

                new_rdkit_smiles = Chem.MolToSmiles(mol, canonical=True)
                new_mw = Descriptors.MolWt(mol)
                new_inchi = MolToInchi(mol) if MolToInchi else None

                if dry_run:
                    self.stdout.write(
                        f'  {compound.formatted_id}: '
                        f'rdkit_smiles={"SET" if new_rdkit_smiles else "FAILED"} '
                        f'mw={new_mw:.2f}'
                        f'{" inchi=SET" if new_inchi else ""}'
                    )
                else:
                    update_fields = []
                    if force or not compound.rdkit_smiles:
                        compound.rdkit_smiles = new_rdkit_smiles
                        update_fields.append('rdkit_smiles')
                    if force or not compound.molecular_weight:
                        compound.molecular_weight = new_mw
                        update_fields.append('molecular_weight')
                    if new_inchi and (force or not compound.inchi):
                        compound.inchi = new_inchi
                        update_fields.append('inchi')

                    if update_fields:
                        compound.save(update_fields=update_fields)

                updated += 1
            except Exception as e:
                self.stderr.write(f'  {compound.formatted_id}: {e}')
                errors += 1

        prefix = '[DRY RUN] ' if dry_run else ''
        self.stdout.write(self.style.SUCCESS(
            f'{prefix}Done: {updated} updated, {errors} errors (of {total} total)'
        ))
