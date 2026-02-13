"""
Management command to deduplicate compounds with identical canonical SMILES.

For each group of compounds sharing the same rdkit_smiles AND equivalent
stereo_comment (unset/achiral/empty are treated as equivalent), keeps the
oldest registration (lowest reg_number) and deletes the newer one(s) — but
ONLY if the newer compound has no references from:
  - Assay DataSeries (compound or batch FK)
  - Design Hypotheses (product_compound FK)
  - CCP4i2 Projects (project name containing NCL-XXXXXXXX)

Usage:
    python manage.py deduplicate_compounds --dry-run    # report only
    python manage.py deduplicate_compounds --delete     # actually delete
"""

from collections import defaultdict

from django.apps import apps
from django.db.models import Count, Q
from django.core.management.base import BaseCommand


# stereo_comment values that all mean "no specific stereochemistry info"
UNSPECIFIED_STEREO = {None, '', 'unset', 'achiral'}


class Command(BaseCommand):
    help = 'Remove duplicate compound registrations that have no external references'

    @staticmethod
    def _normalize_stereo(value):
        """Normalize stereo_comment for equivalence comparison.

        Returns None for all "unspecified" values (None, '', 'unset', 'achiral'),
        otherwise returns the value as-is for exact matching.
        """
        if value in UNSPECIFIED_STEREO:
            return None
        return value

    def add_arguments(self, parser):
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument(
            '--dry-run', action='store_true',
            help='Report duplicates and what would be deleted (no changes made)',
        )
        group.add_argument(
            '--delete', action='store_true',
            help='Actually delete unreferenced duplicate compounds',
        )

    def handle(self, *args, **options):
        dry_run = options['dry_run']
        do_delete = options['delete']

        Compound = apps.get_model('registry', 'Compound')

        # Try to load reference models — they may not be installed
        DataSeries = None
        Hypothesis = None
        Project = None

        try:
            DataSeries = apps.get_model('assays', 'DataSeries')
        except LookupError:
            self.stdout.write('  (assays app not installed — skipping DataSeries check)')

        try:
            Hypothesis = apps.get_model('assays', 'Hypothesis')
        except LookupError:
            pass  # already noted above if assays missing

        try:
            Project = apps.get_model('db', 'Project')
        except LookupError:
            self.stdout.write('  (ccp4i2 db app not installed — skipping Project check)')

        # Find rdkit_smiles with more than one compound
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

        prefix = '[DRY RUN] ' if dry_run else ''
        self.stdout.write(f'{prefix}Processing {total_groups} duplicate group(s)...\n')

        deleted_count = 0
        skipped_count = 0
        blocked_count = 0

        for dup in duplicates:
            compounds = list(
                Compound.objects
                .filter(rdkit_smiles=dup['rdkit_smiles'])
                .select_related('target')
                .order_by('reg_number')
            )

            # Sub-group by normalized stereo_comment so that e.g.
            # R-enantiomer and S-enantiomer are NOT treated as duplicates
            stereo_groups = defaultdict(list)
            for c in compounds:
                key = self._normalize_stereo(c.stereo_comment)
                stereo_groups[key].append(c)

            # Only process sub-groups with actual duplicates (2+ compounds)
            for stereo_key, group in stereo_groups.items():
                if len(group) < 2:
                    continue

                keeper = group[0]
                candidates = group[1:]

                smiles_display = (
                    f'{dup["rdkit_smiles"][:60]}...'
                    if len(dup['rdkit_smiles']) > 60
                    else dup['rdkit_smiles']
                )
                stereo_display = stereo_key or 'unspecified'
                self.stdout.write(
                    f'  Structure: {smiles_display}  '
                    f'[stereo: {stereo_display}]'
                )
                self.stdout.write(
                    f'  Keep: {keeper.formatted_id} (reg #{keeper.reg_number}, '
                    f'target: {keeper.target.name if keeper.target else "?"}, '
                    f'batches: {keeper.batches.count()}, '
                    f'registered: {keeper.registered_at:%Y-%m-%d})'
                )

                # Collect aliases to add to the keeper from deleted candidates
                new_aliases = []

                for candidate in candidates:
                    blockers = self._check_references(
                        candidate, DataSeries, Hypothesis, Project
                    )

                    batch_count = candidate.batches.count()
                    label = (
                        f'{candidate.formatted_id} (reg #{candidate.reg_number}, '
                        f'target: {candidate.target.name if candidate.target else "?"}, '
                        f'batches: {batch_count}, '
                        f'registered: {candidate.registered_at:%Y-%m-%d})'
                    )

                    if blockers:
                        self.stdout.write(self.style.WARNING(
                            f'  SKIP: {label}'
                        ))
                        for reason in blockers:
                            self.stdout.write(f'         - {reason}')
                        blocked_count += 1
                    else:
                        # Collect the deleted compound's identifiers as aliases
                        new_aliases.append(candidate.formatted_id)
                        new_aliases.append(str(candidate.reg_number))
                        if candidate.supplier_ref and candidate.supplier_ref.strip():
                            new_aliases.append(candidate.supplier_ref.strip())

                        if do_delete:
                            candidate.delete()
                            self.stdout.write(self.style.SUCCESS(
                                f'  DELETED: {label}'
                            ))
                            deleted_count += 1
                        else:
                            self.stdout.write(self.style.SUCCESS(
                                f'  WOULD DELETE: {label}'
                            ))
                            deleted_count += 1

                # Add aliases from deleted compounds to the keeper
                if new_aliases:
                    existing_aliases = list(keeper.aliases or [])
                    # Deduplicate: only add aliases not already present
                    existing_set = {str(a).lower() for a in existing_aliases}
                    added = []
                    for alias in new_aliases:
                        if alias.lower() not in existing_set:
                            existing_aliases.append(alias)
                            existing_set.add(alias.lower())
                            added.append(alias)

                    if added:
                        if do_delete:
                            keeper.aliases = existing_aliases
                            keeper.save(update_fields=['aliases'])
                            self.stdout.write(
                                f'  Aliases added to {keeper.formatted_id}: '
                                f'{", ".join(added)}'
                            )
                        else:
                            self.stdout.write(
                                f'  Would add aliases to {keeper.formatted_id}: '
                                f'{", ".join(added)}'
                            )

                self.stdout.write('')

        self.stdout.write(self.style.WARNING(
            f'\n{prefix}Summary: '
            f'{deleted_count} {"deleted" if do_delete else "would delete"}, '
            f'{blocked_count} skipped (have references)'
        ))

    def _check_references(self, compound, DataSeries, Hypothesis, Project):
        """Check if a compound has any external references that block deletion."""
        blockers = []

        # 1. Check DataSeries referencing this compound directly
        if DataSeries is not None:
            ds_count = DataSeries.objects.filter(compound=compound).count()
            if ds_count:
                blockers.append(f'{ds_count} assay data series reference this compound')

        # 2. Check DataSeries referencing any batch of this compound
        if DataSeries is not None:
            batch_ids = list(compound.batches.values_list('id', flat=True))
            if batch_ids:
                batch_ds_count = DataSeries.objects.filter(
                    batch_id__in=batch_ids
                ).count()
                if batch_ds_count:
                    blockers.append(
                        f'{batch_ds_count} assay data series reference '
                        f'batches of this compound'
                    )

        # 3. Check Hypothesis referencing this compound
        if Hypothesis is not None:
            hyp_count = Hypothesis.objects.filter(product_compound=compound).count()
            if hyp_count:
                blockers.append(f'{hyp_count} hypothesis references this compound')

        # 4. Check CCP4i2 Project names containing this compound's NCL-ID
        #    Handles mixed case (ncl, NCL, Ncl) and variable zero-padding
        #    (NCL-31459, NCL-031459, ... NCL-00031459)
        if Project is not None:
            reg_str = str(compound.reg_number)
            # Build all plausible padded variants (5 to 8 digits)
            q = Q()
            for width in range(len(reg_str), 9):
                padded = reg_str.zfill(width)
                # icontains is case-insensitive, covers ncl/NCL/Ncl
                q |= Q(name__icontains=f'NCL-{padded}')
            project_count = Project.objects.filter(q).count()
            if project_count:
                blockers.append(
                    f'{project_count} CCP4i2 project(s) reference '
                    f'{compound.formatted_id} in their name'
                )

        return blockers
