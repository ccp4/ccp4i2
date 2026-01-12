"""
Management command to migrate data from legacy ConstructDatabase.

This command imports plasmids, proteins, cassettes, and related data from
the legacy MNApps ConstructDatabase Django app.

Usage:
    python manage.py migrate_legacy_constructs --legacy-db /path/to/legacy/db.sqlite3

    # Dry run (preview only):
    python manage.py migrate_legacy_constructs --legacy-db /path/to/legacy/db.sqlite3 --dry-run

    # With file migration:
    python manage.py migrate_legacy_constructs --legacy-db /path/to/legacy/db.sqlite3 \
        --legacy-media /path/to/legacy/media
"""

import os
import shutil
from pathlib import Path

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from django.db import connections, transaction


class Command(BaseCommand):
    help = 'Migrate data from legacy ConstructDatabase app'

    def add_arguments(self, parser):
        parser.add_argument(
            '--legacy-db',
            type=str,
            required=True,
            help='Path to the legacy SQLite database file'
        )
        parser.add_argument(
            '--legacy-media',
            type=str,
            help='Path to legacy MEDIA_ROOT for file migration'
        )
        parser.add_argument(
            '--dry-run',
            action='store_true',
            help='Preview migration without making changes'
        )
        parser.add_argument(
            '--skip-files',
            action='store_true',
            help='Skip file migration (GenBank files, sequencing results)'
        )

    def handle(self, *args, **options):
        legacy_db = options['legacy_db']
        legacy_media = options.get('legacy_media')
        dry_run = options['dry_run']
        skip_files = options['skip_files']

        if not os.path.exists(legacy_db):
            raise CommandError(f'Legacy database not found: {legacy_db}')

        if legacy_media and not os.path.exists(legacy_media):
            raise CommandError(f'Legacy media directory not found: {legacy_media}')

        self.stdout.write(f'Migrating from: {legacy_db}')
        if dry_run:
            self.stdout.write(self.style.WARNING('DRY RUN - no changes will be made'))

        # Configure legacy database connection
        connections.databases['legacy'] = {
            'ENGINE': 'django.db.backends.sqlite3',
            'NAME': legacy_db,
        }

        try:
            with transaction.atomic():
                self._migrate_reference_data(dry_run)
                self._migrate_construct_projects(dry_run)
                self._migrate_proteins(dry_run)
                self._migrate_plasmids(dry_run, legacy_media, skip_files)
                self._migrate_cassettes(dry_run)
                self._migrate_cassette_uses(dry_run)
                self._migrate_sequencing_results(dry_run, legacy_media, skip_files)
                self._migrate_expression_tags(dry_run)

                if dry_run:
                    raise CommandError('Dry run complete - rolling back')

        except CommandError:
            if dry_run:
                self.stdout.write(self.style.SUCCESS('Dry run complete'))
            else:
                raise
        finally:
            # Clean up legacy connection
            if 'legacy' in connections.databases:
                del connections.databases['legacy']

    def _get_legacy_cursor(self):
        """Get cursor for legacy database."""
        return connections['legacy'].cursor()

    def _migrate_reference_data(self, dry_run: bool):
        """Migrate reference data (already seeded, just verify)."""
        self.stdout.write('Checking reference data...')
        # Reference data is seeded by migrations, so we just log what's in legacy
        cursor = self._get_legacy_cursor()

        cursor.execute('SELECT COUNT(*) FROM ConstructDatabase_expressiontagtype')
        count = cursor.fetchone()[0]
        self.stdout.write(f'  ExpressionTagType: {count} records in legacy')

        cursor.execute('SELECT COUNT(*) FROM ConstructDatabase_protease')
        count = cursor.fetchone()[0]
        self.stdout.write(f'  Protease: {count} records in legacy')

        cursor.execute('SELECT COUNT(*) FROM ConstructDatabase_expressiontaglocation')
        count = cursor.fetchone()[0]
        self.stdout.write(f'  ExpressionTagLocation: {count} records in legacy')

    def _migrate_construct_projects(self, dry_run: bool):
        """Migrate Project -> ConstructProject."""
        from compounds.constructs.models import ConstructProject

        self.stdout.write('Migrating ConstructProjects...')
        cursor = self._get_legacy_cursor()

        cursor.execute('''
            SELECT id, name, parent_id, createdAt, createdBy_id
            FROM ConstructDatabase_project
        ''')

        count = 0
        for row in cursor.fetchall():
            id_, name, parent_id, created_at, created_by_id = row
            if not dry_run:
                ConstructProject.objects.update_or_create(
                    id=id_,
                    defaults={
                        'name': name,
                        'parent_id': parent_id,
                    }
                )
            count += 1

        self.stdout.write(f'  Migrated {count} ConstructProjects')

    def _migrate_proteins(self, dry_run: bool):
        """Migrate Protein and ProteinSynonym."""
        from compounds.constructs.models import Protein, ProteinSynonym

        self.stdout.write('Migrating Proteins...')
        cursor = self._get_legacy_cursor()

        # Migrate proteins
        cursor.execute('''
            SELECT id, uniprotId, createdAt, createdBy_id
            FROM ConstructDatabase_protein
        ''')

        count = 0
        for row in cursor.fetchall():
            id_, uniprot_id, created_at, created_by_id = row
            if not dry_run:
                Protein.objects.update_or_create(
                    id=id_,
                    defaults={
                        'uniprot_id': uniprot_id,
                    }
                )
            count += 1

        self.stdout.write(f'  Migrated {count} Proteins')

        # Migrate protein synonyms
        cursor.execute('''
            SELECT id, name, protein_id, createdAt, createdBy_id
            FROM ConstructDatabase_proteinsynonym
        ''')

        count = 0
        for row in cursor.fetchall():
            id_, name, protein_id, created_at, created_by_id = row
            if not dry_run:
                ProteinSynonym.objects.update_or_create(
                    id=id_,
                    defaults={
                        'name': name,
                        'protein_id': protein_id,
                    }
                )
            count += 1

        self.stdout.write(f'  Migrated {count} ProteinSynonyms')

    def _migrate_plasmids(self, dry_run: bool, legacy_media: str | None, skip_files: bool):
        """Migrate Plasmid records and files."""
        from compounds.constructs.models import Plasmid

        self.stdout.write('Migrating Plasmids...')
        cursor = self._get_legacy_cursor()

        cursor.execute('''
            SELECT id, ncnId, name, parent_id, snapgeneFile, createdAt, createdBy_id
            FROM ConstructDatabase_plasmid
        ''')

        count = 0
        file_count = 0
        for row in cursor.fetchall():
            id_, ncn_id, name, parent_id, snapgene_file, created_at, created_by_id = row

            if not dry_run:
                plasmid, created = Plasmid.objects.update_or_create(
                    id=id_,
                    defaults={
                        'ncn_id': ncn_id,
                        'name': name,
                        'parent_id': parent_id,
                    }
                )

                # Migrate file if exists
                if snapgene_file and legacy_media and not skip_files:
                    src_path = Path(legacy_media) / snapgene_file
                    if src_path.exists():
                        # Create destination directory
                        dst_dir = Path(settings.MEDIA_ROOT) / 'constructs' / plasmid.formatted_id
                        dst_dir.mkdir(parents=True, exist_ok=True)
                        dst_path = dst_dir / src_path.name

                        if not dst_path.exists():
                            shutil.copy2(src_path, dst_path)
                            file_count += 1

                        # Update plasmid with new file path
                        plasmid.genbank_file.name = str(
                            Path('constructs') / plasmid.formatted_id / src_path.name
                        )
                        plasmid.save(update_fields=['genbank_file'])

            count += 1

        self.stdout.write(f'  Migrated {count} Plasmids, {file_count} files')

    def _migrate_cassettes(self, dry_run: bool):
        """Migrate Cassette records."""
        from compounds.constructs.models import Cassette

        self.stdout.write('Migrating Cassettes...')
        cursor = self._get_legacy_cursor()

        cursor.execute('''
            SELECT id, protein_id, start, "end", createdAt, createdBy_id
            FROM ConstructDatabase_cassette
        ''')

        count = 0
        for row in cursor.fetchall():
            id_, protein_id, start, end, created_at, created_by_id = row
            if not dry_run:
                Cassette.objects.update_or_create(
                    id=id_,
                    defaults={
                        'protein_id': protein_id,
                        'start': start,
                        'end': end,
                    }
                )
            count += 1

        self.stdout.write(f'  Migrated {count} Cassettes')

    def _migrate_cassette_uses(self, dry_run: bool):
        """Migrate CassetteUse records."""
        from compounds.constructs.models import CassetteUse

        self.stdout.write('Migrating CassetteUses...')
        cursor = self._get_legacy_cursor()

        cursor.execute('''
            SELECT id, cassette_id, plasmid_id, alignmentFile, createdAt, createdBy_id
            FROM ConstructDatabase_cassetteuse
        ''')

        count = 0
        for row in cursor.fetchall():
            id_, cassette_id, plasmid_id, alignment_file, created_at, created_by_id = row
            if not dry_run:
                CassetteUse.objects.update_or_create(
                    id=id_,
                    defaults={
                        'cassette_id': cassette_id,
                        'plasmid_id': plasmid_id,
                        # alignment_file migration would need file copy
                    }
                )
            count += 1

        self.stdout.write(f'  Migrated {count} CassetteUses')

    def _migrate_sequencing_results(self, dry_run: bool, legacy_media: str | None, skip_files: bool):
        """Migrate SequencingResult records and files."""
        from compounds.constructs.models import SequencingResult

        self.stdout.write('Migrating SequencingResults...')
        cursor = self._get_legacy_cursor()

        cursor.execute('''
            SELECT id, cassetteUse_id, plasmid_id, file, createdAt, createdBy_id
            FROM ConstructDatabase_sequencingresult
        ''')

        count = 0
        for row in cursor.fetchall():
            id_, cassette_use_id, plasmid_id, file_path, created_at, created_by_id = row
            if not dry_run:
                SequencingResult.objects.update_or_create(
                    id=id_,
                    defaults={
                        'cassette_use_id': cassette_use_id,
                        'plasmid_id': plasmid_id,
                        # file migration would need file copy
                    }
                )
            count += 1

        self.stdout.write(f'  Migrated {count} SequencingResults')

    def _migrate_expression_tags(self, dry_run: bool):
        """Migrate ExpressionTag records."""
        from compounds.constructs.models import (
            ExpressionTag, ExpressionTagType, Protease, ExpressionTagLocation
        )

        self.stdout.write('Migrating ExpressionTags...')
        cursor = self._get_legacy_cursor()

        # Note: Legacy had a bug where cassetteUse referenced ExpressionTagType
        # We need to handle this during migration
        cursor.execute('''
            SELECT id, expressionTagType_id, protease_id, cassetteUse_id, createdAt, createdBy_id
            FROM ConstructDatabase_expressiontag
        ''')

        count = 0
        skipped = 0
        for row in cursor.fetchall():
            id_, tag_type_id, protease_id, cassette_use_id, created_at, created_by_id = row

            # Skip if we can't find valid references (due to legacy FK bug)
            if not dry_run:
                try:
                    # Map legacy IDs to new model IDs
                    # This is a simplified approach - may need adjustment based on actual legacy data
                    tag_type = ExpressionTagType.objects.filter(id=tag_type_id).first()
                    protease = Protease.objects.filter(id=protease_id).first() if protease_id else None

                    if tag_type:
                        ExpressionTag.objects.update_or_create(
                            id=id_,
                            defaults={
                                'expression_tag_type': tag_type,
                                'protease': protease,
                                # cassette_use_id needs mapping - legacy had FK bug
                            }
                        )
                        count += 1
                    else:
                        skipped += 1
                except Exception as e:
                    self.stdout.write(
                        self.style.WARNING(f'  Skipped ExpressionTag {id_}: {e}')
                    )
                    skipped += 1
            else:
                count += 1

        self.stdout.write(f'  Migrated {count} ExpressionTags, skipped {skipped}')
