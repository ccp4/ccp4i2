#!/usr/bin/env python3
"""
Transform legacy MNApps fixtures to new compounds app schema.

This script transforms fixtures from the old RegisterCompounds/AssayCompounds
Django apps to the new compounds.registry/compounds.assays schema.

Usage:
    python transform_fixtures.py <input_dir> <output_dir>

Where input_dir contains:
    - 20260105-14-57-auth.json
    - 20260105-14-57-RegisterCompounds.json
    - 20260105-14-57-AssayCompounds.json

Output will be:
    - auth.json (cleaned auth.user records)
    - compounds_registry.json (suppliers, targets, compounds, batches, etc.)
    - compounds_assays.json (protocols, assays, data_series, etc.)
"""

from __future__ import annotations

import json
import re
import sys
import os
from pathlib import Path
from datetime import datetime
from typing import Optional


# Datetime fields that need timezone info
DATETIME_FIELDS = {
    'created_at', 'updated_at', 'modified_at', 'registered_at', 'uploaded_at',
    'creation_date', 'reg_date', 'modified_date',
}


def make_timezone_aware(value: str) -> str:
    """Add UTC timezone to naive datetime strings."""
    if not value or not isinstance(value, str):
        return value
    # Already has timezone info
    if value.endswith('Z') or '+' in value[-6:]:
        return value
    # Check if it looks like a datetime (YYYY-MM-DD or similar)
    if re.match(r'^\d{4}-\d{2}-\d{2}', value):
        # Append Z for UTC
        return value + 'Z'
    return value


# Model name mappings: old -> new
MODEL_MAPPINGS = {
    # RegisterCompounds -> compounds.registry
    'RegisterCompounds.regprojects': 'registry.target',
    'RegisterCompounds.regsuppliers': 'registry.supplier',
    'RegisterCompounds.regdata': 'registry.compound',
    'RegisterCompounds.regbatch': 'registry.batch',
    'RegisterCompounds.regbatchqcfile': 'registry.batchqcfile',
    'RegisterCompounds.regdatatemplate': 'registry.compoundtemplate',

    # AssayCompounds -> compounds.assays
    'AssayCompounds.dilutionseries': 'assays.dilutionseries',
    'AssayCompounds.protocol': 'assays.protocol',
    'AssayCompounds.experiment': 'assays.assay',
    'AssayCompounds.dataseries': 'assays.dataseries',
    'AssayCompounds.analysis': 'assays.analysisresult',
    'AssayCompounds.hypothesis': 'assays.hypothesis',
    # Note: tabularview is not being migrated (saved queries can be recreated)
}


def clean_json_content(content: str) -> str:
    """Remove Django debug output from fixture files."""
    start = content.find('[')
    if start == -1:
        raise ValueError("No JSON array found in file")
    return content[start:]


def transform_field_name(old_name: str) -> str:
    """Transform camelCase field names to snake_case."""
    mappings = {
        # RegisterCompounds field mappings
        'project_name': 'name',  # RegProjects -> Target
        'creation_date': 'created_at',
        'reg_id': 'reg_number',
        'project_id': 'target',  # FK reference
        'reg_date': 'registered_at',
        'labbook_id': 'labbook_number',
        'page_id': 'page_number',
        'compound_id': 'compound_number',
        'user_name': 'legacy_registered_by',
        'mw': 'molecular_weight',
        'inchi1': 'inchi',
        'inchi1_qualifier': None,  # Drop this field
        'stereo_comments': 'stereo_comment',
        'pklFile': None,  # Drop pickle files
        'svgFile': 'svg_file',
        'rdkitSmiles': 'rdkit_smiles',
        'aliases': None,  # Drop aliases for now
        'original_id': None,  # Legacy import tracking
        'modified_date': 'modified_at',
        'regdata_id': 'compound',  # FK to compound
        'batch_number': 'batch_number',
        'mwt': 'molecular_weight',
        'batch_id': 'batch',  # FK to batch
        'mol2d': 'mol2d',

        # AssayCompounds field mappings
        'shortName': None,  # Drop - not needed
        'created_by': 'created_by',  # Will need special handling for natural keys
        'created_at': 'created_at',
        'analysisMethod': 'analysis_method',
        'preferredDilutions': 'preferred_dilutions',
        'pherastarTableChoice': 'pherastar_table',
        'dataFile': 'data_file',
        'experiment': 'assay',  # FK renamed
        'compoundName': 'compound_name',
        'startColumn': 'start_column',
        'endColumn': 'end_column',
        'dilutionSeries': 'dilution_series',
        'extractedData': 'extracted_data',
        'skipPoints': 'skip_points',
        'svgFile': 'svg_file',
        'modelMinyURL': 'model_url',
        'updated_at': 'updated_at',
        'productRegId': 'product_compound',  # FK to compound
        'completionNotes': 'completion_notes',
        'project': 'target',  # FK renamed
    }
    return mappings.get(old_name, old_name)


def transform_file_path(old_path: str, model_type: str) -> str:
    """Transform old file paths to new structure."""
    if not old_path:
        return old_path

    # Map old app prefixes to new
    path_mappings = {
        'RegisterCompounds/svg/': 'compounds/registry/svg/',
        'RegBatchQCFile_': 'compounds/registry/qc/',
        'AssayCompounds/Experiments/': 'compounds/assays/data/',
        'AssayCompounds/svg/': 'compounds/assays/svg/',
    }

    for old_prefix, new_prefix in path_mappings.items():
        if old_prefix in old_path:
            return old_path.replace(old_prefix, new_prefix)

    return old_path


def transform_record(record: dict, user_lookup: dict) -> dict | None:
    """Transform a single fixture record to new schema."""
    old_model = record['model']

    # Skip models we're not migrating
    if old_model not in MODEL_MAPPINGS:
        return None

    new_model = MODEL_MAPPINGS[old_model]
    pk = record['pk']
    old_fields = record['fields']
    new_fields = {}

    for old_key, value in old_fields.items():
        new_key = transform_field_name(old_key)

        # Skip dropped fields
        if new_key is None:
            continue

        # Transform file paths
        if 'file' in old_key.lower() or old_key in ('svgFile', 'dataFile'):
            value = transform_file_path(value, new_model)

        # Handle natural key references (created_by is [username])
        # Note: We remove created_by from fixtures because:
        # 1. Django User natural key is username, not email
        # 2. The legacy system stored display names, not usernames
        # 3. Users will be re-linked via email when they next log in (JIT)
        if old_key == 'created_by' and isinstance(value, list):
            if value and len(value) > 0:
                # Store original value for reference in legacy field if it exists
                if 'legacy_registered_by' not in new_fields and new_model in ('registry.compound',):
                    new_fields['legacy_registered_by'] = value[0]
            # Don't include created_by in fixture - will be null
            value = None

        # Normalize stereo_comment to lowercase to match model choices
        if new_key == 'stereo_comment':
            if value:
                value = value.lower().replace(' ', '_').replace('-', '_')
            else:
                value = 'unset'  # Default for null values

        # Normalize analysis_method to lowercase with underscores
        if new_key == 'analysis_method' and value:
            value = value.lower().replace('-', '_').replace(' ', '_')

        # Normalize status fields to lowercase
        if new_key == 'status' and value:
            value = value.lower()

        # Unwrap legacy {"data": ...} wrapper from JSONFields
        # The old system wrapped arrays/objects in a {"data": ...} container
        if new_key in ('extracted_data', 'skip_points') and isinstance(value, dict):
            if 'data' in value:
                value = value['data']
            # Also handle extracted_data which may be {"data": {"0": val, "1": val, ...}}
            # Convert indexed object to array for extracted_data
            if new_key == 'extracted_data' and isinstance(value, dict):
                # Convert {"0": 1.2, "1": 3.4, ...} to [1.2, 3.4, ...]
                try:
                    max_idx = max(int(k) for k in value.keys()) if value else -1
                    value = [value.get(str(i)) for i in range(max_idx + 1)]
                except (ValueError, TypeError):
                    # If keys aren't numeric, leave as-is
                    pass

        # Make datetime fields timezone-aware
        if new_key in DATETIME_FIELDS or old_key in DATETIME_FIELDS:
            value = make_timezone_aware(value)

        if new_key:
            new_fields[new_key] = value

    return {
        'model': new_model,
        'pk': pk,
        'fields': new_fields
    }


def build_user_lookup(auth_data: list) -> dict:
    """Build a lookup from username to email for user resolution."""
    lookup = {}
    for record in auth_data:
        if record['model'] == 'auth.user':
            fields = record['fields']
            username = fields.get('username', '')
            email = fields.get('email', '')
            lookup[username] = {
                'email': email,
                'first_name': fields.get('first_name', ''),
                'last_name': fields.get('last_name', ''),
            }
    return lookup


def clean_auth_fixtures(auth_data: list) -> list:
    """Clean auth fixtures - remove passwords, normalize emails."""
    cleaned = []
    for record in auth_data:
        if record['model'] == 'auth.user':
            fields = record['fields'].copy()
            # Clear passwords - users will re-authenticate via Azure AD
            fields['password'] = ''
            # Ensure email is set (use username as fallback if it looks like email)
            if not fields.get('email') and '@' in fields.get('username', ''):
                fields['email'] = fields['username']

            # Make datetime fields timezone-aware
            for dt_field in ('date_joined', 'last_login'):
                if fields.get(dt_field):
                    fields[dt_field] = make_timezone_aware(fields[dt_field])

            cleaned.append({
                'model': 'auth.user',
                'fields': fields
            })
    return cleaned


def transform_fixtures(input_dir: Path, output_dir: Path):
    """Main transformation function."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load auth data first for user lookup
    auth_file = list(input_dir.glob('*auth.json'))[0]
    print(f"Loading auth data from {auth_file}")
    with open(auth_file, 'r') as f:
        auth_content = clean_json_content(f.read())
        auth_data = json.loads(auth_content)

    user_lookup = build_user_lookup(auth_data)
    print(f"Built user lookup with {len(user_lookup)} users")

    # Clean and save auth fixtures
    cleaned_auth = clean_auth_fixtures(auth_data)
    with open(output_dir / 'auth.json', 'w') as f:
        json.dump(cleaned_auth, f, indent=2)
    print(f"Wrote {len(cleaned_auth)} auth records to auth.json")

    # Transform RegisterCompounds
    reg_file = list(input_dir.glob('*RegisterCompounds.json'))[0]
    print(f"\nLoading RegisterCompounds from {reg_file}")
    with open(reg_file, 'r') as f:
        reg_content = clean_json_content(f.read())
        reg_data = json.loads(reg_content)

    registry_records = []
    for record in reg_data:
        transformed = transform_record(record, user_lookup)
        if transformed:
            registry_records.append(transformed)

    # Sort by model to ensure dependencies load in order
    model_order = ['registry.supplier', 'registry.target', 'registry.compound',
                   'registry.batch', 'registry.batchqcfile', 'registry.compoundtemplate']
    registry_records.sort(key=lambda r: model_order.index(r['model']) if r['model'] in model_order else 999)

    with open(output_dir / 'compounds_registry.json', 'w') as f:
        json.dump(registry_records, f, indent=2)
    print(f"Wrote {len(registry_records)} registry records to compounds_registry.json")

    # Transform AssayCompounds
    assay_file = list(input_dir.glob('*AssayCompounds.json'))[0]
    print(f"\nLoading AssayCompounds from {assay_file}")
    with open(assay_file, 'r') as f:
        assay_content = clean_json_content(f.read())
        assay_data = json.loads(assay_content)

    assay_records = []
    for record in assay_data:
        transformed = transform_record(record, user_lookup)
        if transformed:
            assay_records.append(transformed)

    # Sort by model to ensure dependencies load in order
    model_order = ['assays.dilutionseries', 'assays.protocol', 'assays.assay',
                   'assays.analysisresult', 'assays.dataseries', 'assays.hypothesis']
    assay_records.sort(key=lambda r: model_order.index(r['model']) if r['model'] in model_order else 999)

    with open(output_dir / 'compounds_assays.json', 'w') as f:
        json.dump(assay_records, f, indent=2)
    print(f"Wrote {len(assay_records)} assay records to compounds_assays.json")

    # Summary
    print("\n=== Transformation Summary ===")
    print(f"Auth users: {len(cleaned_auth)}")

    registry_counts = {}
    for r in registry_records:
        registry_counts[r['model']] = registry_counts.get(r['model'], 0) + 1
    for model, count in sorted(registry_counts.items()):
        print(f"{model}: {count}")

    assay_counts = {}
    for r in assay_records:
        assay_counts[r['model']] = assay_counts.get(r['model'], 0) + 1
    for model, count in sorted(assay_counts.items()):
        print(f"{model}: {count}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    input_dir = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])

    if not input_dir.exists():
        print(f"Error: Input directory {input_dir} does not exist")
        sys.exit(1)

    transform_fixtures(input_dir, output_dir)
    print("\nTransformation complete!")
