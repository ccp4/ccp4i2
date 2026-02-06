#!/usr/bin/env python3
"""
Campaign Uploader - Automated upload of XChem fragment campaign data via i2remote.

This script automates the process of uploading locally-retrieved XChem data
(from filepicker pull) to a CCP4i2 server campaign.

Workflow stages:
    1. scan    - Scan local data, generate campaign_scan.yaml
    2. plan    - Create upload plan, generate campaign_plan.yaml
    3. execute - Upload files, create jobs, generate campaign_status.yaml

Usage:
    # Stage 1: Scan the data
    campaign_uploader scan /path/to/local_data -o campaign_scan.yaml

    # Stage 2: Plan the upload (review campaign_scan.yaml first)
    campaign_uploader plan campaign_scan.yaml -o campaign_plan.yaml

    # Stage 3: Execute the upload
    campaign_uploader execute campaign_plan.yaml --dry-run           # Preview
    campaign_uploader execute campaign_plan.yaml --dry-run --verbose # Show i2remote commands
    campaign_uploader execute campaign_plan.yaml --campaign-id 123   # Actually upload

Design notes:
    - Project names include upload date for uniqueness across sessions
    - NCL-00000000 is used for APO crystals (no compound)
    - NCL-XXXXXXXX compounds: looked up directly by reg_number
    - Non-NCL compounds (Z*, POB*): resolved via SMILES→InChI registry lookup
    - Unresolved compounds are saved to a remaining file for later registration
"""

import argparse
import glob
import os
import re
import sys
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any

import json

from server.ccp4i2.pipelines.servalcat_pipe.script.servalcat_pipe import servalcat_pipe

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False
    yaml = None


# ─────────────────────────────────────────────────────────────────────────────
# Data Classes for YAML Configuration
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class CompoundInfo:
    """Information about a compound from the crystal directory."""
    identifier: str                    # NCL-XXXXXXXX, Z*, POB*, or None
    identifier_type: str               # 'ncl', 'enamine', 'internal', 'unknown', 'apo'
    cif_path: Optional[str] = None     # Path to compound CIF
    smiles: Optional[str] = None       # SMILES string if available
    smiles_path: Optional[str] = None  # Path to .smiles file
    resolved_ncl: Optional[str] = None # NCL number if resolved (future)


@dataclass
class CrystalData:
    """Data for a single crystal directory."""
    crystal_name: str                  # e.g., CDK4CyclinD1-x0002
    crystal_index: int                 # e.g., 2
    dir_path: str                      # Full path to crystal directory

    # Files
    merged_mtz: Optional[str] = None   # Main processed MTZ
    free_mtz: Optional[str] = None     # FreeR MTZ
    unmerged_mtz: Optional[str] = None # Unmerged data
    dimple_pdb: Optional[str] = None   # Dimple output structure
    dimple_mtz: Optional[str] = None   # Dimple output MTZ

    # Compound
    compound: Optional[CompoundInfo] = None

    # Status flags
    has_dimple: bool = False
    is_complete: bool = False          # Has MTZ + compound CIF
    is_ncl_identified: bool = False    # Has NCL-* compound


@dataclass
class ScanResult:
    """Result of scanning local data directory."""
    scan_date: str
    source_dir: str
    target_prefix: str                 # e.g., CDK4CyclinD1
    total_crystals: int

    # Breakdown
    ncl_identified: int                # Ready for MVP
    other_compounds: int               # Z*, POB* - need resolution
    apo_crystals: int                  # No compound (NCL-00000000)
    incomplete: int                    # Missing files

    crystals: List[CrystalData] = field(default_factory=list)

    # Reference candidates (crystals with dimple output for parent project)
    reference_candidates: List[str] = field(default_factory=list)


@dataclass
class UploadAction:
    """A single upload action in the plan."""
    crystal_name: str
    project_name: str                  # Generated project name with date
    compound_id: str                   # NCL-XXXXXXXX or NCL-00000000 for APO

    # Files for SubstituteLigand (unmerged data)
    unmerged_mtz_path: str

    # Files for servalcat_pipe (dimple harvest)
    dimple_pdb: Optional[str] = None   # final.pdb from dimple (XYZIN)
    dimple_mtz: Optional[str] = None   # Output MTZ from dimple (HKLIN)
    compound_cif: Optional[str] = None # Compound CIF (DICT_LIST[0])

    # SMILES for SubstituteLigand (looked up from registry)
    smiles: Optional[str] = None

    # Status
    status: str = 'pending'            # pending, skipped, completed, failed
    error: Optional[str] = None

    # Results (filled in during execute)
    project_id: Optional[int] = None
    sublig_job_id: Optional[int] = None
    servalcat_pipe_job_id: Optional[int] = None


@dataclass
class UploadPlan:
    """Complete upload plan for execution."""
    plan_date: str
    campaign_name: str
    upload_date_tag: str               # YYYYMMDD for project names
    parent_project_name: str

    # Filter settings
    compound_pattern: Optional[str] = None  # Regex used to filter compounds

    # Parent project files
    reference_pdb: Optional[str] = None
    reference_freer: Optional[str] = None

    # Actions
    actions: List[UploadAction] = field(default_factory=list)

    # Summary
    total_actions: int = 0
    ncl_actions: int = 0
    apo_actions: int = 0


# ─────────────────────────────────────────────────────────────────────────────
# Scanning Functions
# ─────────────────────────────────────────────────────────────────────────────

# Patterns for crystal directory and compound files
CRYSTAL_DIR_PATTERN = re.compile(r'^(.+)-x(\d+)$')
NCL_PATTERN = re.compile(r'^NCL-(\d{8})\.cif$', re.IGNORECASE)
ENAMINE_PATTERN = re.compile(r'^Z\d+\.cif$', re.IGNORECASE)
INTERNAL_PATTERN = re.compile(r'^POB\d+\.cif$', re.IGNORECASE)


def find_compound_info(crystal_dir: Path) -> Optional[CompoundInfo]:
    """
    Find compound information from a crystal's compound/ subdirectory.

    Priority:
    1. NCL-XXXXXXXX.cif (DDU registry - ready to use)
    2. Z*.cif (Enamine - needs resolution)
    3. POB*.cif (internal - needs resolution)
    4. merged.cif (fallback, unknown source)
    5. None (APO)
    """
    compound_dir = crystal_dir / 'compound'
    if not compound_dir.exists():
        return None

    cif_files = list(compound_dir.glob('*.cif'))
    if not cif_files:
        return None

    # Look for NCL pattern first
    for cif in cif_files:
        match = NCL_PATTERN.match(cif.name)
        if match:
            # Found NCL identifier
            ncl_num = match.group(1)
            smiles = _read_smiles(compound_dir / f'NCL-{ncl_num}.smiles')
            return CompoundInfo(
                identifier=f'NCL-{ncl_num}',
                identifier_type='ncl',
                cif_path=str(cif),
                smiles=smiles,
                smiles_path=str(compound_dir / f'NCL-{ncl_num}.smiles') if smiles else None,
            )

    # Look for Enamine (Z*) pattern
    for cif in cif_files:
        match = ENAMINE_PATTERN.match(cif.name)
        if match:
            identifier = cif.stem
            smiles = _read_smiles(compound_dir / f'{identifier}.smiles')
            return CompoundInfo(
                identifier=identifier,
                identifier_type='enamine',
                cif_path=str(cif),
                smiles=smiles,
                smiles_path=str(compound_dir / f'{identifier}.smiles') if smiles else None,
            )

    # Look for internal (POB*) pattern
    for cif in cif_files:
        match = INTERNAL_PATTERN.match(cif.name)
        if match:
            identifier = cif.stem
            smiles = _read_smiles(compound_dir / f'{identifier}.smiles')
            return CompoundInfo(
                identifier=identifier,
                identifier_type='internal',
                cif_path=str(cif),
                smiles=smiles,
                smiles_path=str(compound_dir / f'{identifier}.smiles') if smiles else None,
            )

    # Fallback to merged.cif if it exists but no identifiable compound
    merged_cif = compound_dir / 'merged.cif'
    if merged_cif.exists():
        # Try to find any .smiles file
        smiles_files = list(compound_dir.glob('*.smiles'))
        smiles = _read_smiles(smiles_files[0]) if smiles_files else None
        return CompoundInfo(
            identifier='unknown',
            identifier_type='unknown',
            cif_path=str(merged_cif),
            smiles=smiles,
        )

    return None


def _read_smiles(smiles_path: Path) -> Optional[str]:
    """Read SMILES string from file."""
    if smiles_path and smiles_path.exists():
        try:
            return smiles_path.read_text().strip()
        except Exception:
            pass
    return None


def scan_crystal_directory(crystal_dir: Path) -> Optional[CrystalData]:
    """
    Scan a single crystal directory and extract all available data.
    """
    dir_name = crystal_dir.name
    match = CRYSTAL_DIR_PATTERN.match(dir_name)
    if not match:
        return None

    prefix, index_str = match.groups()
    index = int(index_str)

    # Initialize crystal data
    crystal = CrystalData(
        crystal_name=dir_name,
        crystal_index=index,
        dir_path=str(crystal_dir),
    )

    # Find main processed MTZ (crystal_name.mtz, NOT .free.mtz)
    main_mtz_candidates = [
        crystal_dir / f'{dir_name}.mtz',
    ]
    for mtz in main_mtz_candidates:
        if mtz.exists():
            crystal.merged_mtz = str(mtz)
            break

    # Find FreeR MTZ
    free_mtz = crystal_dir / f'{dir_name}.free.mtz'
    if free_mtz.exists():
        crystal.free_mtz = str(free_mtz)

    # Find unmerged MTZ in autoprocessing/
    autoprocessing = crystal_dir / 'autoprocessing'
    if autoprocessing.exists():
        unmerged = list(autoprocessing.glob('*_scaled_unmerged.mtz'))
        if unmerged:
            crystal.unmerged_mtz = str(unmerged[0])

    # Find dimple output
    dimple_pdb = crystal_dir / 'final.pdb'
    dimple_mtz = crystal_dir / 'final.mtz'
    if dimple_pdb.exists():
        crystal.dimple_pdb = str(dimple_pdb)
        crystal.has_dimple = True
    if dimple_mtz.exists():
        crystal.dimple_mtz = str(dimple_mtz)

    # Find compound info
    crystal.compound = find_compound_info(crystal_dir)

    # Set status flags
    crystal.is_complete = bool(crystal.merged_mtz and crystal.compound)
    crystal.is_ncl_identified = (
        crystal.compound is not None and
        crystal.compound.identifier_type == 'ncl'
    )

    return crystal


def scan_data_directory(data_dir: str, prefix_filter: Optional[str] = None) -> ScanResult:
    """
    Scan the entire local data directory and build a ScanResult.
    """
    data_path = Path(data_dir)
    if not data_path.exists():
        raise ValueError(f"Data directory does not exist: {data_dir}")

    crystals = []
    target_prefix = None

    # Find all crystal directories
    for item in sorted(data_path.iterdir()):
        if not item.is_dir():
            continue

        crystal = scan_crystal_directory(item)
        if crystal is None:
            continue

        # Extract target prefix from first crystal
        if target_prefix is None:
            match = CRYSTAL_DIR_PATTERN.match(crystal.crystal_name)
            if match:
                target_prefix = match.group(1)

        # Apply prefix filter if specified
        if prefix_filter and not crystal.crystal_name.startswith(prefix_filter):
            continue

        crystals.append(crystal)

    # Build statistics
    ncl_identified = sum(1 for c in crystals if c.is_ncl_identified)
    apo_crystals = sum(1 for c in crystals if c.has_dimple and not c.compound)
    other_compounds = sum(
        1 for c in crystals
        if c.compound and c.compound.identifier_type in ('enamine', 'internal', 'unknown')
    )
    incomplete = sum(1 for c in crystals if not c.merged_mtz)

    # Find reference candidates (have dimple output)
    reference_candidates = [
        c.crystal_name for c in crystals
        if c.has_dimple and c.dimple_pdb
    ]

    return ScanResult(
        scan_date=datetime.now().isoformat(),
        source_dir=str(data_path.absolute()),
        target_prefix=target_prefix or 'Unknown',
        total_crystals=len(crystals),
        ncl_identified=ncl_identified,
        other_compounds=other_compounds,
        apo_crystals=apo_crystals,
        incomplete=incomplete,
        crystals=crystals,
        reference_candidates=reference_candidates[:10],  # First 10 as suggestions
    )


# ─────────────────────────────────────────────────────────────────────────────
# Planning Functions
# ─────────────────────────────────────────────────────────────────────────────

def create_upload_plan(
    scan_result: ScanResult,
    campaign_name: str,
    reference_crystal: Optional[str] = None,
    compound_pattern: Optional[str] = r'NCL-\d{8}',
    include_apo: bool = True,
) -> UploadPlan:
    """
    Create an upload plan from scan results.

    Args:
        scan_result: Output from scan stage
        campaign_name: Name for the campaign
        reference_crystal: Crystal to use for parent project reference files
        compound_pattern: Regex pattern to filter compound identifiers (None = all)
        include_apo: Include APO crystals with NCL-00000000
    """
    upload_date = datetime.now().strftime('%Y%m%d')

    plan = UploadPlan(
        plan_date=datetime.now().isoformat(),
        campaign_name=campaign_name,
        upload_date_tag=upload_date,
        parent_project_name=f'{campaign_name}_parent',
        compound_pattern=compound_pattern,
    )

    # Find reference crystal for parent project files
    if reference_crystal:
        ref_data = next(
            (c for c in scan_result.crystals if c.crystal_name == reference_crystal),
            None
        )
    else:
        # Auto-select first crystal with dimple output
        ref_data = next(
            (c for c in scan_result.crystals if c.has_dimple and c.dimple_pdb),
            None
        )

    if ref_data:
        plan.reference_pdb = ref_data.dimple_pdb
        # Use any free.mtz as the collective FreeR
        plan.reference_freer = ref_data.free_mtz

    # Compile pattern if provided
    pattern_re = re.compile(compound_pattern) if compound_pattern else None

    # Build actions for each crystal
    for crystal in scan_result.crystals:
        # Skip if no unmerged MTZ (required for SubstituteLigand)
        if not crystal.unmerged_mtz:
            continue

        # Determine compound ID
        if crystal.compound:
            compound_id = crystal.compound.identifier
        elif crystal.has_dimple:
            # APO crystal (dimple output but no compound)
            if not include_apo:
                continue
            compound_id = 'NCL-00000000'
        else:
            # No compound info at all
            compound_id = 'unknown'

        # Apply compound pattern filter (APO always passes if include_apo is True)
        if pattern_re and compound_id != 'NCL-00000000':
            if not pattern_re.match(compound_id):
                continue

        # Generate project name: {prefix}_x{index}_{compound}_{date}
        # Uses just the crystal index (x0002) to avoid redundancy with prefix
        project_name = f'{scan_result.target_prefix}_x{crystal.crystal_index:04d}_{compound_id}_{upload_date}'

        # Get compound CIF and SMILES if available (from local XChem data)
        compound_cif = crystal.compound.cif_path if crystal.compound else None
        smiles = crystal.compound.smiles if crystal.compound else None

        action = UploadAction(
            crystal_name=crystal.crystal_name,
            project_name=project_name,
            compound_id=compound_id,
            unmerged_mtz_path=crystal.unmerged_mtz,
            # Dimple harvest files
            dimple_pdb=crystal.dimple_pdb,
            dimple_mtz=crystal.dimple_mtz,
            compound_cif=compound_cif,
            # SMILES from local .smiles file (no API lookup needed)
            smiles=smiles,
        )
        plan.actions.append(action)

    # Update summary counts
    plan.total_actions = len(plan.actions)
    plan.ncl_actions = sum(1 for a in plan.actions if a.compound_id.startswith('NCL-') and a.compound_id != 'NCL-00000000')
    plan.apo_actions = sum(1 for a in plan.actions if a.compound_id == 'NCL-00000000')

    return plan


# ─────────────────────────────────────────────────────────────────────────────
# YAML I/O
# ─────────────────────────────────────────────────────────────────────────────

def dataclass_to_dict(obj) -> Any:
    """Recursively convert dataclass instances to dicts for YAML serialization."""
    if hasattr(obj, '__dataclass_fields__'):
        return {k: dataclass_to_dict(v) for k, v in asdict(obj).items()}
    elif isinstance(obj, list):
        return [dataclass_to_dict(item) for item in obj]
    elif isinstance(obj, dict):
        return {k: dataclass_to_dict(v) for k, v in obj.items()}
    return obj


def save_config(data: Any, output_path: str):
    """Save data to YAML or JSON file (auto-detect from extension)."""
    data_dict = dataclass_to_dict(data)

    if output_path.endswith('.yaml') or output_path.endswith('.yml'):
        if not HAS_YAML:
            # Fallback to JSON
            output_path = output_path.rsplit('.', 1)[0] + '.json'
            print(f"Note: PyYAML not installed, using JSON format instead")

    if output_path.endswith('.json'):
        with open(output_path, 'w') as f:
            json.dump(data_dict, f, indent=2)
    else:
        # YAML format
        with open(output_path, 'w') as f:
            yaml.dump(data_dict, f, default_flow_style=False, sort_keys=False)

    print(f"Saved: {output_path}")


def load_config(input_path: str) -> dict:
    """Load data from YAML or JSON file (auto-detect from extension)."""
    with open(input_path) as f:
        if input_path.endswith('.json'):
            return json.load(f)
        elif HAS_YAML:
            return yaml.safe_load(f)
        else:
            raise ImportError(
                f"Cannot load YAML file without PyYAML. "
                f"Install with: pip install pyyaml"
            )


# ─────────────────────────────────────────────────────────────────────────────
# CLI Commands
# ─────────────────────────────────────────────────────────────────────────────

def cmd_scan(args):
    """Stage 1: Scan local data directory."""
    print(f"Scanning: {args.data_dir}")

    result = scan_data_directory(args.data_dir, prefix_filter=args.prefix)

    print(f"\n=== Scan Summary ===")
    print(f"Target: {result.target_prefix}")
    print(f"Total crystals: {result.total_crystals}")
    print(f"  NCL-identified (ready): {result.ncl_identified}")
    print(f"  Other compounds (need resolution): {result.other_compounds}")
    print(f"  APO (no compound): {result.apo_crystals}")
    print(f"  Incomplete (missing MTZ): {result.incomplete}")

    if result.reference_candidates:
        print(f"\nReference candidates (have dimple output):")
        for name in result.reference_candidates[:5]:
            print(f"  - {name}")

    output = args.output or 'campaign_scan.yaml'
    save_config(result, output)

    print(f"\nNext step: Review {output}, then run:")
    print(f"  campaign_uploader plan {output} -o campaign_plan.yaml")


def cmd_plan(args):
    """Stage 2: Create upload plan from scan results."""
    print(f"Loading scan: {args.scan_file}")

    scan_data = load_config(args.scan_file)

    # Reconstruct ScanResult (simplified - just use dict)
    campaign_name = args.campaign or scan_data.get('target_prefix', 'Campaign')

    # Build crystals list from scan data
    crystals = []
    for c in scan_data.get('crystals', []):
        compound = None
        if c.get('compound'):
            compound = CompoundInfo(**c['compound'])
        crystals.append(CrystalData(
            crystal_name=c['crystal_name'],
            crystal_index=c['crystal_index'],
            dir_path=c['dir_path'],
            merged_mtz=c.get('merged_mtz'),
            free_mtz=c.get('free_mtz'),
            unmerged_mtz=c.get('unmerged_mtz'),
            dimple_pdb=c.get('dimple_pdb'),
            dimple_mtz=c.get('dimple_mtz'),
            compound=compound,
            has_dimple=c.get('has_dimple', False),
            is_complete=c.get('is_complete', False),
            is_ncl_identified=c.get('is_ncl_identified', False),
        ))

    scan_result = ScanResult(
        scan_date=scan_data['scan_date'],
        source_dir=scan_data['source_dir'],
        target_prefix=scan_data['target_prefix'],
        total_crystals=scan_data['total_crystals'],
        ncl_identified=scan_data['ncl_identified'],
        other_compounds=scan_data['other_compounds'],
        apo_crystals=scan_data['apo_crystals'],
        incomplete=scan_data['incomplete'],
        crystals=crystals,
        reference_candidates=scan_data.get('reference_candidates', []),
    )

    # Determine compound pattern
    if args.all_compounds:
        compound_pattern = None  # No filter
    elif args.compound_pattern:
        compound_pattern = args.compound_pattern
    else:
        compound_pattern = r'NCL-\d{8}'  # Default: NCL only

    plan = create_upload_plan(
        scan_result,
        campaign_name=campaign_name,
        reference_crystal=args.reference,
        compound_pattern=compound_pattern,
        include_apo=args.include_apo,
    )

    print(f"\n=== Upload Plan ===")
    print(f"Campaign: {plan.campaign_name}")
    print(f"Upload date tag: {plan.upload_date_tag}")
    print(f"Compound filter: {plan.compound_pattern or '(all compounds)'}")
    print(f"Total actions: {plan.total_actions}")
    print(f"  NCL compounds: {plan.ncl_actions}")
    print(f"  APO crystals: {plan.apo_actions}")

    if plan.reference_pdb:
        print(f"\nParent project reference:")
        print(f"  PDB: {plan.reference_pdb}")
        print(f"  FreeR: {plan.reference_freer}")

    output = args.output or 'campaign_plan.yaml'
    save_config(plan, output)

    print(f"\nNext step: Review {output}, then run:")
    print(f"  campaign_uploader execute {output} --dry-run  # Preview")
    print(f"  campaign_uploader execute {output}            # Execute")


def get_i2remote_client():
    """
    Import and instantiate the CCP4i2Client from i2remote.

    Returns:
        CCP4i2Client instance configured from ~/.ccp4i2remote.json
    """
    # Try to import from the Docker/cli directory
    import sys
    cli_path = Path(__file__).parent.parent.parent.parent / 'cli'
    if str(cli_path) not in sys.path:
        sys.path.insert(0, str(cli_path))

    try:
        from i2remote import CCP4i2Client, load_config as load_i2config
        config = load_i2config()
        return CCP4i2Client(
            api_url=config['api_url'],
            token=config.get('api_token'),
            timeout=config.get('timeout', 30),
            verify_ssl=config.get('verify_ssl', True),
        )
    except ImportError as e:
        raise ImportError(
            f"Could not import i2remote: {e}\n"
            "Ensure Docker/cli/i2remote.py is accessible."
        )


def lookup_smiles(client, ncl_id: str) -> Optional[str]:
    """
    Look up SMILES for an NCL compound from the compounds registry.

    The compounds registry API is at a different base URL than the ccp4i2 API:
    - ccp4i2:    https://ddudatabase.ncl.ac.uk/api/proxy/ccp4i2
    - compounds: https://ddudatabase.ncl.ac.uk/api/proxy/compounds

    This function derives the compounds URL from the ccp4i2 URL and uses
    the same authentication token.

    Args:
        client: CCP4i2Client instance (authenticated with Azure AD)
        ncl_id: NCL identifier (e.g., 'NCL-00031069')

    Returns:
        SMILES string if found, None otherwise
    """
    import requests

    # Extract numeric part from NCL-XXXXXXXX
    match = re.match(r'NCL-(\d+)', ncl_id)
    if not match:
        return None

    reg_id = match.group(1).lstrip('0') or '0'  # Remove leading zeros

    try:
        # Derive compounds API URL from ccp4i2 URL
        # e.g., https://ddudatabase.ncl.ac.uk/api/proxy/ccp4i2 -> .../compounds
        compounds_url = client.api_url.replace('/ccp4i2', '/compounds')
        endpoint = f"{compounds_url}/compounds/?reg_number={reg_id}"

        # Use the same authentication token
        auth_header = client.session.headers.get('Authorization', '')
        headers = {'Authorization': auth_header} if auth_header else {}

        response = requests.get(
            endpoint,
            headers=headers,
            timeout=client.timeout,
            verify=client.verify_ssl,
        )

        if response.ok:
            data = response.json()
            # DRF returns list for filter queries
            if isinstance(data, list) and len(data) > 0:
                return data[0].get('smiles')
            # Or could be paginated results
            elif isinstance(data, dict) and data.get('results'):
                results = data['results']
                if len(results) > 0:
                    return results[0].get('smiles')
        else:
            print(f"    Warning: Compounds API returned {response.status_code} for {ncl_id}")

    except Exception as e:
        print(f"    Warning: Could not look up SMILES for {ncl_id}: {e}")

    return None


def resolve_compound_by_smiles(client, smiles: str) -> Optional[dict]:
    """
    Resolve a non-NCL compound by SMILES via canonical SMILES matching in the registry.

    Uses the compounds registry's resolve_by_smiles endpoint to find a
    registered compound that matches the given SMILES. The server converts
    SMILES to InChI and performs an exact lookup.

    Args:
        client: CCP4i2Client instance (authenticated with Azure AD)
        smiles: SMILES string from local XChem data

    Returns:
        Dict with 'smiles', 'reg_number', 'formatted_id' if found, None otherwise.
        The returned SMILES is the registry's canonical SMILES.
    """
    import requests
    from urllib.parse import quote

    try:
        # Derive compounds API URL from ccp4i2 URL
        compounds_url = client.api_url.replace('/ccp4i2', '/compounds')
        endpoint = f"{compounds_url}/compounds/resolve_by_smiles/?smiles={quote(smiles)}"

        # Use the same authentication token
        auth_header = client.session.headers.get('Authorization', '')
        headers = {'Authorization': auth_header} if auth_header else {}

        response = requests.get(
            endpoint,
            headers=headers,
            timeout=client.timeout,
            verify=client.verify_ssl,
        )

        if not response.ok:
            print(f"    Resolution endpoint returned HTTP {response.status_code}")
            try:
                print(f"    Response: {response.json()}")
            except Exception:
                print(f"    Response: {response.text[:200]}")
            return None

        data = response.json()
        if data.get('found') and data.get('compound'):
            compound = data['compound']
            return {
                'smiles': compound.get('smiles'),
                'reg_number': compound.get('reg_number'),
                'formatted_id': compound.get('formatted_id'),
            }
        else:
            print(f"    No match. query_inchi={data.get('query_inchi', 'N/A')}")
            return None

    except Exception as e:
        print(f"    Warning: Could not resolve compound by SMILES: {e}")

    return None


def find_project_by_name(client, project_name: str) -> Optional[dict]:
    """
    Find an existing project by exact name match.

    Args:
        client: CCP4i2Client instance
        project_name: Exact project name to search for

    Returns:
        Project dict if found, None otherwise
    """
    try:
        projects = client.list_projects()
        for project in projects:
            if project.get('name') == project_name:
                return project
    except Exception as e:
        print(f"    Warning: Could not list projects: {e}")
    return None


def execute_upload(
    client,
    plan_data: dict,
    campaign_id: int,
    coords_file_id: int,
    dry_run: bool = False,
    verbose: bool = False,
    status_file: Optional[str] = None,
    remaining_file: Optional[str] = None,
    force_smiles_fallback: bool = False,
    rerun: bool = False,
    rerun_create_missing: bool = False,
    run_sublig: bool = True,
    run_servalcat: bool = True,
) -> dict:
    """
    Execute the upload plan using i2remote client.

    This follows the same pattern as BatchImportDialog:
    1. Create project and add to campaign (or find existing in rerun mode)
    2. Create SubstituteLigand job
    3. Upload reference coords (XYZIN) from parent
    4. Set SMILES parameter (or LIGANDAS=NONE for APO)
    5. Set PIPELINE to DIMPLE
    6. Upload unmerged data to UNMERGEDFILES[0].file
    7. Run the job

    Args:
        client: CCP4i2Client instance
        plan_data: Loaded plan dictionary
        campaign_id: ID of pre-created campaign
        coords_file_id: File ID of reference coordinates from parent project
        dry_run: If True, only show what would be done
        verbose: If True (with dry_run), print i2remote commands
        status_file: Path to write status updates
        remaining_file: Path to save unresolved non-NCL compounds for later
            registration. If None, unresolved compounds are only reported.
        force_smiles_fallback: If True, fall back to local .smiles files when
            registry lookup fails. If False (default), fail the action instead
            to avoid potential SMILES mismatches.
        rerun: If True, find existing projects by name instead of creating new ones.
            Creates new jobs in existing projects (old jobs are preserved).
        rerun_create_missing: If True (with rerun), create projects that don't exist.
            If False (default with rerun), skip projects that don't exist.
        run_sublig: If True (default), create SubstituteLigand jobs.
        run_servalcat: If True (default), create servalcat_pipe jobs.

    Returns:
        Status dictionary with results for each action
    """
    results = {
        'campaign_id': campaign_id,
        'coords_file_id': coords_file_id,
        'plan_date': plan_data['plan_date'],
        'execution_start': datetime.now().isoformat(),
        'actions': [],
        'summary': {
            'total': 0,
            'completed': 0,
            'failed': 0,
            'skipped': 0,
        }
    }

    unresolved_compounds = []
    actions = plan_data.get('actions', [])
    results['summary']['total'] = len(actions)

    # Download reference coords once (reused for all jobs)
    coords_content = None
    coords_filename = 'reference.pdb'
    if not dry_run:
        print(f"\nDownloading reference coordinates (file ID: {coords_file_id})...")
        try:
            coords_content = client.download_file(coords_file_id)
            # Try to get filename from file metadata
            file_info = client.get(f"/files/{coords_file_id}")
            coords_filename = file_info.get('filename', 'reference.pdb')
            print(f"  Downloaded: {coords_filename} ({len(coords_content)} bytes)")
        except Exception as e:
            print(f"  ERROR: Failed to download coords: {e}")
            results['error'] = f"Failed to download reference coords: {e}"
            return results

    print(f"\nProcessing {len(actions)} actions...")

    for i, action in enumerate(actions):
        print(f"\n[{i+1}/{len(actions)}] {action['crystal_name']}")
        print(f"    Project: {action['project_name']}")
        print(f"    Compound: {action['compound_id']}")

        action_result = {
            'crystal_name': action['crystal_name'],
            'project_name': action['project_name'],
            'compound_id': action['compound_id'],
            'resolved_compound_id': None,
            'status': 'pending',
            'project_id': None,
            'sublig_job_id': None,
            'servalcat_pipe_job_id': None,
            'error': None,
        }

        if dry_run:
            if verbose:
                # Print verbose i2remote commands that would be executed
                print(f"    # === i2remote commands for {action['crystal_name']} ===")
                print(f"    client.create_project('{action['project_name']}')")
                print(f"    client.post('/projectgroups/{campaign_id}/add_member', {{'project_id': <project_id>, 'type': 'member'}})")

                # SubstituteLigand job commands
                if run_sublig:
                    print(f"    # --- SubstituteLigand job ---")
                    print(f"    client.create_job(<project_id>, 'SubstituteLigand')")
                    print(f"    client.upload_file_content(<job_id>, 'SubstituteLigand.inputData.XYZIN', <coords_content>, 'reference.pdb')")

                    compound_id = action['compound_id']
                    if compound_id == 'NCL-00000000':
                        print(f"    client.set_job_parameter(<job_id>, 'SubstituteLigand.controlParameters.LIGANDAS', 'NONE')")
                    elif compound_id.startswith('NCL-'):
                        print(f"    smiles = lookup_smiles(client, '{compound_id}')  # Registry lookup")
                        print(f"    client.set_job_parameter(<job_id>, 'SubstituteLigand.inputData.SMILESIN', smiles)")
                    else:
                        local_smiles = action.get('smiles', '')
                        print(f"    # Resolve non-NCL compound {compound_id} via canonical SMILES")
                        print(f"    resolved = resolve_compound_by_smiles(client, '{local_smiles[:50]}...')")
                        print(f"    # If resolved: use registry SMILES")
                        print(f"    # If not resolved: SKIP (saved to --remaining file)")

                    print(f"    client.set_job_parameter(<job_id>, 'SubstituteLigand.inputData.PIPELINE', 'DIMPLE')")
                    print(f"    client.upload_file_param(<job_id>, 'SubstituteLigand.inputData.UNMERGEDFILES[0].file', '{action['unmerged_mtz_path']}')")
                    print(f"    client.run_job(<job_id>)")
                else:
                    print(f"    # --- SubstituteLigand skipped (--only-servalcat) ---")

                # Servalcat_pipe job commands
                if run_servalcat:
                    dimple_pdb = action.get('dimple_pdb')
                    dimple_mtz = action.get('dimple_mtz')
                    compound_cif = action.get('compound_cif')

                    if dimple_pdb and dimple_mtz:
                        print(f"    # --- servalcat_pipe job (dimple harvest) ---")
                        print(f"    client.create_job(<project_id>, 'servalcat_pipe')")
                        print(f"    client.upload_file_param(<servalcat_pipe_job_id>, 'servalcat_pipe.inputData.HKLIN', '{dimple_mtz}')")
                        print(f"    client.upload_file_param(<servalcat_pipe_job_id>, 'servalcat_pipe.inputData.XYZIN', '{dimple_pdb}')")
                        if compound_cif:
                            print(f"    client.upload_file_param(<servalcat_pipe_job_id>, 'servalcat_pipe.inputData.DICT_LIST[0]', '{compound_cif}')")
                        print(f"    client.run_job(<servalcat_pipe_job_id>)")
                    else:
                        print(f"    # --- servalcat_pipe skipped (no dimple output) ---")
                else:
                    print(f"    # --- servalcat_pipe skipped (--only-sublig) ---")
                print()
            else:
                print(f"    [DRY RUN] Would create project and upload files")

            action_result['status'] = 'skipped'
            results['summary']['skipped'] += 1
            results['actions'].append(action_result)
            continue

        try:
            # Step 0: Resolve non-NCL compounds via canonical SMILES before creating anything
            compound_id = action['compound_id']
            resolved_smiles = None  # Will hold registry SMILES if resolved

            if not compound_id.startswith('NCL-'):
                local_smiles = action.get('smiles')
                if not local_smiles:
                    # No SMILES at all - can't resolve
                    print(f"    No SMILES available for {compound_id} - skipping")
                    unresolved_compounds.append({
                        'crystal_name': action['crystal_name'],
                        'compound_id': compound_id,
                        'smiles': None,
                        'compound_cif': action.get('compound_cif'),
                        'unmerged_mtz_path': action.get('unmerged_mtz_path'),
                        'dimple_pdb': action.get('dimple_pdb'),
                        'dimple_mtz': action.get('dimple_mtz'),
                    })
                    action_result['status'] = 'skipped'
                    action_result['error'] = f'No SMILES available for {compound_id}'
                    results['actions'].append(action_result)
                    results['summary']['skipped'] += 1
                    continue
                else:
                    print(f"    Resolving {compound_id} via canonical SMILES lookup...")
                    resolved = resolve_compound_by_smiles(client, local_smiles)

                    if resolved:
                        resolved_id = resolved['formatted_id']
                        resolved_smiles = resolved['smiles']
                        action_result['resolved_compound_id'] = resolved_id
                        # Update project name with resolved NCL ID
                        action['project_name'] = action['project_name'].replace(
                            compound_id, resolved_id
                        )
                        print(f"    Resolved to {resolved_id}")
                    else:
                        # Not in registry - skip and save for later registration
                        print(f"    Not resolved - skipping (will be saved to remaining file)")
                        unresolved_compounds.append({
                            'crystal_name': action['crystal_name'],
                            'compound_id': compound_id,
                            'smiles': local_smiles,
                            'compound_cif': action.get('compound_cif'),
                            'unmerged_mtz_path': action.get('unmerged_mtz_path'),
                            'dimple_pdb': action.get('dimple_pdb'),
                            'dimple_mtz': action.get('dimple_mtz'),
                        })
                        action_result['status'] = 'skipped'
                        action_result['error'] = f'Compound {compound_id} not found in registry'
                        results['actions'].append(action_result)
                        results['summary']['skipped'] += 1
                        continue

            # Step 1: Create or find project
            if rerun:
                # Rerun mode: look for existing project by name
                print(f"    Looking for existing project...")
                existing_project = find_project_by_name(client, action['project_name'])

                if existing_project:
                    project_id = existing_project.get('id')
                    action_result['project_id'] = project_id
                    print(f"    Found existing project ID: {project_id}")
                    # Skip adding to campaign - already a member
                elif rerun_create_missing:
                    # Project doesn't exist, but we're allowed to create it
                    print(f"    Project not found, creating...")
                    project_resp = client.create_project(action['project_name'])
                    project_id = project_resp.get('id')
                    action_result['project_id'] = project_id
                    print(f"    Created project ID: {project_id}")

                    # Add new project to campaign
                    print(f"    Adding to campaign {campaign_id}...")
                    client.post(f"/projectgroups/{campaign_id}/add_member", {
                        'project_id': project_id,
                        'type': 'member'
                    })
                else:
                    # Project doesn't exist and we're not creating - skip
                    print(f"    Project not found, skipping (use --rerun-create-missing to create)")
                    action_result['status'] = 'skipped'
                    action_result['error'] = 'Project not found in rerun mode'
                    results['actions'].append(action_result)
                    results['summary']['skipped'] += 1
                    continue
            else:
                # Normal mode: create new project
                print(f"    Creating project...")
                project_resp = client.create_project(action['project_name'])
                project_id = project_resp.get('id')
                action_result['project_id'] = project_id
                print(f"    Created project ID: {project_id}")

                # Step 2: Add project to campaign as member
                print(f"    Adding to campaign {campaign_id}...")
                client.post(f"/projectgroups/{campaign_id}/add_member", {
                    'project_id': project_id,
                    'type': 'member'
                })

            # Step 3: Create SubstituteLigand job (if enabled)
            if run_sublig:
                print(f"    Creating SubstituteLigand job...")
                sublig_job = client.create_job(project_id, 'SubstituteLigand', auto_context=False)
                # i2remote.create_job returns unwrapped job data: {id, uuid, ...}
                sublig_job_id = sublig_job.get('id')
                print(f"    Created SubstituteLigand job ID: {sublig_job_id}")

                # Step 4: Upload reference coordinates (XYZIN)
                print(f"    Uploading reference coords...")
                client.upload_file_content(
                    sublig_job_id,
                    'SubstituteLigand.inputData.XYZIN',
                    coords_content,
                    coords_filename
                )

                # Step 5: Set compound parameters
                compound_id = action['compound_id']

                if compound_id == 'NCL-00000000':
                    # APO crystal - no ligand
                    print(f"    Setting LIGANDAS=NONE (APO)...")
                    client.set_job_parameter(
                        sublig_job_id,
                        'SubstituteLigand.controlParameters.LIGANDAS',
                        'NONE'
                    )
                elif compound_id.startswith('NCL-'):
                    # Look up SMILES from the compounds registry (registerer's preferred SMILES)
                    print(f"    Looking up SMILES for {compound_id}...")
                    smiles = lookup_smiles(client, compound_id)

                    # Handle registry lookup failure
                    if not smiles:
                        local_smiles = action.get('smiles')  # From local XChem data
                        if force_smiles_fallback and local_smiles:
                            print(f"    WARNING: Registry lookup failed, using local SMILES (--force-smiles-fallback)")
                            print(f"    Local SMILES: {local_smiles[:50]}...")
                            smiles = local_smiles
                        elif local_smiles:
                            # Local SMILES exists but fallback not enabled - fail safely
                            error_msg = (
                                f"Registry lookup failed for {compound_id}. "
                                f"Local SMILES exists but --force-smiles-fallback not set. "
                                f"Skipping to avoid potential SMILES mismatch."
                            )
                            print(f"    ERROR: {error_msg}")
                            action_result['status'] = 'failed'
                            action_result['error'] = error_msg
                            results['actions'].append(action_result)
                            results['summary']['failed'] += 1
                            continue
                        else:
                            print(f"    Warning: No SMILES available for {compound_id} (neither registry nor local)")

                    if smiles:
                        print(f"    Setting SMILES: {smiles[:50]}...")
                        client.set_job_parameter(
                            sublig_job_id,
                            'SubstituteLigand.inputData.SMILESIN',
                            smiles
                        )
                    else:
                        print(f"    Warning: No SMILES available for {compound_id}")
                else:
                    # Non-NCL compound (Z*, POB*, etc.) - resolved earlier via canonical SMILES
                    smiles = resolved_smiles
                    if smiles:
                        print(f"    Setting SMILES (from registry): {smiles[:50]}...")
                        client.set_job_parameter(
                            sublig_job_id,
                            'SubstituteLigand.inputData.SMILESIN',
                            smiles
                        )
                    else:
                        print(f"    Warning: No resolved SMILES for {compound_id}")

                # Step 6: Set pipeline to DIMPLE
                print(f"    Setting PIPELINE=DIMPLE...")
                client.set_job_parameter(
                    sublig_job_id,
                    'SubstituteLigand.inputData.PIPELINE',
                    'DIMPLE'
                )

                # Step 7: Upload unmerged reflection data
                unmerged_path = action['unmerged_mtz_path']
                print(f"    Uploading unmerged data: {Path(unmerged_path).name}...")
                client.upload_file_param(
                    sublig_job_id,
                    'SubstituteLigand.inputData.UNMERGEDFILES[0].file',
                    unmerged_path
                )

                # Step 8: Queue SubstituteLigand job
                print(f"    Queueing SubstituteLigand job...")
                client.run_job(sublig_job_id)
                action_result['sublig_job_id'] = sublig_job_id
            else:
                print(f"    Skipping SubstituteLigand (--only-servalcat)")

            # ─────────────────────────────────────────────────────────────
            # Create servalcat_pipe job to harvest dimple results
            # ─────────────────────────────────────────────────────────────
            if run_servalcat:
                dimple_pdb = action.get('dimple_pdb')
                dimple_mtz = action.get('dimple_mtz')

                if dimple_pdb and dimple_mtz and Path(dimple_pdb).exists() and Path(dimple_mtz).exists():
                    print(f"    Creating servalcat_pipe job (dimple harvest)...")
                    servalcat_pipe_job = client.create_job(project_id, 'servalcat_pipe', auto_context=False)
                    # i2remote.create_job returns unwrapped job data: {id, uuid, ...}
                    servalcat_pipe_job_id = servalcat_pipe_job.get('id')
                    action_result['servalcat_pipe_job_id'] = servalcat_pipe_job_id
                    print(f"    Created servalcat_pipe job ID: {servalcat_pipe_job_id}")

                    # Upload dimple MTZ as HKLIN
                    print(f"    Uploading HKLIN: {Path(dimple_mtz).name}...")
                    client.upload_file_param(
                        servalcat_pipe_job_id,
                        'servalcat_pipe.inputData.HKLIN',
                        dimple_mtz
                    )

                    # Upload dimple PDB as XYZIN
                    print(f"    Uploading XYZIN: {Path(dimple_pdb).name}...")
                    client.upload_file_param(
                        servalcat_pipe_job_id,
                        'servalcat_pipe.inputData.XYZIN',
                        dimple_pdb
                    )

                    # Upload compound CIF as DICT_LIST[0] if available
                    compound_cif = action.get('compound_cif')
                    if compound_cif and Path(compound_cif).exists():
                        print(f"    Uploading DICT: {Path(compound_cif).name}...")
                        client.upload_file_param(
                            servalcat_pipe_job_id,
                            'servalcat_pipe.inputData.DICT_LIST[0]',
                            compound_cif
                        )

                    # Queue servalcat_pipe job
                    print(f"    Queueing servalcat_pipe job...")
                    client.run_job(servalcat_pipe_job_id)
                else:
                    print(f"    Skipping servalcat_pipe (no dimple output)")
            else:
                print(f"    Skipping servalcat_pipe (--only-sublig)")

            action_result['status'] = 'completed'
            results['summary']['completed'] += 1
            print(f"    Done!")

        except Exception as e:
            action_result['status'] = 'failed'
            action_result['error'] = str(e)
            results['summary']['failed'] += 1
            print(f"    ERROR: {e}")

        results['actions'].append(action_result)

        # Save intermediate status
        if status_file:
            save_config(results, status_file)

    # Add resolution summary
    results['summary']['resolved'] = sum(
        1 for a in results['actions']
        if a.get('resolved_compound_id')
    )
    results['summary']['unresolved'] = len(unresolved_compounds)

    # Save unresolved compounds to remaining file
    if unresolved_compounds and remaining_file:
        remaining_data = {
            'campaign_name': plan_data.get('campaign_name'),
            'execution_date': datetime.now().isoformat(),
            'total_unresolved': len(unresolved_compounds),
            'note': 'Register these compounds in the registry, then re-run with --rerun',
            'compounds': unresolved_compounds,
        }
        save_config(remaining_data, remaining_file)
        print(f"\nUnresolved compounds saved to: {remaining_file}")
    elif unresolved_compounds:
        print(f"\nWarning: {len(unresolved_compounds)} unresolved compounds "
              f"(use --remaining to save them)")

    results['execution_end'] = datetime.now().isoformat()
    return results


def cmd_execute(args):
    """Stage 3: Execute upload plan using i2remote."""
    print(f"Loading plan: {args.plan_file}")

    plan_data = load_config(args.plan_file)

    # Validate pipeline selection flags
    if args.only_sublig and args.only_servalcat:
        print("\nERROR: Cannot use --only-sublig and --only-servalcat together.")
        return 1

    # Determine which pipelines to run
    run_sublig = not args.only_servalcat
    run_servalcat = not args.only_sublig

    if args.only_sublig:
        print("\n[MODE] Running only SubstituteLigand (skipping servalcat_pipe)")
    elif args.only_servalcat:
        print("\n[MODE] Running only servalcat_pipe (skipping SubstituteLigand)")

    if args.dry_run:
        print("\n=== DRY RUN - No changes will be made ===")
        if args.verbose:
            print("=== VERBOSE MODE - Showing i2remote commands ===")
        print()

    print(f"Campaign: {plan_data['campaign_name']}")
    print(f"Actions to execute: {plan_data['total_actions']}")

    # Validate campaign ID is provided (not required for dry-run)
    if not args.campaign_id and not args.dry_run:
        print("\nERROR: Campaign ID required.")
        print("Create a campaign first via the web UI, then provide its ID:")
        print(f"  campaign_uploader execute {args.plan_file} --campaign-id <ID>")
        print("\nTo find campaign IDs:")
        print("  1. Open the CCP4i2 web app")
        print("  2. Go to Campaigns")
        print("  3. Create or select a campaign")
        print("  4. The ID is shown in the URL or campaign details")
        return 1

    # Preview actions
    print("\n--- Actions Preview ---")
    for i, action in enumerate(plan_data.get('actions', [])[:5]):
        print(f"\n[{i+1}] {action['crystal_name']}")
        print(f"    Project: {action['project_name']}")
        print(f"    Compound: {action['compound_id']}")
        print(f"    Unmerged MTZ: {Path(action['unmerged_mtz_path']).name}")

    if len(plan_data.get('actions', [])) > 5:
        print(f"\n... and {len(plan_data['actions']) - 5} more actions")

    if args.dry_run:
        if args.verbose:
            # In verbose dry-run, show detailed commands for all actions
            print("\n--- i2remote Commands Preview ---")
            campaign_id = args.campaign_id or '<campaign_id>'
            coords_file_id = '<coords_file_id>'
            results = execute_upload(
                client=None,
                plan_data=plan_data,
                campaign_id=campaign_id,
                coords_file_id=coords_file_id,
                dry_run=True,
                verbose=True,
                run_sublig=run_sublig,
                run_servalcat=run_servalcat,
            )
        print("\n[DRY RUN] No changes made.")
        return 0

    # Execute the upload
    try:
        print("\n--- Initializing i2remote client ---")
        client = get_i2remote_client()

        # Verify campaign exists (also verifies authentication)
        campaign = client.get(f"/projectgroups/{args.campaign_id}")
        print(f"Campaign: {campaign.get('name')} (ID: {args.campaign_id})")

        # Get parent project files (for reference coordinates)
        print(f"Getting parent project files...")
        parent_files = client.get(f"/projectgroups/{args.campaign_id}/parent_files")
        coords_files = parent_files.get('coordinates', [])

        if not coords_files:
            print("\nERROR: No reference coordinates found in parent project.")
            print("Upload coordinates to the parent project first, then re-run.")
            return 1

        # Use the most recent coords file
        coords_file = coords_files[-1]  # Assuming sorted by date
        coords_file_id = coords_file.get('id')
        print(f"Using reference coords: {coords_file.get('filename')} (ID: {coords_file_id})")

    except Exception as e:
        print(f"\nERROR: Failed to connect: {e}")
        print("\nMake sure you're logged in:")
        print("  cd Docker/cli && python i2remote.py login")
        return 1

    # Confirm execution
    if not args.yes:
        if args.rerun:
            print(f"\n[RERUN MODE] This will create new jobs in {plan_data['total_actions']} existing projects.")
            if args.rerun_create_missing:
                print("Projects that don't exist will be created.")
            else:
                print("Projects that don't exist will be skipped.")
        else:
            print(f"\nThis will create {plan_data['total_actions']} projects and jobs.")
        response = input("Continue? [y/N]: ")
        if response.lower() != 'y':
            print("Aborted.")
            return 0

    # Execute
    status_file = args.output or 'campaign_status.json'
    results = execute_upload(
        client,
        plan_data,
        campaign_id=args.campaign_id,
        coords_file_id=coords_file_id,
        dry_run=False,
        status_file=status_file,
        remaining_file=args.remaining,
        force_smiles_fallback=args.force_smiles_fallback,
        rerun=args.rerun,
        rerun_create_missing=args.rerun_create_missing,
        run_sublig=run_sublig,
        run_servalcat=run_servalcat,
    )

    # Final summary
    print(f"\n=== Execution Complete ===")
    print(f"Total: {results['summary']['total']}")
    print(f"Completed: {results['summary']['completed']}")
    print(f"Failed: {results['summary']['failed']}")
    if results['summary'].get('skipped', 0) > 0:
        print(f"Skipped: {results['summary']['skipped']}")
    if results['summary'].get('resolved', 0) > 0:
        print(f"Resolved (InChI match): {results['summary']['resolved']}")
    if results['summary'].get('unresolved', 0) > 0:
        print(f"Unresolved (remaining): {results['summary']['unresolved']}")

    save_config(results, status_file)
    print(f"\nStatus saved to: {status_file}")

    return 0 if results['summary']['failed'] == 0 else 1


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        prog='campaign_uploader',
        description='Upload XChem fragment campaign data via i2remote',
    )

    subparsers = parser.add_subparsers(dest='command', help='Command')

    # Scan command
    scan_parser = subparsers.add_parser('scan', help='Stage 1: Scan local data')
    scan_parser.add_argument('data_dir', help='Path to local_data directory')
    scan_parser.add_argument('-o', '--output', help='Output YAML file')
    scan_parser.add_argument('-p', '--prefix', help='Filter by crystal prefix')

    # Plan command
    plan_parser = subparsers.add_parser('plan', help='Stage 2: Create upload plan')
    plan_parser.add_argument('scan_file', help='Scan YAML file from stage 1')
    plan_parser.add_argument('-o', '--output', help='Output YAML file')
    plan_parser.add_argument('-c', '--campaign', help='Campaign name')
    plan_parser.add_argument('-r', '--reference', help='Reference crystal name')
    plan_parser.add_argument('-p', '--compound-pattern',
                            help='Regex pattern for compound IDs to include (default: NCL-\\d{8})')
    plan_parser.add_argument('--all-compounds', action='store_true',
                            help='Include all compounds (no pattern filter)')
    plan_parser.add_argument('--include-apo', action='store_true', default=True,
                            help='Include APO crystals (default: True)')
    plan_parser.add_argument('--no-apo', dest='include_apo', action='store_false',
                            help='Exclude APO crystals')

    # Execute command
    exec_parser = subparsers.add_parser('execute', help='Stage 3: Execute upload')
    exec_parser.add_argument('plan_file', help='Plan YAML file from stage 2')
    exec_parser.add_argument('--campaign-id', type=int, required=False,
                            help='ID of pre-created campaign (required unless --dry-run)')
    exec_parser.add_argument('-n', '--dry-run', action='store_true',
                            help='Show what would be done without doing it')
    exec_parser.add_argument('-v', '--verbose', action='store_true',
                            help='With --dry-run, print detailed i2remote commands')
    exec_parser.add_argument('-y', '--yes', action='store_true',
                            help='Skip confirmation prompt')
    exec_parser.add_argument('-o', '--output', help='Output status JSON file')
    exec_parser.add_argument('--force-smiles-fallback', action='store_true',
                            help='Fall back to local .smiles files if registry lookup fails '
                                 '(by default, missing registry SMILES causes an error)')
    exec_parser.add_argument('--rerun', action='store_true',
                            help='Re-run mode: find existing projects by name instead of creating new ones. '
                                 'Creates new jobs in existing projects (old jobs preserved).')
    exec_parser.add_argument('--rerun-create-missing', action='store_true',
                            help='With --rerun: create projects that do not exist (default: skip missing)')
    exec_parser.add_argument('--only-sublig', action='store_true',
                            help='Run only SubstituteLigand pipeline (skip servalcat_pipe)')
    exec_parser.add_argument('--only-servalcat', action='store_true',
                            help='Run only servalcat_pipe pipeline (skip SubstituteLigand)')
    exec_parser.add_argument('--remaining',
                            help='Path to save unresolved non-NCL compounds for later registration '
                                 '(default: campaign_remaining.yaml)',
                            default='campaign_remaining.yaml')

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 1

    if args.command == 'scan':
        cmd_scan(args)
    elif args.command == 'plan':
        cmd_plan(args)
    elif args.command == 'execute':
        cmd_execute(args)

    return 0


if __name__ == '__main__':
    sys.exit(main())
