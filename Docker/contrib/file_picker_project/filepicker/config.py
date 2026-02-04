"""
Configuration and file template definitions for filepicker.

Defines what files to extract per crystal, including special handling
for symlink-dereferenced files and sibling file discovery.
"""

from typing import List, Dict, Any, Optional, NamedTuple
import os


class FileSpec(NamedTuple):
    """Specification for a file to extract from a crystal directory."""

    # Logical name for this file in the output
    output_name = None  # type: Optional[str]

    # Glob pattern relative to crystal directory
    pattern = None  # type: str

    # If True, follow symlinks and copy actual file content
    dereference = True  # type: bool

    # If True, also grab sibling files from the resolved symlink's directory
    # (used for getting unmerged data from the same autoprocessing dir)
    grab_siblings = None  # type: Optional[List[str]]

    # Whether this file is required (error if missing) or optional
    required = False  # type: bool


# Can't use NamedTuple with defaults in Python 3.6, so use a factory function
def file_spec(pattern,  # type: str
              output_name=None,  # type: Optional[str]
              dereference=True,  # type: bool
              grab_siblings=None,  # type: Optional[List[str]]
              required=False  # type: bool
              ):
    # type: (...) -> Dict[str, Any]
    """Create a file specification dictionary."""
    return {
        'pattern': pattern,
        'output_name': output_name,
        'dereference': dereference,
        'grab_siblings': grab_siblings,
        'required': required,
    }


# Default file template: the "key files" for crystallography data
# All files are optional - crystals without dimple results will just get fewer files
DEFAULT_FILE_TEMPLATE = [
    # Final refined model and data (optional - not all crystals have dimple results)
    file_spec(
        pattern='dimple/dimple/final.pdb',
        output_name='final.pdb',
        required=False,
    ),
    file_spec(
        pattern='dimple/dimple/final.mtz',
        output_name='final.mtz',
        required=False,
    ),

    # Compound information
    file_spec(
        pattern='compound/*.cif',
        output_name=None,  # Keep original name
        required=False,
    ),
    file_spec(
        pattern='compound/*.smiles',
        output_name=None,
        required=False,
    ),
    file_spec(
        pattern='compound/*.png',
        output_name=None,
        required=False,
    ),
    file_spec(
        pattern='compound/merged.cif',
        output_name='merged.cif',
        required=False,
    ),

    # The selected/best merged MTZ (top-level symlink)
    # Also grab unmerged data and xia2.mmcif.bz2 from the same autoprocessing dir
    file_spec(
        pattern='*.mtz',
        output_name=None,
        dereference=True,
        grab_siblings=['*_scaled_unmerged.mtz', 'xia2.mmcif.bz2'],
        required=False,
    ),

    # Electron density maps
    file_spec(
        pattern='2fofc.map',
        output_name='2fofc.map',
        dereference=True,
        required=False,
    ),
    file_spec(
        pattern='fofc.map',
        output_name='fofc.map',
        dereference=True,
        required=False,
    ),
]


# Minimal template - just the essentials (also optional)
MINIMAL_FILE_TEMPLATE = [
    file_spec(pattern='dimple/dimple/final.pdb', output_name='final.pdb', required=False),
    file_spec(pattern='dimple/dimple/final.mtz', output_name='final.mtz', required=False),
    file_spec(pattern='compound/*.cif', required=False),
    file_spec(pattern='compound/*.smiles', required=False),
]


# Strict template - requires dimple results (will skip crystals without them)
STRICT_FILE_TEMPLATE = [
    file_spec(
        pattern='dimple/dimple/final.pdb',
        output_name='final.pdb',
        required=True,
    ),
    file_spec(
        pattern='dimple/dimple/final.mtz',
        output_name='final.mtz',
        required=True,
    ),
    file_spec(pattern='compound/*.cif', required=False),
    file_spec(pattern='compound/*.smiles', required=False),
    file_spec(pattern='compound/*.png', required=False),
    file_spec(
        pattern='*.mtz',
        dereference=True,
        grab_siblings=['*_scaled_unmerged.mtz', 'xia2.mmcif.bz2'],
        required=True,
    ),
    file_spec(pattern='2fofc.map', dereference=True, required=False),
    file_spec(pattern='fofc.map', dereference=True, required=False),
]


# Full template - everything including intermediate files
FULL_FILE_TEMPLATE = DEFAULT_FILE_TEMPLATE + [
    # Additional dimple outputs
    file_spec(pattern='dimple/dimple/final.mmcif', required=False),
    file_spec(pattern='dimple/dimple/*.log', required=False),

    # Processing logs
    file_spec(pattern='*.log', dereference=True, required=False),
]


# Named template presets
TEMPLATES = {
    'default': DEFAULT_FILE_TEMPLATE,
    'minimal': MINIMAL_FILE_TEMPLATE,
    'full': FULL_FILE_TEMPLATE,
    'strict': STRICT_FILE_TEMPLATE,
}


def get_template(name):
    # type: (str) -> List[Dict[str, Any]]
    """Get a named file template."""
    if name not in TEMPLATES:
        raise ValueError("Unknown template '{}'. Available: {}".format(
            name, ', '.join(TEMPLATES.keys())))
    return TEMPLATES[name]


def parse_custom_template(spec_string):
    # type: (str) -> List[Dict[str, Any]]
    """
    Parse a custom template from a comma-separated string of patterns.

    Example: "final.pdb,compound/*.cif,*.mtz"
    """
    patterns = [p.strip() for p in spec_string.split(',') if p.strip()]
    return [file_spec(pattern=p) for p in patterns]
