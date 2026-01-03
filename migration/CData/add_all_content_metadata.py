#!/usr/bin/env python3
"""
Add CONTENT_FLAGS and SUBTYPES metadata for ALL CDataFile descendants.

Extracted from CCP4i2 source files:
- CCP4XtalData.py: CObsDataFile, CPhsDataFile, CMapCoeffsDataFile, CFreeRDataFile
- CCP4ModelData.py: CPdbDataFile
- CCP4CootData.py: CCootHistoryDataFile
"""

import json
from pathlib import Path

# Complete metadata for ALL classes with CONTENT_FLAGS/SUBTYPES

ALL_CLASS_METADATA = {
    # ========================================================================
    # MTZ File Classes (CCP4XtalData.py)
    # ========================================================================

    'CObsDataFile': {
        'SUBTYPES': {
            'SUBTYPE_OBSERVED': {
                'value': 1,
                'description': 'observed data'
            },
            'SUBTYPE_DERIVED': {
                'value': 2,
                'description': 'derived data'
            },
            'SUBTYPE_REFERENCE': {
                'value': 3,
                'description': 'reference data'
            }
        },
        'CONTENT_FLAGS': {
            'CONTENT_FLAG_IPAIR': {
                'value': 1,
                'annotation': 'Anomalous Is',
                'columns': ['Iplus', 'SIGIplus', 'Iminus', 'SIGIminus']
            },
            'CONTENT_FLAG_FPAIR': {
                'value': 2,
                'annotation': 'Anomalous SFs',
                'columns': ['Fplus', 'SIGFplus', 'Fminus', 'SIGFminus']
            },
            'CONTENT_FLAG_IMEAN': {
                'value': 3,
                'annotation': 'Mean Is',
                'columns': ['I', 'SIGI']
            },
            'CONTENT_FLAG_FMEAN': {
                'value': 4,
                'annotation': 'Mean SFs',
                'columns': ['F', 'SIGF']
            }
        }
    },

    'CPhsDataFile': {
        'SUBTYPES': {
            'SUBTYPE_UNBIASED': {
                'value': 1,
                'description': 'unbiased data'
            },
            'SUBTYPE_BIASED': {
                'value': 2,
                'description': 'biased data'
            }
        },
        'CONTENT_FLAGS': {
            'CONTENT_FLAG_HL': {
                'value': 1,
                'annotation': 'Hendrickson-Lattmann coeffs',
                'columns': ['HLA', 'HLB', 'HLC', 'HLD']
            },
            'CONTENT_FLAG_PHIFOM': {
                'value': 2,
                'annotation': 'Phi,FOM',
                'columns': ['PHI', 'FOM']
            }
        }
    },

    'CMapCoeffsDataFile': {
        'SUBTYPES': {
            'SUBTYPE_NORMAL': {
                'value': 1,
                'description': 'normal map'
            },
            'SUBTYPE_DIFFERENCE': {
                'value': 2,
                'description': 'difference map'
            },
            'SUBTYPE_ANOM_DIFFERENCE': {
                'value': 3,
                'description': 'anomalous difference map'
            }
        },
        'CONTENT_FLAGS': {
            'CONTENT_FLAG_FPHI': {
                'value': 1,
                'annotation': 'FPhi',
                'columns': ['F', 'PHI']
            }
        }
    },

    'CFreeRDataFile': {
        'SUBTYPES': {},  # No subtypes defined
        'CONTENT_FLAGS': {
            'CONTENT_FLAG_FREER': {
                'value': 1,
                'annotation': 'FreeR',
                'columns': ['FREER']
            }
        }
    },

    # ========================================================================
    # Model/Coordinate File Classes (CCP4ModelData.py)
    # ========================================================================

    'CPdbDataFile': {
        'SUBTYPES': {
            'SUBTYPE_UNKNOWN': {
                'value': 0,
                'description': 'unknown'
            },
            'SUBTYPE_MODEL': {
                'value': 1,
                'description': 'model'
            },
            'SUBTYPE_HOMOLOG': {
                'value': 2,
                'description': 'homolog'
            },
            'SUBTYPE_FRAGMENT': {
                'value': 3,
                'description': 'fragment'
            },
            'SUBTYPE_HEAVY_ATOMS': {
                'value': 4,
                'description': 'heavy atoms'
            }
        },
        'CONTENT_FLAGS': {
            'CONTENT_FLAG_PDB': {
                'value': 1,
                'annotation': 'PDB format',
                'columns': None  # Not applicable for non-MTZ files
            },
            'CONTENT_FLAG_MMCIF': {
                'value': 2,
                'annotation': 'mmCIF format',
                'columns': None  # Not applicable for non-MTZ files
            }
        }
    },

    # ========================================================================
    # Coot Data Classes (CCP4CootData.py)
    # ========================================================================

    'CCootHistoryDataFile': {
        'SUBTYPES': {
            'SUBTYPE_INITIAL': {
                'value': 1,
                'description': 'Coot 0-state.scm'
            },
            'SUBTYPE_HISTORY': {
                'value': 2,
                'description': 'Coot history.scm'
            }
        },
        'CONTENT_FLAGS': {}  # No content flags for Coot history files
    }
}


def add_all_metadata():
    """Add all CONTENT_FLAGS/SUBTYPES metadata to cdata.json."""

    cdata_path = Path(__file__).parent / 'cdata.json'

    print("="*70)
    print("Adding CONTENT_FLAGS and SUBTYPES metadata to cdata.json")
    print("="*70)

    print(f"\nReading {cdata_path}...")
    with open(cdata_path, 'r') as f:
        data = json.load(f)

    stats = {
        'found': 0,
        'not_found': 0,
        'subtypes_added': 0,
        'content_flags_added': 0
    }

    # Add metadata to each class
    for class_name, metadata in ALL_CLASS_METADATA.items():
        if class_name in data['classes']:
            stats['found'] += 1
            print(f"\n✓ {class_name} found in cdata.json")

            # Add SUBTYPES
            if metadata['SUBTYPES']:
                data['classes'][class_name]['SUBTYPES'] = metadata['SUBTYPES']
                stats['subtypes_added'] += 1
                print(f"  ✓ Added {len(metadata['SUBTYPES'])} subtypes")

            # Add CONTENT_FLAGS
            if metadata['CONTENT_FLAGS']:
                data['classes'][class_name]['CONTENT_FLAGS'] = metadata['CONTENT_FLAGS']
                stats['content_flags_added'] += 1
                print(f"  ✓ Added {len(metadata['CONTENT_FLAGS'])} content flags")
        else:
            stats['not_found'] += 1
            print(f"\n⚠  {class_name} NOT found in cdata.json (skipping)")

    # Write back to file
    print(f"\nWriting updated data to {cdata_path}...")
    with open(cdata_path, 'w') as f:
        json.dump(data, f, indent=2)

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Classes processed: {stats['found']}/{len(ALL_CLASS_METADATA)}")
    print(f"Classes with SUBTYPES added: {stats['subtypes_added']}")
    print(f"Classes with CONTENT_FLAGS added: {stats['content_flags_added']}")
    print(f"Classes not found: {stats['not_found']}")

    print("\n" + "="*70)
    print("COMPLETE METADATA")
    print("="*70)

    for class_name in sorted(ALL_CLASS_METADATA.keys()):
        if class_name in data['classes']:
            print(f"\n{class_name}:")

            if 'SUBTYPES' in data['classes'][class_name] and data['classes'][class_name]['SUBTYPES']:
                print("  Subtypes:")
                for name, info in sorted(data['classes'][class_name]['SUBTYPES'].items(),
                                        key=lambda x: x[1]['value']):
                    print(f"    {name} = {info['value']} # {info['description']}")

            if 'CONTENT_FLAGS' in data['classes'][class_name] and data['classes'][class_name]['CONTENT_FLAGS']:
                print("  Content Flags:")
                for name, info in sorted(data['classes'][class_name]['CONTENT_FLAGS'].items(),
                                        key=lambda x: x[1]['value']):
                    annotation = info.get('annotation', 'N/A')
                    print(f"    {name} = {info['value']} # {annotation}")
                    if info.get('columns'):
                        columns_str = ', '.join(info['columns'])
                        print(f"      Columns: [{columns_str}]")

    print("\n✓ Successfully updated cdata.json with all metadata")


if __name__ == '__main__':
    add_all_metadata()
