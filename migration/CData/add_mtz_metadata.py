#!/usr/bin/env python3
"""
Add MTZ-specific metadata extracted from old CCP4XtalData.py to cdata.json.

This script adds CONTENT_FLAGS and SUBTYPES metadata to CMiniMtzDataFile subclasses.
"""

import json
from pathlib import Path

# Metadata extracted from /Users/nmemn/Developer/ccp4i2/core/CCP4XtalData.py

MTZ_CLASS_METADATA = {
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
    }
}


def add_mtz_metadata():
    """Add MTZ metadata to cdata.json."""

    cdata_path = Path(__file__).parent / 'cdata.json'

    print(f"Reading {cdata_path}...")
    with open(cdata_path, 'r') as f:
        data = json.load(f)

    # Add metadata to each class
    for class_name, metadata in MTZ_CLASS_METADATA.items():
        if class_name in data['classes']:
            print(f"\nAdding metadata to {class_name}:")

            # Add SUBTYPES
            if metadata['SUBTYPES']:
                data['classes'][class_name]['SUBTYPES'] = metadata['SUBTYPES']
                print(f"  ✓ Added {len(metadata['SUBTYPES'])} subtypes")

            # Add CONTENT_FLAGS
            if metadata['CONTENT_FLAGS']:
                data['classes'][class_name]['CONTENT_FLAGS'] = metadata['CONTENT_FLAGS']
                print(f"  ✓ Added {len(metadata['CONTENT_FLAGS'])} content flags")
        else:
            print(f"\n⚠ Warning: {class_name} not found in cdata.json")

    # Write back to file
    print(f"\nWriting updated data to {cdata_path}...")
    with open(cdata_path, 'w') as f:
        json.dump(data, f, indent=2)

    print("✓ Successfully updated cdata.json")

    # Print summary
    print("\n" + "="*60)
    print("SUMMARY OF ADDED METADATA")
    print("="*60)

    for class_name in MTZ_CLASS_METADATA:
        if class_name in data['classes']:
            print(f"\n{class_name}:")

            if 'SUBTYPES' in data['classes'][class_name]:
                print("  Subtypes:")
                for name, info in data['classes'][class_name]['SUBTYPES'].items():
                    print(f"    {name} = {info['value']} ({info['description']})")

            if 'CONTENT_FLAGS' in data['classes'][class_name]:
                print("  Content Flags:")
                for name, info in data['classes'][class_name]['CONTENT_FLAGS'].items():
                    columns_str = ', '.join(info['columns'])
                    print(f"    {name} = {info['value']}: {info['annotation']}")
                    print(f"      Columns: [{columns_str}]")


if __name__ == '__main__':
    add_mtz_metadata()
