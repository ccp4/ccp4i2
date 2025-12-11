#!/usr/bin/env python
"""Test i2run with WEIGHT parameter that IS in children()."""
import os, sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / 'server'))
os.environ['CCP4I2_ROOT'] = str(PROJECT_ROOT)
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'ccp4i2.config.settings')

import django
django.setup()

from ccp4i2.i2run.CCP4i2RunnerBase import CCP4i2RunnerBase

# Test with WEIGHT (which we know is in children())
args = ['prosmart_refmac', '--WEIGHT', 'AUTO']
runner = CCP4i2RunnerBase(the_args=args, parent=None)

print('Testing WEIGHT parameter...')
parsed_args = runner.parseArgs(arguments_parsed=False)
print(f'✓ Parsed WEIGHT: {parsed_args.WEIGHT}')

plugin = runner.pluginWithArgs(parsed_args, workDirectory='/tmp')
weight_value = plugin.container.controlParameters.WEIGHT.value
print(f'✓ WEIGHT in container: {weight_value}')

if weight_value == 'AUTO':
    print('\n✅ SUCCESS! Parameter setting works with modern approach!')
else:
    print(f'\n⚠️  Value mismatch: expected AUTO, got {weight_value}')
