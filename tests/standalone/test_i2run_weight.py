#!/usr/bin/env python
"""Test i2run with WEIGHT parameter that IS in children()."""
import os, sys
sys.path.insert(0, '/Users/nmemn/Developer/cdata-codegen')
sys.path.insert(0, '/Users/nmemn/Developer/cdata-codegen/server')
os.environ['CCP4I2_ROOT'] = '/Users/nmemn/Developer/cdata-codegen'
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'ccp4x.config.settings')

import django
django.setup()

from ccp4x.i2run.CCP4i2RunnerBase import CCP4i2RunnerBase

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
