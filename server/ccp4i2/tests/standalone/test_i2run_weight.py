#!/usr/bin/env python
"""Test i2run with WEIGHT parameter that IS in children()."""
from ccp4i2.cli.i2run.CCP4i2RunnerBase import CCP4i2RunnerBase

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
