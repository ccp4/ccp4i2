#!/usr/bin/env python
"""
Test i2run parameter parsing and setting with modern approach.
"""
import sys

from ccp4i2.cli.i2run.CCP4i2RunnerBase import CCP4i2RunnerBase

def test_keyword_extraction():
    """Test that keywords are extracted from plugin."""
    print('Test 1: Keyword extraction for prosmart_refmac')
    print('=' * 60)

    keywords = CCP4i2RunnerBase.keywordsOfTaskName('prosmart_refmac', parent=None)
    print(f'✓ Found {len(keywords)} keywords')

    # Check for NCYCLES
    ncycles_kw = [kw for kw in keywords if kw['path'].endswith('.NCYCLES')]
    if ncycles_kw:
        kw = ncycles_kw[0]
        print(f'✓ NCYCLES keyword found:')
        print(f'  Path: {kw["path"]}')
        print(f'  Minimum path: {kw["minimumPath"]}')
        print(f'  Type: {kw["className"]}')
        return True
    else:
        print('✗ NCYCLES keyword not found')
        return False

def test_parameter_setting():
    """Test that parameters can be set via command line."""
    print('\nTest 2: Parameter setting via command line')
    print('=' * 60)

    # Simulate command line: prosmart_refmac --NCYCLES 20
    args = ['prosmart_refmac', '--NCYCLES', '20']
    print(f'Command line: {" ".join(args)}')

    # Create runner
    runner = CCP4i2RunnerBase(the_args=args, parent=None)

    # Parse arguments
    print('\nParsing arguments...')
    parsed_args = runner.parseArgs(arguments_parsed=False)
    print(f'✓ Parsed NCYCLES: {parsed_args.NCYCLES}')

    # Create plugin with arguments
    print('\nCreating plugin and setting parameters...')
    plugin = runner.pluginWithArgs(parsed_args, workDirectory='/tmp')

    # Check value
    ncycles_value = plugin.container.controlParameters.NCYCLES.value
    print(f'✓ NCYCLES in container: {ncycles_value}')

    if ncycles_value == 20:
        print('\n✅ SUCCESS! Parameter was set correctly')
        return True
    elif ncycles_value == 10:
        print('\n⚠️  WARNING: Parameter still at default (10), not set to 20')
        return False
    else:
        print(f'\n❓ UNEXPECTED: Value is {ncycles_value}')
        return False

if __name__ == '__main__':
    success1 = test_keyword_extraction()
    success2 = test_parameter_setting()

    print('\n' + '=' * 60)
    print('SUMMARY:')
    print(f'  Keyword extraction: {"✅ PASS" if success1 else "❌ FAIL"}')
    print(f'  Parameter setting:  {"✅ PASS" if success2 else "❌ FAIL"}')

    sys.exit(0 if (success1 and success2) else 1)
