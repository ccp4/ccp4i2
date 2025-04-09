#!/usr/bin/env ccp4-python

import os

from . import common, parse


#if not "CRANK2" in os.environ:
#  sys.exit('Error: Environment variable CRANK2 not found.')

#if not os.path.isdir(os.getenv('CRANK2')):
#  sys.exit('Error: Environment variable CRANK2 points to non-existing directory {0}'.format(
#    os.getenv('CRANK2')))

version_file=os.path.join(os.path.dirname(parse.__file__), 'VERSION')
if not os.path.isfile(version_file):
  common.Error('VERSION file could not be found in {0}'.format(os.path.dirname(parse.__file__)))
with open(version_file) as f:
  common.Info('CRANK2, version {0}'.format(f.readline()))

if __name__ == "__main__":
  parse.parse().ParseAndRun()
