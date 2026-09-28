"""Allow ``ccp4-python -m ccp4i2.cli.i2run`` as well as the ``i2run`` script.

The bare module path is the form the job panel's "i2run command" button
renders: the ``i2run`` console script is shadowed on a normal CCP4 setup by
the legacy Qt ``$CCP4/bin/i2run``, whereas ``-m`` names this package
unambiguously.
"""

from .main import main

main()
