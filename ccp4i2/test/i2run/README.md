# i2run tests

Testing i2 using i2run command line scripts.

Examples:

```
ccp4-python -m pytest test/i2run
ccp4-python -m pytest test/i2run/test_8xfm.py
ccp4-python -m pytest test/i2run/test_8xfm.py -k servalcat
CCP4I2=/path/to/ccp4i2 ccp4-python -m pytest test/i2run/test_8xfm.py -k servalcat
```
