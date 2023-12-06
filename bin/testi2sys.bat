@echo off
call "%CCP4%\ccp4.setup.bat"

ccp4-python %CCP4I2%\bin\testi2sys.py %*
