@echo off
call "C:\CCP4-nightly\6.4\ccp4.setup.bat"
set CCP4I2_TOP = C:\Users\lizp\DEV\ccp4i2
set PYTHONSTARTUP=C:\Users\lizp\DEV\ccp4i2\bin\ccp4i2.pythonrc
ccp4-python
