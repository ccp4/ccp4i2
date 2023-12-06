@echo off
REM call "%CCP4%\ccp4.setup.bat"
set CCP4I2="%~dp0\.."
echo %CCP4I2%
ccp4-python ../bin/browser.py

