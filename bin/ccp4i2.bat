@echo off
call "%CCP4%\ccp4.setup.bat"
set CCP4I2="%~dp0\.."
ccp4-python "%~dp0\browser.py"

