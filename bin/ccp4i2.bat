@echo off
setlocal

REM ===================================================

SET caller=%~dpf0
CALL "%~dp0..\ccp4.setup.bat"

SET CCP4I2=%CCP4%\Python39\Lib\site-packages\ccp4i2
SET CCP4I2_TOP=%CCP4I2%

SET PYTHONSTARTUP=%CCP4I2%\bin\ccp4i2.pythonrc

REM LIBTBX_BUILD is used in cctbx to find libtbx_env file (pickle)
REM SET LIBTBX_BUILD=%CCP4%\lib\cctbx

ccp4-python.bat "%CCP4I2%\bin\browser.py" %*
