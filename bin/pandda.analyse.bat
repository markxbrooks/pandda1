@echo off
REM Adjust paths as needed

SET CBIN=C:\Users\MBrooks\CCP4-9\CCP4-9\CCP4\bin
SET PYTHONPATH=C:\Users\MBrooks\projects\pandda1\lib-python;%PYTHONPATH%

"%CBIN%\ccp4-python.bat" -m pandda.analyse %*
