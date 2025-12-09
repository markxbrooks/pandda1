@echo off
REM Adjust paths as needed

SET VENV=C:\Users\MBrooks\projects\pandda-bitbucket\venv
SET PYTHONPATH=C:\Users\MBrooks\projects\pandda-bitbucket\lib-python;%PYTHONPATH%

"ccp4-python" -m pandda.analyse %*
