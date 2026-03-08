@echo off
set EXP=%1
if "%EXP%"=="" set EXP=hall
echo Recording response for experiment: %EXP%
.\build\WASAbiApp.exe --experiment %EXP% --mode sim-record-response
pause
