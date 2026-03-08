@echo off
set EXP=%1
if "%EXP%"=="" set EXP=hall
echo Recording field for experiment: %EXP%
.\build\WASAbiApp.exe --experiment %EXP% --mode sim-record-field
pause
