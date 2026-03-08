@echo off
set EXP=%1
if "%EXP%"=="" set EXP=my_experiment
echo Running experiment: %EXP%
.\build\WASAbiApp.exe --experiment %EXP% --mode sim-viz
pause
