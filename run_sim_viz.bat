@echo off
set EXP=%1
if "%EXP%"=="" (
    if exist experiment_name.txt (
        set /p EXP=<experiment_name.txt
    ) else (
        set EXP=my_experiment
    )
)
echo Running experiment: %EXP%
.\build\WASAbiApp.exe --experiment %EXP% --mode sim-viz
pause
