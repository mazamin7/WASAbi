@echo off
set EXP=%1
if "%EXP%"=="" (
    if exist experiment_name.txt (
        set /p EXP=<experiment_name.txt
    ) else (
        set EXP=hall
    )
)
echo Recording response for experiment: %EXP%
.\build\WASAbiApp.exe --experiment %EXP% --mode sim-record-response
pause
