@echo off
set EXP=%1
if "%EXP%"=="" (
    if exist experiment_name.txt (
        set /p EXP=<experiment_name.txt
    ) else (
        set EXP=hall
    )
)
set PFILE=%2
if "%PFILE%"=="" set PFILE=./experiments/%EXP%/output/record_data_0.bin

echo Visualizing playback for experiment: %EXP%
echo File: %PFILE%
.\build\WASAbiApp.exe --experiment %EXP% --mode viz-record --playback-file %PFILE%
pause
