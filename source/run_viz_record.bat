@echo off
if "%~1"=="" (
    echo Usage: run_viz_record.bat [path_to_playback_file] [delay_ms]
    echo Example: run_viz_record.bat ./output/0.500000_1.000000_0.000000_0.000000/out_0.txt 50
    pause
    exit /b
)
set PLAYBACK_FILE=%~1
set DELAY=%~2
if "%DELAY%"=="" set DELAY=50

.\build\WASAbiApp.exe --mode viz-record --playback-file "%PLAYBACK_FILE%" --playback-delay %DELAY% --config ./config/config.ini %*
pause
