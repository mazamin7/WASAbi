@echo off
set /p EXP_NAME=<experiment_name.txt
set EXP_DIR=experiments\%EXP_NAME%

echo Precomputing IR for experiment: %EXP_NAME%
echo Experiment directory: %EXP_DIR%

:: Define potential conda paths
set "CONDA_PATH="
if exist "C:\ProgramData\anaconda3\condabin\conda.bat" set "CONDA_PATH=C:\ProgramData\anaconda3\condabin\conda.bat"
if not defined CONDA_PATH if exist "%USERPROFILE%\anaconda3\condabin\conda.bat" set "CONDA_PATH=%USERPROFILE%\anaconda3\condabin\conda.bat"
if not defined CONDA_PATH if exist "C:\ProgramData\miniconda3\condabin\conda.bat" set "CONDA_PATH=C:\ProgramData\miniconda3\condabin\conda.bat"
if not defined CONDA_PATH if exist "%USERPROFILE%\miniconda3\condabin\conda.bat" set "CONDA_PATH=%USERPROFILE%\miniconda3\condabin\conda.bat"

if defined CONDA_PATH (
    echo Found Conda at: %CONDA_PATH%
    echo Activating 'wasabi' environment...
    
    if "%1"=="" (
        call "%CONDA_PATH%" run -n wasabi --no-capture-output python tools/auralizer/precompute_ir.py %EXP_DIR%
    ) else (
        call "%CONDA_PATH%" run -n wasabi --no-capture-output python tools/auralizer/precompute_ir.py %EXP_DIR% %1
    )
) else (
    echo Conda not found in standard locations, trying system PATH...
    where conda >nul 2>nul
    if %ERRORLEVEL% EQU 0 (
        if "%1"=="" (
            conda run -n wasabi --no-capture-output python tools/auralizer/precompute_ir.py %EXP_DIR%
        ) else (
            conda run -n wasabi --no-capture-output python tools/auralizer/precompute_ir.py %EXP_DIR% %1
        )
    ) else (
        echo Conda not found, falling back to system python...
        if "%1"=="" (
            python tools/auralizer/precompute_ir.py %EXP_DIR%
        ) else (
            python tools/auralizer/precompute_ir.py %EXP_DIR% %1
        )
    )
)

pause
