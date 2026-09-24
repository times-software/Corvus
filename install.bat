setlocal EnableDelayedExpansion

REM ===========================================================================
REM Environment Selection
REM ===========================================================================

where conda >nul 2>nul

if %ERRORLEVEL%==0 (
    echo.
    echo Conda detected.
    echo.
    set /p USE_CONDA="Use Conda (recommended)? [Y/n] "

    if /I "!USE_CONDA!"=="n" (
        set USE_CONDA=0
    ) else (
        set USE_CONDA=1
    )
) else (
    echo.
    echo Conda was not found.
    echo.

    set /p USE_VENV="Use a Python venv instead? [Y/n] "

    if /I "!USE_VENV!"=="n" (
        echo.
        echo Please install Miniforge:
        echo https://conda-forge.org/miniforge/
        exit /b 1
    )

    set USE_CONDA=0
)

REM ===========================================================================
REM Get environment name from current directory
REM ===========================================================================

REM for %%I in ("%CD%") do set ENV_NAME=%%~nxI
set "ENV_NAME=Corvus"

REM ===========================================================================
REM Conda path
REM ===========================================================================


if "!USE_CONDA!"=="1" (
    echo
    echo Creating Conda environment "!ENV_NAME!"...

    call conda create -n "!ENV_NAME!" "python>=3.12,<3.14" pip
    if errorlevel 1 (
        echo ERROR: Conda environment creation failed.
        
        exit /b 1
    )
    @echo on
    echo Activating environment "!ENV_NAME!"
    call conda activate "!ENV_NAME!"

    call conda env list

    if errorlevel 1 (
        echo ERROR: Conda activation failed.
        
        exit /b 1
    )
    echo Before
    set PYTHON=python
    echo "!PYTHON!"
    
    "!PYTHON!" -c "import sys; exit(0 if (3,12) <= sys.version_info[:2] < (3,14) else 1)"


    if errorlevel 1 (
        echo "ERROR: Python version is not in the range [3.12, 3.14)."
        
        exit /b 1
    )
    echo "After"
    

) else (

REM ===========================================================================
REM venv path
REM ===========================================================================

    where python >nul 2>nul

    if errorlevel 1 (
        echo ERROR: Python not found.
        
        exit /b 1
    )

    if exist .venv (
        echo ERROR: .venv already exists.
        
        exit /b 1
    )

    echo.
    echo Creating virtual environment...

    python -m venv .venv

    if errorlevel 1 (
        echo ERROR: venv creation failed.
        
        exit /b 1
    )

    if not exist .venv\Scripts\python.exe (
        echo ERROR: venv not created correctly.
        
        exit /b 1
    )

    set PYTHON=.venv\Scripts\python.exe
)

REM ===========================================================================
REM Upgrade packaging tools
REM ===========================================================================

echo.
echo Upgrading packaging tools...

"!PYTHON!" -m pip install --upgrade pip setuptools wheel

if errorlevel 1 (
    echo ERROR: pip upgrade failed.
    
    exit /b 1
)

REM ===========================================================================
REM Install current package
REM ===========================================================================

echo.
echo Installing current package...

"!PYTHON!" -m pip install .

if errorlevel 1 (
    echo ERROR: package installation failed.
    
    exit /b 1
)
echo.
echo.
REM INSTALL examples folder
@echo off
setlocal enabledelayedexpansion

echo ===================================================
echo   Locating Desktop and Deploying Examples Folder
echo ===================================================

REM 1. Define the source folder name (relative to where this script runs)
set "SOURCE_DIR=%~dp0examples"
set "TARGET_FOLDER_NAME=corvus_examples"

REM 2. Query the registry for the true Desktop path (handles Network/OneDrive/Local)
for /f "tokens=2*" %%A in ('reg query "HKCU\Software\Microsoft\Windows\CurrentVersion\Explorer\User Shell Folders" /v Desktop 2^>nul') do (
    set "RAW_DESKTOP_PATH=%%B"
)

REM 3. Expand any nested environment variables in the registry path
for /f "delims=" %%I in ('echo !RAW_DESKTOP_PATH!') do set "TRUE_DESKTOP=%%I"

REM 4. Construct the absolute final destination path
set "DESTINATION_DIR=!TRUE_DESKTOP!\!TARGET_FOLDER_NAME!"

echo Source Directory:      "%SOURCE_DIR%"
echo Target Desktop Path:   "%TRUE_DESKTOP%"
echo Destination Absolute:  "%DESTINATION_DIR%"
echo ---------------------------------------------------

REM 5. Validation Check: Does the source "examples" folder actually exist?
if not exist "%SOURCE_DIR%" (
    echo [ERROR] Source folder "examples" not found at: "%SOURCE_DIR%"
    echo Please ensure the "examples" directory is next to this installer script.
    goto :EXIT_FAILURE
)

REM 6. Perform the safe copy operation using Robocopy
echo Copying files, please wait...
REM /E      = Copies subdirectories, including empty ones.
REM /R:3    = Retry 3 times on failed file locks (crucial for network glitches).
REM /W:5    = Wait 5 seconds between retries.
REM /MT:16  = Multithreaded copy (faster over networks).
robocopy "%SOURCE_DIR%" "%DESTINATION_DIR%" /E /R:3 /W:5 /MT:16 >nul

REM 7. Evaluate the Robocopy Exit Code
REM Robocopy exit codes below 8 indicate success (0=no changes, 1=files copied, 2/3=tweaks)
if %ERRORLEVEL% LSS 8 (
    echo [SUCCESS] Successfully deployed "%TARGET_FOLDER_NAME%" to the Desktop.
    goto :EXIT_SUCCESS
) else (
    echo [ERROR] Robocopy failed with Exit Code %ERRORLEVEL%. 
    echo This could be due to a disconnected network home drive or missing permissions.
)

:EXIT_FAILURE
:EXIT_SUCCESS

REM ===========================================================================
REM Optional SciGUI install
REM ===========================================================================

echo.
set /p INSTALL_SCIGUI="Install SciGUI? [y/N] "

if /I "!INSTALL_SCIGUI!"=="y" goto INSTALL_SCIGUI
if /I "!INSTALL_SCIGUI!"=="yes" goto INSTALL_SCIGUI
goto FINISH

:INSTALL_SCIGUI

set TMPDIR=%TEMP%\SciGUI_%RANDOM%

mkdir "%TMPDIR%"

echo.
echo Downloading SciGUI...

powershell -NoProfile -ExecutionPolicy Bypass ^
    -Command "Invoke-WebRequest -Uri 'https://github.com/times-software/SciGUI/archive/refs/heads/main.zip' -OutFile '%TMPDIR%\scigui.zip'"

if errorlevel 1 (
    echo ERROR: download failed.
    rmdir /s /q "%TMPDIR%"
    
    exit /b 1
)

echo Extracting SciGUI...

powershell -NoProfile -ExecutionPolicy Bypass ^
    -Command "Expand-Archive '%TMPDIR%\scigui.zip' '%TMPDIR%\src'"

if errorlevel 1 (
    echo ERROR: extraction failed.
    rmdir /s /q "%TMPDIR%"
    
    exit /b 1
)

for /f "delims=" %%F in (
    'dir /s /b "%TMPDIR%\src\setup.py"'
) do (
    set SETUPPY=%%F
    goto FOUND_SETUP
)

echo ERROR: setup.py not found.
rmdir /s /q "%TMPDIR%"

exit /b 1

:FOUND_SETUP

for %%F in ("!SETUPPY!") do set SCIGUI_DIR=%%~dpF

echo Installing SciGUI...

pushd "!SCIGUI_DIR!"

"!PYTHON!" -m pip install .

if errorlevel 1 (
    popd
    rmdir /s /q "%TMPDIR%"
    echo ERROR: SciGUI installation failed.
    
    exit /b 1
)

popd

rmdir /s /q "%TMPDIR%"

echo SciGUI installed successfully.

:FINISH

echo.
echo Setup complete.
echo.

if "!USE_CONDA!"=="1" (
    REM ===========================================================================
    REM Create desktop launcher
    REM ===========================================================================

    set "DESKTOP=!TRUE_DESKTOP!

    set "CONDA_BAT=%CONDA_PREFIX%\condabin\conda.bat"
    set "LAUNCHER=%DESKTOP%\%ENV_NAME%.bat"

    if not exist "%DESKTOP%" (
        echo ERROR: Desktop directory not found:
        echo   %DESKTOP%
        exit /b 1
    )
    
    REM Check if a network home share exists (e.g., \\server\share)
    if not "%HOMESHARE%"=="" (
        set "TRUE_HOME=%HOMESHARE%"
    ) else if not "%HOMEDRIVE%"=="" (
        REM If it's a mapped drive letter (e.g., H:\path)
        set "TRUE_HOME=%HOMEDRIVE%%HOMEPATH%"
    ) else (
        REM Fallback to local profile if no network home is defined
        set "TRUE_HOME=%USERPROFILE%"
    )
    set "PROJECT_DIR=!TRUE_HOME!\"

    echo @echo off
    echo cd /d "%PROJECT_DIR%"
    echo call "%CONDA_BAT%" activate "%ENV_NAME%"
    echo.
    echo if errorlevel 1 ^(
    echo     echo ERROR: Failed to activate environment %ENV_NAME%
    echo     pause
    echo     exit /b 1
    echo ^)
    echo.
    echo title %ENV_NAME%
    echo echo Activated conda environment: %ENV_NAME%
    echo echo.
    echo cmd /k
    ) > "%LAUNCHER%"
    
    if errorlevel 1 (
        echo ERROR: Failed to create launcher:
        echo   %LAUNCHER%
        exit /b 1
    )
    
    if not exist "%LAUNCHER%" (
        echo ERROR: Launcher was not created:
        echo   %LAUNCHER%
        exit /b 1
    )
    
    echo.
    echo Created desktop launcher:
    echo   %LAUNCHER%
    echo.
) else (
    echo To activate later:
    echo     .venv\Scripts\activate
)


endlocal

