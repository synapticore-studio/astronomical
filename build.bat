@echo off
REM Simple build script for Windows
REM Following project standards - direct approach only

echo Build Astronomical Extension
echo.

REM Check Python
python --version >nul 2>&1
if errorlevel 1 (
    echo Python not found
    exit /b 1
)

echo Python version:
python --version
echo.

echo Choose build:
echo 1. Quick build
echo 2. Full build
echo 3. Clean build
echo.

set /p choice="Choice (1-3): "

if "%choice%"=="1" (
    python quick_build.py
) else if "%choice%"=="2" (
    python build_extension.py
) else if "%choice%"=="3" (
    python build_extension.py --clean
) else (
    echo Invalid choice
    exit /b 1
)

echo.
echo Build complete
pause
