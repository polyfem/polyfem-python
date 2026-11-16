@echo off
set CMAKE_GENERATOR=Ninja
set CMAKE_GENERATOR_PLATFORM=
set CMAKE_GENERATOR_TOOLSET=
set CONDA_PREFIX=%PREFIX%
"%PYTHON%" -m pip install . --no-deps --no-build-isolation -v
if errorlevel 1 exit /b 1
