@echo off
REM ----------------------------------------------------------------------------
REM make_tiles.bat YEAR1 YEAR2
REM   e.g. make_tiles.bat 2001 2020
REM ----------------------------------------------------------------------------

if "%~1"=="" (
  echo Usage: %~nx0 YEAR1 YEAR2
  exit /b 1
)
if "%~2"=="" (
  echo Usage: %~nx0 YEAR1 YEAR2
  exit /b 1
)

set "YEAR1=%~1"
set "YEAR2=%~2"

REM build all the intermediate filenames
set "INPUT=ocean_elevation_%YEAR2%_vs_%YEAR1%.tif"
set "TILEDIR=tiles_%YEAR2%_vs_%YEAR1%"

echo Processing %INPUT% into tile set %TILEDIR%\...

REM 6) remove any old tiles dir
if exist "%TILEDIR%" rmdir /s /q "%TILEDIR%"

REM 7) generate XYZ tiles zoom 0-6
gdal2tiles -z 0-6 -p mercator --processes=4 "%INPUT%" "%TILEDIR%"

echo Done! Tiles are in .\%TILEDIR%\{z}\{x}\{y}.png
