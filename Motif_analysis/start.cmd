@echo off
setlocal
title Motif Analysis
echo Motif Analysis
cd /d "%~dp0"
docker info >nul 2>&1
if errorlevel 1 (
  echo Please start Docker Desktop and try again.
  pause
  exit /b 1
)
docker build -t motif-analysis .
if errorlevel 1 (
  echo Build failed. Check the message above.
  pause
  exit /b 1
)
echo.
echo Open http://localhost:8765 in your browser.
echo Stop with Ctrl+C. Keep this window open while using the tool.
echo.
docker run --rm -p 127.0.0.1:8765:8765 motif-analysis
if errorlevel 1 (
  echo Startup failed. Port 8765 may already be in use.
  pause
)
