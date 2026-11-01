@echo off
rem Install VarSight from Command Prompt. This only starts install.ps1, which does the work.
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0install.ps1" %*
