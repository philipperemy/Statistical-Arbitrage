@echo off
if "%1"=="" echo Missing file specification>&2 & exit /b 1
PowerShell -C "[System.BitConverter]::ToString([System.Security.Cryptography.MD5]::Create().ComputeHash([System.IO.File]::ReadAllBytes('%1'))).ToLowerInvariant().Replace('-', '')