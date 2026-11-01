<#
.SYNOPSIS
    Install the VarSight command on Windows, for PowerShell and Command Prompt.

.DESCRIPTION
    Needs only Python 3.10 or newer from python.org: no pip and no administrator rights.
    It writes a small varsight launcher to a folder of yours, and with your agreement adds that folder
    to your user PATH.

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File install.ps1

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File install.ps1 -NoSetup -NoPath
#>
[CmdletBinding()]
param(
    [string] $BinDir = (Join-Path $env:LOCALAPPDATA 'VarSight\bin'),
    [switch] $NoDemo,
    [switch] $NoSetup,
    [switch] $NoPath,
    [switch] $Test
)

$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $MyInvocation.MyCommand.Path

function Find-Python {
    $check = 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)'
    foreach ($candidate in @(@('py', '-3'), @('python'), @('python3'))) {
        $command = Get-Command $candidate[0] -ErrorAction SilentlyContinue
        if (-not $command) { continue }
        # The Microsoft Store placeholder for python exists but does nothing useful.
        if ($command.Source -like '*WindowsApps*') { continue }
        $arguments = @($candidate | Select-Object -Skip 1) + @('-c', $check)
        & $candidate[0] @arguments 2>$null
        if ($LASTEXITCODE -eq 0) {
            $found = & $candidate[0] @($candidate | Select-Object -Skip 1) -c 'import sys; print(sys.executable)'
            return $found.Trim()
        }
    }
    return $null
}

$python = Find-Python
if (-not $python) {
    Write-Error 'VarSight needs Python 3.10 or newer. Install it from https://www.python.org/downloads/ and tick "Add python.exe to PATH".'
}
Write-Output "Python: $(& $python -V) at $python"

New-Item -ItemType Directory -Force -Path $BinDir | Out-Null
$launcher = Join-Path $BinDir 'varsight.exe'
$cmdLauncher = Join-Path $BinDir 'varsight.cmd'
Remove-Item -Path $launcher, $cmdLauncher -Force -ErrorAction SilentlyContinue
& $python (Join-Path $repo 'scripts\windows\make_launcher.py') $BinDir "$repo\src" 2>$null | Out-Null
if ($LASTEXITCODE -ne 0 -or -not (Test-Path $launcher)) {
    # No pip to supply the launcher. A .cmd works too, but > in a variant needs care: see the README.
    $lines = @(
        '@echo off',
        'setlocal',
        "set `"PYTHONPATH=$repo\src`"",
        "`"$python`" -m varsight %*"
    )
    Set-Content -Path $cmdLauncher -Value $lines -Encoding ASCII
    $launcher = $cmdLauncher
}
Write-Output "Command: $launcher"
& $launcher version

if ($Test) {
    Write-Output ''
    Write-Output 'Running the unit tests'
    $env:PYTHONPATH = "$repo\src"
    $env:PYTHONDONTWRITEBYTECODE = '1'
    Push-Location $repo
    try {
        & $python -m unittest discover -s tests -t .
        if ($LASTEXITCODE -ne 0) { Write-Error 'Some tests failed.' }
    }
    finally {
        Pop-Location
    }
}

if (-not $NoDemo) {
    Write-Output ''
    Write-Output 'Running the demo (no download)'
    & $launcher demo "demodir=$repo\demo" "out=$((Get-Location).Path)\demo_output" replace
}

$userPath = [Environment]::GetEnvironmentVariable('Path', 'User')
$onPath = ($userPath -split ';') -contains $BinDir
if (-not $onPath -and -not $NoPath) {
    Write-Output ''
    $answer = Read-Host "Add $BinDir to your PATH so varsight works from any window? yes or no [yes]"
    if ($answer -eq '' -or $answer -match '^(y|yes)$') {
        [Environment]::SetEnvironmentVariable('Path', "$userPath;$BinDir", 'User')
        $env:Path = "$env:Path;$BinDir"
        Write-Output 'Added. Open a new window for it to take effect everywhere.'
    }
}

if (-not $NoSetup) {
    Write-Output ''
    $answer = Read-Host 'Set up the databases now? It asks whether you want hg19, hg38 or both. yes or no [no]'
    if ($answer -match '^(y|yes)$') { & $launcher setup }
}

Write-Output ''
Write-Output 'Installed. Next steps:'
Write-Output '  varsight setup                               choose hg19, hg38 or both, and download'
Write-Output '  varsight annotate "BRAF V600E" build=hg19    annotate a variant'
Write-Output '  varsight help                                see every command'
