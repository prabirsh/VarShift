<#
.SYNOPSIS
    Register a weekly VarShift data update in Windows Task Scheduler.

.DESCRIPTION
    If the computer was off at the scheduled time, Windows runs the task at the next chance. Registering
    a task for your own user needs no administrator rights.

    Native Windows install: give -DataDir only.
    WSL install: also give -WslCommand, because cron inside WSL only runs while WSL is awake.

.EXAMPLE
    .\register-weekly-update.ps1 -DataDir C:\Users\me\varshift-data

.EXAMPLE
    .\register-weekly-update.ps1 -DataDir /home/me/varshift-data -WslCommand /home/me/.local/bin/varshift

.EXAMPLE
    .\register-weekly-update.ps1 -Unregister
#>
[CmdletBinding(SupportsShouldProcess = $true)]
param(
    [string] $DataDir,
    [string] $WslCommand,
    [string] $Distribution,
    [ValidateSet('Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday')]
    [string] $Day = 'Sunday',
    [string] $Time = '09:00',
    [switch] $Unregister
)

$ErrorActionPreference = 'Stop'
$taskName = 'VarShift weekly data update'

if ($Unregister) {
    if ($PSCmdlet.ShouldProcess($taskName, 'Remove scheduled task')) {
        Unregister-ScheduledTask -TaskName $taskName -Confirm:$false
        Write-Output "Removed the task: $taskName"
    }
    return
}

if (-not $DataDir) {
    throw 'Please give -DataDir, the folder you chose at setup.'
}

if ($WslCommand) {
    $wslArguments = @()
    if ($Distribution) { $wslArguments += @('-d', $Distribution) }
    $wslArguments += @('-e', 'sh', '-lc', ('"' + "'$WslCommand' update data='$DataDir'" + '"'))
    $action = New-ScheduledTaskAction -Execute 'wsl.exe' -Argument ($wslArguments -join ' ')
}
else {
    $folder = Join-Path $env:LOCALAPPDATA 'VarShift\bin'
    $launcher = @('varshift.exe', 'varshift.cmd') | ForEach-Object { Join-Path $folder $_ } | Where-Object { Test-Path $_ } | Select-Object -First 1
    if (-not $launcher) {
        throw "Cannot find the varshift command in $folder. Run install.ps1 first, or give -WslCommand for a WSL install."
    }
    $action = New-ScheduledTaskAction -Execute $launcher -Argument "update `"data=$DataDir`""
}

$trigger = New-ScheduledTaskTrigger -Weekly -DaysOfWeek $Day -At $Time
$settings = New-ScheduledTaskSettingsSet -StartWhenAvailable -RunOnlyIfNetworkAvailable -DontStopIfGoingOnBatteries -ExecutionTimeLimit (New-TimeSpan -Hours 6)

if ($PSCmdlet.ShouldProcess($taskName, 'Register scheduled task')) {
    Register-ScheduledTask -TaskName $taskName -Action $action -Trigger $trigger -Settings $settings `
        -Description 'Runs one VarShift data update. See docs/updates.md.' -Force | Out-Null
    Write-Output "Registered: $taskName, every $Day at $Time."
    Write-Output "Try it now with: Start-ScheduledTask -TaskName '$taskName'"
    Write-Output "Remove it with:  .\register-weekly-update.ps1 -Unregister"
}
