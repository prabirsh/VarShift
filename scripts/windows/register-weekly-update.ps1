<#
.SYNOPSIS
    Register a weekly VarSight data update in Windows Task Scheduler, for use with WSL.

.DESCRIPTION
    Cron inside WSL only runs while the WSL distribution is awake, so a Windows scheduled task is the
    dependable way to run the update on a laptop. The task starts WSL, runs one update, and exits. If
    the laptop was off at the scheduled time, Windows runs it at the next chance.

    Registering a task for your own user needs no administrator rights.

.EXAMPLE
    .\register-weekly-update.ps1 -RepoPath /mnt/i/VIP-Annotations/VarSight -DataDir /home/me/varsight-data

.EXAMPLE
    .\register-weekly-update.ps1 -Unregister
#>
[CmdletBinding(SupportsShouldProcess = $true)]
param(
    [string] $RepoPath,
    [string] $DataDir,
    [string] $VarsightCommand,
    [string] $Distribution,
    [ValidateSet('Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday')]
    [string] $Day = 'Sunday',
    [string] $Time = '09:00',
    [switch] $Unregister
)

$ErrorActionPreference = 'Stop'
$taskName = 'VarSight weekly data update'

if ($Unregister) {
    if ($PSCmdlet.ShouldProcess($taskName, 'Remove scheduled task')) {
        Unregister-ScheduledTask -TaskName $taskName -Confirm:$false
        Write-Output "Removed the task: $taskName"
    }
    return
}

foreach ($required in @{ RepoPath = $RepoPath; DataDir = $DataDir }.GetEnumerator()) {
    if (-not $required.Value) {
        throw "Please give -$($required.Key). Both paths are WSL paths, such as /home/me/varsight-data."
    }
}

function ConvertTo-ShellArgument([string] $Text) {
    return "'" + $Text.Replace("'", "'\''") + "'"
}

$script = ConvertTo-ShellArgument "$RepoPath/scripts/weekly-update.sh"
$dataArgument = ConvertTo-ShellArgument $DataDir
$command = "sh $script $dataArgument"
if ($VarsightCommand) {
    $command += ' ' + (ConvertTo-ShellArgument $VarsightCommand)
}

$wslArguments = @()
if ($Distribution) {
    $wslArguments += @('-d', $Distribution)
}
$wslArguments += @('-e', 'sh', '-lc', ('"' + $command.Replace('"', '\"') + '"'))

$action = New-ScheduledTaskAction -Execute 'wsl.exe' -Argument ($wslArguments -join ' ')
$trigger = New-ScheduledTaskTrigger -Weekly -DaysOfWeek $Day -At $Time
$settings = New-ScheduledTaskSettingsSet -StartWhenAvailable -RunOnlyIfNetworkAvailable -DontStopIfGoingOnBatteries -ExecutionTimeLimit (New-TimeSpan -Hours 6)

if ($PSCmdlet.ShouldProcess($taskName, 'Register scheduled task')) {
    Register-ScheduledTask -TaskName $taskName -Action $action -Trigger $trigger -Settings $settings `
        -Description 'Runs one VarSight data update. See docs/updates.md.' -Force | Out-Null
    Write-Output "Registered: $taskName, every $Day at $Time."
    Write-Output "Try it now with: Start-ScheduledTask -TaskName '$taskName'"
    Write-Output "Remove it with:  .\register-weekly-update.ps1 -Unregister"
}
