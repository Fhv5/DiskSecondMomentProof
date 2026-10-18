function Invoke-LoggedNativeCommand {
    [CmdletBinding()]
    param(
        [Parameter(Mandatory = $true)][string]$Executable,
        [Parameter(Mandatory = $true)][string[]]$Arguments,
        [Parameter(Mandatory = $true)][string]$LogPath,
        [Parameter(Mandatory = $true)][string]$Description
    )

    $application = Get-Command -Name $Executable -CommandType Application -ErrorAction Stop |
        Select-Object -First 1
    # Windows PowerShell turns redirected native stderr into ErrorRecords.
    # A message on stderr is not itself a failed process; inspect its exit code.
    # These preferences are local to this function, leaving the caller on Stop.
    $ErrorActionPreference = 'Continue'
    $PSNativeCommandUseErrorActionPreference = $false
    $global:LASTEXITCODE = $null
    & $application.Source @Arguments 2>&1 |
        ForEach-Object { $_.ToString() } |
        Tee-Object -FilePath $LogPath -ErrorAction Stop
    $nativeExitCode = $global:LASTEXITCODE
    if ($null -eq $nativeExitCode) {
        throw "$Description did not return an exit code. See $LogPath."
    }
    if ($nativeExitCode -ne 0) {
        throw "$Description failed with exit code $nativeExitCode. See $LogPath."
    }
}
