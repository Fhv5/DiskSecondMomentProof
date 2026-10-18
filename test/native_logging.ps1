param(
    [string]$HelperPath = (Join-Path (Split-Path -Parent $PSScriptRoot) 'scripts/native_logging.ps1')
)

$ErrorActionPreference = 'Stop'
. $HelperPath
$childShell = Join-Path $env:SystemRoot 'System32/WindowsPowerShell/v1.0/powershell.exe'
$testDirectory = Join-Path ([System.IO.Path]::GetTempPath()) ('dsmp-native-logging-' + [Guid]::NewGuid())
New-Item -ItemType Directory -Path $testDirectory | Out-Null

function Assert-LoggingCondition([bool]$Condition, [string]$Message) {
    if (-not $Condition) { throw $Message }
}

try {
    $successLog = Join-Path $testDirectory 'success.log'
    $successOutput = @(Invoke-LoggedNativeCommand -Executable $childShell -Arguments @(
        '-NoProfile', '-Command',
        "[Console]::Error.WriteLine('Precompiling packages...'); [Console]::Out.WriteLine('completed'); exit 0"
    ) -LogPath $successLog -Description 'Successful test process')
    $successText = Get-Content -LiteralPath $successLog -Raw
    Assert-LoggingCondition ($successText.Contains('Precompiling packages...')) 'Native stderr is missing from the log.'
    Assert-LoggingCondition ($successText.Contains('completed')) 'Native stdout is missing from the log.'
    Assert-LoggingCondition ($successOutput -contains 'Precompiling packages...') 'Native stderr is missing from console output.'
    Assert-LoggingCondition ($successOutput -contains 'completed') 'Native stdout is missing from console output.'
    Assert-LoggingCondition ($ErrorActionPreference -eq 'Stop') 'The helper changed the caller error preference.'

    $failureLog = Join-Path $testDirectory 'failure.log'
    $failure = $null
    try {
        Invoke-LoggedNativeCommand -Executable $childShell -Arguments @(
            '-NoProfile', '-Command', "[Console]::Error.WriteLine('genuine failure'); exit 17"
        ) -LogPath $failureLog -Description 'Failing test process' | Out-Null
    } catch { $failure = $_.Exception.Message }
    Assert-LoggingCondition ($null -ne $failure -and $failure -like '*failed with exit code 17*') 'A nonzero native exit code did not stop execution.'
    Assert-LoggingCondition ((Get-Content -LiteralPath $failureLog -Raw).Contains('genuine failure')) 'Failure stderr is missing from the log.'

    $missingExecutableFailed = $false
    try {
        Invoke-LoggedNativeCommand -Executable (Join-Path $testDirectory 'missing.exe') -Arguments @('unused') -LogPath $successLog -Description 'Missing process' | Out-Null
    } catch { $missingExecutableFailed = $true }
    Assert-LoggingCondition $missingExecutableFailed 'An unavailable executable did not stop execution.'

    $invalidLogFailed = $false
    try {
        Invoke-LoggedNativeCommand -Executable $childShell -Arguments @(
            '-NoProfile', '-Command', "[Console]::Out.WriteLine('cannot save this'); exit 0"
        ) -LogPath $testDirectory -Description 'Invalid log destination' | Out-Null
    } catch { $invalidLogFailed = $true }
    Assert-LoggingCondition $invalidLogFailed 'A log write failure did not stop execution.'
    Write-Output "Native logging checks passed in PowerShell $($PSVersionTable.PSVersion)."
} finally {
    Get-ChildItem -LiteralPath $testDirectory -File | ForEach-Object {
        Remove-Item -LiteralPath $_.FullName -Force
    }
    Remove-Item -LiteralPath $testDirectory
}
