param(
    [string]$OutputRoot = '',
    [ValidateSet('2', '3', '4', '5', '6')]
    [ValidateCount(1, 5)]
    [int[]]$DiskCounts = @(2, 3, 4, 5, 6)
)

$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'native_logging.ps1')
if (@($DiskCounts | Select-Object -Unique).Count -ne $DiskCounts.Count) {
    throw 'DiskCounts must not contain duplicates.'
}
$proofRoot = Split-Path -Parent $PSScriptRoot
$proofScript = Join-Path $proofRoot 'scripts/run_proof.jl'
$censusScript = Join-Path $proofRoot 'scripts/filter_census.jl'
$threadSetup = Join-Path $proofRoot 'scripts/single_threaded.jl'
$juliaRuntimeArguments = @(
    '--startup-file=no', '--threads=1,0', '--gcthreads=1',
    "--project=$proofRoot", "--load=$threadSetup"
)
$revision = (& git -C $proofRoot rev-parse HEAD).Trim()
if ($LASTEXITCODE -ne 0) { throw 'Unable to identify the Git revision.' }
$dirty = & git -C $proofRoot status --porcelain
if ($LASTEXITCODE -ne 0 -or $dirty) {
    throw 'Commit all source changes before creating a final evidence run.'
}
if (-not $OutputRoot) { $OutputRoot = Join-Path $proofRoot 'output/releases' }
$runId = ([DateTime]::UtcNow.ToString('yyyyMMddTHHmmssZ')) + '-' + $revision.Substring(0, 12)
$runDirectory = Join-Path $OutputRoot $runId
if (Test-Path -LiteralPath $runDirectory) { throw "Run directory already exists: $runDirectory" }
New-Item -ItemType Directory -Path $runDirectory | Out-Null
$sourceArchive = Join-Path $runDirectory 'source.zip'
& git -C $proofRoot archive --format=zip "--output=$sourceArchive" HEAD
if ($LASTEXITCODE -ne 0) { throw 'Unable to archive the source revision.' }

foreach ($diskCount in $DiskCounts) {
    $planarInput = Join-Path $proofRoot "data/planar_graphs/planar$diskCount.g6"
    $censusLog = Join-Path $runDirectory "n$diskCount-census.log"
    Invoke-LoggedNativeCommand -Executable 'julia' -Arguments ($juliaRuntimeArguments + @(
        $censusScript,
        '--n', "$diskCount", '--input', $planarInput
    )) -LogPath $censusLog -Description "Census verification for n=$diskCount"
    $proofLog = Join-Path $runDirectory "n$diskCount-proof.log"
    Invoke-LoggedNativeCommand -Executable 'julia' -Arguments ($juliaRuntimeArguments + @(
        $proofScript,
        '--n', "$diskCount", '--precision', '512', '--record-evidence', '--output', $runDirectory
    )) -LogPath $proofLog -Description "Proof run for n=$diskCount"
}

$digests = @{}
Get-ChildItem -LiteralPath $runDirectory -File | ForEach-Object {
    $digests[$_.Name] = (Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
}
foreach ($diskCount in $DiskCounts) {
    $path = Join-Path $runDirectory "n$diskCount/run.json"
    $digests["n$diskCount/run.json"] = (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant()
}
$releaseRecord = [ordered]@{
    schema_version = 1
    disk_counts = @($DiskCounts)
    run_id = $runId
    git_revision = $revision
    created_utc = [DateTime]::UtcNow.ToString('o')
    source_archive = 'source.zip'
    thread_policy = [ordered]@{ julia = 1; interactive = 0; gc = 1; blas = 1; flint = 1 }
    file_sha256 = $digests
}
$releaseJson = $releaseRecord | ConvertTo-Json -Depth 5
[System.IO.File]::WriteAllText((Join-Path $runDirectory 'release.json'),
    $releaseJson, [System.Text.UTF8Encoding]::new($false))
Write-Output "Final evidence run: $runDirectory"
