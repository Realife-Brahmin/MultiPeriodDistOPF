param(
    [switch]$IncludeLarge10kT12,
    [string]$TagPrefix = "optimized_timing_matrix"
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
$profileShort = "C:\Users\ARYANR~1"
$repoNative = if ($repo.StartsWith($env:USERPROFILE, [StringComparison]::OrdinalIgnoreCase)) {
    $profileShort + $repo.Substring($env:USERPROFILE.Length)
} else { $repo }

$julia = Join-Path $repoNative ".tools\julia-1.12.6\bin\julia.exe"
$depot = Join-Path $repoNative ".tools\julia_depot"
$extractor = Join-Path $repo "scripts\extract_filterddp_trace.ps1"
$resultDir = Join-Path $repo "ddp\results\network_filterddp"
$logDir = Join-Path $repo "logs"
$statusLog = Join-Path $logDir "$TagPrefix.status.log"
$summary = Join-Path $resultDir "$TagPrefix.csv"

New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$env:JULIA_DEPOT_PATH = $depot
$env:FILTERDDP_MAX_ITERATIONS = "200"
$env:FILTERDDP_OPTIMALITY_TOLERANCE = "1e-7"
$env:FILTERDDP_FACTOR_BACKED_POLICY = "1"
Remove-Item Env:FILTERDDP_BLOCKED_VALUE_RHS -ErrorAction SilentlyContinue
Remove-Item Env:FILTERDDP_VALUE_BLOCK_WIDTH -ErrorAction SilentlyContinue
Remove-Item Env:FILTERDDP_MEMORY_DIAGNOSTIC -ErrorAction SilentlyContinue
Remove-Item Env:FILTERDDP_TIMING_DIAGNOSTIC -ErrorAction SilentlyContinue

$cases = [Collections.Generic.List[object]]::new()
foreach ($system in @("ieee123C_1ph", "ieee2522C_1ph")) {
    foreach ($horizon in @(3, 6, 12, 24, 48, 96)) {
        $cases.Add([pscustomobject]@{ System = $system; Horizon = $horizon })
    }
}
foreach ($horizon in @(3, 6)) {
    $cases.Add([pscustomobject]@{ System = "large10kC_1ph"; Horizon = $horizon })
}
if ($IncludeLarge10kT12) {
    $cases.Add([pscustomobject]@{ System = "large10kC_1ph"; Horizon = 12 })
}

if (-not (Test-Path -LiteralPath $summary)) {
    "system,horizon,iterations,status,wall_s,solve_s,peak_working_set_MiB,trace" |
        Set-Content -LiteralPath $summary
}
"$(Get-Date -Format o) QUEUE_START cases=$($cases.Count) factor_backed=1 tolerance=1e-7" |
    Set-Content -LiteralPath $statusLog

foreach ($case in $cases) {
    $system = $case.System
    $horizon = $case.Horizon
    $tag = "${TagPrefix}_${system}_T${horizon}"
    $stdout = Join-Path $logDir "$tag.out.log"
    $stderr = Join-Path $logDir "$tag.err.log"
    $trace = Join-Path $resultDir "${tag}_trace.csv"
    $traceRelative = "ddp/results/network_filterddp/${tag}_trace.csv"

    $existing = Import-Csv -LiteralPath $summary | Where-Object {
        $_.system -eq $system -and [int]$_.horizon -eq $horizon -and [int]$_.status -eq 0
    }
    if ($existing -and (Test-Path -LiteralPath $trace)) {
        Add-Content -LiteralPath $statusLog -Value `
            "$(Get-Date -Format o) SKIP_COMPLETE system=$system T=$horizon"
        continue
    }

    Add-Content -LiteralPath $statusLog -Value `
        "$(Get-Date -Format o) RUN_START system=$system T=$horizon"
    $arguments = "--startup-file=no --project=envs/ddp2026 ddp/examples/power_system/ieee123c_filterddp.jl $system $horizon solve"
    $timer = [Diagnostics.Stopwatch]::StartNew()
    $process = Start-Process -FilePath $julia -ArgumentList $arguments `
        -WorkingDirectory $repoNative -WindowStyle Hidden -PassThru `
        -RedirectStandardOutput $stdout -RedirectStandardError $stderr
    $peakBytes = [int64]0
    while (-not $process.HasExited) {
        $process.Refresh()
        if ($process.WorkingSet64 -gt $peakBytes) { $peakBytes = $process.WorkingSet64 }
        Start-Sleep -Seconds 1
    }
    $process.WaitForExit()
    $process.Refresh()
    if ($process.WorkingSet64 -gt $peakBytes) { $peakBytes = $process.WorkingSet64 }
    $timer.Stop()

    $text = Get-Content -LiteralPath $stdout -Raw
    $match = [regex]::Match($text, 'solve complete:\s*([0-9.]+)\s*s,\s*iterations=([0-9]+),\s*status=([0-9]+)')
    $exitCode = $process.ExitCode
    # Windows PowerShell can return a null ExitCode after polling a redirected
    # native process.  A parsed FilterDDP status is the authoritative fallback.
    if ($null -eq $exitCode -and $match.Success) {
        $exitCode = if ([int]$match.Groups[3].Value -eq 0) { 0 } else { 1 }
    }
    if ($exitCode -ne 0 -or -not $match.Success) {
        Add-Content -LiteralPath $statusLog -Value `
            "$(Get-Date -Format o) RUN_FAIL system=$system T=$horizon exit=$exitCode parsed=$($match.Success)"
        throw "Optimized timing run failed: system=$system T=$horizon"
    }

    & $extractor -InputLog $stdout -OutputCsv $trace |
        Out-File (Join-Path $logDir "${tag}_trace.log")
    $solveSeconds = [double]$match.Groups[1].Value
    $iterations = [int]$match.Groups[2].Value
    $status = [int]$match.Groups[3].Value
    $peakMiB = $peakBytes / 1MB

    # Remove a stale row for the same case before appending a recovered run.
    $rows = Import-Csv -LiteralPath $summary | Where-Object {
        -not ($_.system -eq $system -and [int]$_.horizon -eq $horizon)
    }
    "system,horizon,iterations,status,wall_s,solve_s,peak_working_set_MiB,trace" |
        Set-Content -LiteralPath $summary
    foreach ($row in $rows) {
        Add-Content -LiteralPath $summary -Value ("{0},{1},{2},{3},{4},{5},{6},{7}" -f `
            $row.system,$row.horizon,$row.iterations,$row.status,$row.wall_s,$row.solve_s,
            $row.peak_working_set_MiB,$row.trace)
    }
    Add-Content -LiteralPath $summary -Value ("{0},{1},{2},{3},{4:F3},{5:F3},{6:F3},{7}" -f `
        $system,$horizon,$iterations,$status,$timer.Elapsed.TotalSeconds,$solveSeconds,$peakMiB,$traceRelative)
    Add-Content -LiteralPath $statusLog -Value `
        "$(Get-Date -Format o) RUN_END system=$system T=$horizon iterations=$iterations status=$status wall_s=$($timer.Elapsed.TotalSeconds) solve_s=$solveSeconds peak_MiB=$peakMiB"
}

Add-Content -LiteralPath $statusLog -Value "$(Get-Date -Format o) QUEUE_END"
