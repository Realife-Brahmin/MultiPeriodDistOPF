param([switch]$SkipIeee2522T12)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
$profileShort = "C:\Users\ARYANR~1"
$repoNative = if ($repo.StartsWith($env:USERPROFILE, [StringComparison]::OrdinalIgnoreCase)) {
    $profileShort + $repo.Substring($env:USERPROFILE.Length)
} else { $repo }
$toolRoot = $repoNative
if (-not (Test-Path -LiteralPath (Join-Path $toolRoot ".tools\julia-1.12.6\bin\julia.exe"))) {
    # Linked worktrees share the main checkout's untracked tool installation.
    $toolRoot = [IO.Path]::GetFullPath((Join-Path $repoNative "..\..\.."))
}
$julia = Join-Path $toolRoot ".tools\julia-1.12.6\bin\julia.exe"
$env:JULIA_DEPOT_PATH = Join-Path $toolRoot ".tools\julia_depot"
$env:FILTERDDP_MAX_ITERATIONS = "200"
$env:FILTERDDP_OPTIMALITY_TOLERANCE = "1e-7"
$env:FILTERDDP_FACTOR_BACKED_POLICY = "1"
$env:FILTERDDP_TIMING_DIAGNOSTIC = "1"
$env:FILTERDDP_SKIP_SOLUTION_WRITE = "1"
Remove-Item Env:FILTERDDP_BLOCKED_VALUE_RHS -ErrorAction SilentlyContinue
Remove-Item Env:FILTERDDP_MEMORY_DIAGNOSTIC -ErrorAction SilentlyContinue

$logDir = Join-Path $repo "logs"
$resultDir = Join-Path $repo "ddp\results\network_filterddp"
$statusLog = Join-Path $logDir "iteration_timing_queue.status.log"
$summary = Join-Path $resultDir "iteration_timing_runs.csv"
$cases = @()
if (-not $SkipIeee2522T12) { $cases += ,@("ieee2522C_1ph", 12) }
$cases += ,@("large10kC_1ph", 3)

New-Item -ItemType Directory -Force -Path $logDir,$resultDir | Out-Null
"system,horizon,iterations,status,solve_s,timed_passes,timing_csv,trace_csv" | Set-Content -LiteralPath $summary
"$(Get-Date -Format o) QUEUE_START cases=$($cases.Count)" | Set-Content -LiteralPath $statusLog

foreach ($case in $cases) {
    $system = $case[0]; $horizon = [int]$case[1]
    $tag = "iteration_timing_${system}_T${horizon}"
    $stdout = Join-Path $logDir "$tag.out.log"
    $stderr = Join-Path $logDir "$tag.err.log"
    $timingCsv = Join-Path $resultDir "$tag.csv"
    $traceCsv = Join-Path $resultDir "${tag}_trace.csv"
    Add-Content -LiteralPath $statusLog -Value "$(Get-Date -Format o) RUN_START system=$system T=$horizon"
    $arguments = "--startup-file=no --project=envs/ddp2026 ddp/examples/power_system/ieee123c_filterddp.jl $system $horizon solve"
    $process = Start-Process -FilePath $julia -ArgumentList $arguments -WorkingDirectory $repoNative `
        -WindowStyle Hidden -PassThru -Wait -RedirectStandardOutput $stdout -RedirectStandardError $stderr
    $text = Get-Content -LiteralPath $stdout -Raw
    $match = [regex]::Match($text, 'solve complete:\s*([0-9.]+)\s*s,\s*iterations=([0-9]+),\s*status=([0-9]+)')
    if (-not $match.Success -or [int]$match.Groups[3].Value -ne 0) {
        Add-Content -LiteralPath $statusLog -Value "$(Get-Date -Format o) RUN_FAIL system=$system T=$horizon exit=$($process.ExitCode)"
        throw "Timing run failed: $system T=$horizon"
    }
    & (Join-Path $repo "scripts\extract_filterddp_trace.ps1") -InputLog $stdout -OutputCsv $traceCsv | Out-Null
    & (Join-Path $repo "scripts\extract_filterddp_timing.ps1") -InputLog $stdout -OutputCsv $timingCsv | Out-Null
    $timedPasses = (Import-Csv -LiteralPath $timingCsv).Count
    Add-Content -LiteralPath $summary -Value ("{0},{1},{2},{3},{4},{5},{6},{7}" -f `
        $system,$horizon,$match.Groups[2].Value,$match.Groups[3].Value,$match.Groups[1].Value,$timedPasses,
        "ddp/results/network_filterddp/$tag.csv","ddp/results/network_filterddp/${tag}_trace.csv")
    Add-Content -LiteralPath $statusLog -Value "$(Get-Date -Format o) RUN_END system=$system T=$horizon solve_s=$($match.Groups[1].Value) passes=$timedPasses"
}
Add-Content -LiteralPath $statusLog -Value "$(Get-Date -Format o) QUEUE_END"
