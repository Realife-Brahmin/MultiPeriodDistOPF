param(
    [int]$Horizon = 3,
    [int]$Repetitions = 2,
    [string]$TagPrefix = "idle_policy"
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
# Windows PowerShell 5 can split native-command paths at spaces even when the
# same invocation works in pwsh. Use the equivalent 8.3 path for unattended
# background runs launched through powershell.exe.
$profileShort = "C:\Users\ARYANR~1"
$repoNative = if ($repo.StartsWith($env:USERPROFILE, [StringComparison]::OrdinalIgnoreCase)) {
    $profileShort + $repo.Substring($env:USERPROFILE.Length)
} else {
    $repo
}
$julia = Join-Path $repoNative ".tools\julia-1.12.6\bin\julia.exe"
$depot = Join-Path $repoNative ".tools\julia_depot"
$project = Join-Path $repo "envs\ddp2026"
$driver = Join-Path $repo "ddp\examples\power_system\ieee123c_filterddp.jl"
$extractor = Join-Path $repo "scripts\extract_filterddp_trace.ps1"
$resultDir = Join-Path $repo "ddp\results\network_filterddp"
$logDir = Join-Path $repo "logs"
$statusLog = Join-Path $logDir "$TagPrefix.status.log"
$summary = Join-Path $resultDir "${TagPrefix}_T${Horizon}_runs.csv"

New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$env:JULIA_DEPOT_PATH = $depot
$env:FILTERDDP_MAX_ITERATIONS = "200"
$env:FILTERDDP_OPTIMALITY_TOLERANCE = "1e-7"
Remove-Item Env:FILTERDDP_MEMORY_DIAGNOSTIC -ErrorAction SilentlyContinue
Remove-Item Env:FILTERDDP_TIMING_DIAGNOSTIC -ErrorAction SilentlyContinue
Remove-Item Env:FILTERDDP_CAPTURE_KKT -ErrorAction SilentlyContinue
Remove-Item Env:FILTERDDP_CAPTURE_STAGE -ErrorAction SilentlyContinue

"mode,repetition,horizon,wall_s,peak_working_set_MiB,exit_code,trace" | Set-Content -LiteralPath $summary
"$(Get-Date -Format o) QUEUE_START T=$Horizon repetitions=$Repetitions" | Set-Content -LiteralPath $statusLog

# Reverse the order on alternating repetitions to reduce systematic thermal or
# background-load bias without overlapping any solves.
for ($rep = 1; $rep -le $Repetitions; $rep++) {
    $modes = if ($rep % 2 -eq 1) { @("dense", "factor", "blocked") } else { @("blocked", "factor", "dense") }
    foreach ($mode in $modes) {
        if ($mode -eq "dense") {
            Remove-Item Env:FILTERDDP_FACTOR_BACKED_POLICY -ErrorAction SilentlyContinue
            Remove-Item Env:FILTERDDP_BLOCKED_VALUE_RHS -ErrorAction SilentlyContinue
            Remove-Item Env:FILTERDDP_VALUE_BLOCK_WIDTH -ErrorAction SilentlyContinue
        } elseif ($mode -eq "factor") {
            $env:FILTERDDP_FACTOR_BACKED_POLICY = "1"
            Remove-Item Env:FILTERDDP_BLOCKED_VALUE_RHS -ErrorAction SilentlyContinue
            Remove-Item Env:FILTERDDP_VALUE_BLOCK_WIDTH -ErrorAction SilentlyContinue
        } else {
            $env:FILTERDDP_FACTOR_BACKED_POLICY = "1"
            $env:FILTERDDP_BLOCKED_VALUE_RHS = "1"
            $env:FILTERDDP_VALUE_BLOCK_WIDTH = "128"
        }

        $tag = "${TagPrefix}_${mode}_T${Horizon}_r${rep}"
        $stdout = Join-Path $logDir "$tag.out.log"
        $stderr = Join-Path $logDir "$tag.err.log"
        $trace = Join-Path $resultDir "${tag}_trace.csv"
        Add-Content -LiteralPath $statusLog -Value "$(Get-Date -Format o) RUN_START mode=$mode rep=$rep T=$Horizon"
        # Start-Process on Windows can lose embedded quoting around absolute
        # paths containing spaces.  The process already starts at the repo
        # root, so use equivalent relative paths here.
        $timer = [Diagnostics.Stopwatch]::StartNew()
        Push-Location $repoNative
        try {
            & $julia --startup-file=no --project=envs/ddp2026 `
                ddp/examples/power_system/ieee123c_filterddp.jl `
                ieee2522C_1ph $Horizon solve 1> $stdout 2> $stderr
            $exitCode = $LASTEXITCODE
        } finally {
            Pop-Location
        }
        $timer.Stop()
        & $extractor -InputLog $stdout -OutputCsv $trace | Out-File (Join-Path $logDir "${tag}_trace.log")
        # Direct invocation avoids Start-Process corrupting paths containing
        # spaces. Peak RAM is established by the dedicated memory probes.
        $peakMiB = [double]::NaN
        $traceRelative = "ddp/results/network_filterddp/${tag}_trace.csv"
        Add-Content -LiteralPath $summary -Value ("{0},{1},{2},{3:F3},{4:F3},{5},{6}" -f `
            $mode, $rep, $Horizon, $timer.Elapsed.TotalSeconds, $peakMiB, $exitCode, $traceRelative)
        Add-Content -LiteralPath $statusLog -Value "$(Get-Date -Format o) RUN_END mode=$mode rep=$rep T=$Horizon exit=$exitCode wall_s=$($timer.Elapsed.TotalSeconds) peak_MiB=$peakMiB"
        if ($exitCode -ne 0) {
            throw "Benchmark failed: mode=$mode rep=$rep T=$Horizon"
        }
    }
}

Add-Content -LiteralPath $statusLog -Value "$(Get-Date -Format o) QUEUE_END T=$Horizon"
