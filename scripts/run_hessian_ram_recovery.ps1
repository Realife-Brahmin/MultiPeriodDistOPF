param(
    [string]$SourceWorktree = "C:\Users\Aryan Ritwajeet Jha\Documents\documents_general\MultiPeriodDistOPF\.claude\worktrees\new-session-a7b7a2"
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
$julia = "C:\Users\Aryan Ritwajeet Jha\Documents\documents_general\MultiPeriodDistOPF\.tools\julia-1.12.6\bin\julia.exe"
$outDir = Join-Path $repo "ddp\results\hessian_rewrites\ram_recovery"
$sourceData = Join-Path $SourceWorktree "ddp\results\network_filterddp"
$references = Import-Csv (Join-Path $SourceWorktree "ddp\results\matched_ipopt_race\matched_race.csv")
$summary = Join-Path $outDir "peak_ram.csv"
$status = Join-Path $outDir "status.log"
New-Item -ItemType Directory -Force -Path $outDir | Out-Null
if (-not (Test-Path $summary)) { Set-Content $summary "system,T,peak_working_set_gb,exit_code,completed_at" }

$cases = @(
    @{System="ieee123C_1ph";T=6}, @{System="ieee123C_1ph";T=12},
    @{System="ieee123C_1ph";T=24}, @{System="ieee123C_1ph";T=48},
    @{System="ieee123C_1ph";T=96}, @{System="ieee2522C_1ph";T=6},
    @{System="ieee2522C_1ph";T=12}, @{System="ieee2522C_1ph";T=24}
)

Push-Location $repo
try {
    foreach ($case in $cases) {
        $system=$case.System; $T=$case.T
        if ((Import-Csv $summary | Where-Object { $_.system -eq $system -and [int]$_.T -eq $T })) { continue }
        $ref=($references | Where-Object { $_.system -eq $system -and [int]$_.T -eq $T } | Select-Object -First 1).ipopt_objective
        Copy-Item (Join-Path $sourceData "network_data_${system}_T${T}_periodic.jls") (Join-Path $repo "ddp\results\network_filterddp\network_data_${system}_T${T}_periodic.jls") -Force
        $env:REDUCED_PROFILE="periodic"; $env:REDUCED_CB="1e-3"; $env:TERMINAL_SOC_SOFT="1"
        $env:FILTERDDP_SKIP_SOLUTION_WRITE="1"; $env:FILTERDDP_MAX_ITERATIONS="400"
        $env:FILTERDDP_DIAG_HESSIAN="1"; $env:FILTERDDP_DIAG_HESSIAN_FLOOR="1e-8"
        $env:FILTERDDP_DIRECT_DIAG_HESSIAN="1"; $env:FILTERDDP_TRIPLET_SECOND_DERIVATIVES="1"
        $env:FILTERDDP_CACHE_KKT_PATTERN="1"; $env:FILTERDDP_NEAR_OPT_REFERENCE=$ref
        $env:FILTERDDP_NEAR_OPT_GAP="0.005"; $env:FILTERDDP_NEAR_OPT_PRIMAL="1e-6"
        $env:FILTERDDP_FEASIBILITY_DIAGNOSTIC="0"; $env:FILTERDDP_TIMING_DIAGNOSTIC="0"
        Add-Content $status "$(Get-Date -Format o) START $system T=$T"
        $stdout=Join-Path $outDir "${system}_T${T}.out.log"; $stderr=Join-Path $outDir "${system}_T${T}.err.log"
        $args="--project=envs/ddp2026 ddp/examples/power_system/ieee123c_filterddp.jl $system $T solve"
        $p=Start-Process $julia -ArgumentList $args -WorkingDirectory $repo -RedirectStandardOutput $stdout -RedirectStandardError $stderr -PassThru -WindowStyle Hidden
        $peak=0.0
        while (-not $p.HasExited) {
            $p.Refresh(); $peak=[math]::Max($peak,$p.PeakWorkingSet64/1GB); Start-Sleep -Seconds 2
        }
        $p.WaitForExit()
        $p.Refresh(); $peak=[math]::Max($peak,$p.PeakWorkingSet64/1GB)
        $exitCode = $p.ExitCode
        if ($null -eq $exitCode) {
            $exitCode = if (Select-String -LiteralPath $stdout -Pattern 'solve complete: .* status=9' -Quiet) { 0 } else { 1 }
        }
        Add-Content $summary ("{0},{1},{2:F6},{3},{4}" -f $system,$T,$peak,$exitCode,(Get-Date -Format o))
        Add-Content $status ("{0} COMPLETE {1} T={2} peak_gb={3:F3} exit={4}" -f (Get-Date -Format o),$system,$T,$peak,$exitCode)
        if ($exitCode -ne 0) { throw "RAM recovery failed for $system T=$T" }
    }
} finally { Pop-Location }
