param(
    [string]$SourceWorktree = "C:\Users\Aryan Ritwajeet Jha\Documents\documents_general\MultiPeriodDistOPF\.claude\worktrees\new-session-a7b7a2"
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
$julia = "C:\Users\Aryan Ritwajeet Jha\Documents\documents_general\MultiPeriodDistOPF\.tools\julia-1.12.6\bin\julia.exe"
$resultDir = Join-Path $repo "ddp\results\hessian_rewrites\overnight"
$sourceData = Join-Path $SourceWorktree "ddp\results\network_filterddp"
$referenceCsv = Join-Path $SourceWorktree "ddp\results\matched_ipopt_race\matched_race.csv"
$statusLog = Join-Path $resultDir "queue.status.log"
$summaryCsv = Join-Path $resultDir "queue_summary.csv"
$lockFile = Join-Path $resultDir "queue.lock"

New-Item -ItemType Directory -Force -Path $resultDir | Out-Null
Set-Content -LiteralPath $lockFile -Value $PID
if (-not (Test-Path -LiteralPath $summaryCsv)) {
    Set-Content -LiteralPath $summaryCsv -Value "system,T,reference,primal_threshold,state,trace,completed_at"
}

function Write-Status([string]$message) {
    $line = "$(Get-Date -Format o) $message"
    Add-Content -LiteralPath $statusLog -Value $line
    Write-Output $line
}

$cases = @(
    @{ System = "ieee123C_1ph"; T = 6;  Primal = "1e-6"; Factor = "0" },
    @{ System = "ieee123C_1ph"; T = 12; Primal = "1e-6"; Factor = "0" },
    @{ System = "ieee123C_1ph"; T = 24; Primal = "1e-6"; Factor = "0" },
    @{ System = "ieee123C_1ph"; T = 48; Primal = "1e-6"; Factor = "0" },
    @{ System = "ieee123C_1ph"; T = 96; Primal = "1e-6"; Factor = "0" },
    @{ System = "ieee2522C_1ph"; T = 6;  Primal = "1e-6"; Factor = "0" },
    @{ System = "ieee2522C_1ph"; T = 12; Primal = "1e-6"; Factor = "0" },
    @{ System = "ieee2522C_1ph"; T = 48; Primal = "1e-6"; Factor = "0" },
    @{ System = "ieee2522C_1ph"; T = 96; Primal = "1e-6"; Factor = "0" },
    @{ System = "large10kC_1ph"; T = 24; Primal = "1e-4"; Factor = "1" }
)

$references = Import-Csv -LiteralPath $referenceCsv
Push-Location $repo
try {
    foreach ($case in $cases) {
        $system = $case.System
        $T = $case.T
        $tag = "combined_${system}_T${T}"
        $trace = Join-Path $resultDir "${tag}_trace.csv"
        $log = Join-Path $resultDir "${tag}.out.log"
        $err = Join-Path $resultDir "${tag}.err.log"
        if (Test-Path -LiteralPath $trace) {
            Write-Status "SKIP $system T=$T trace already exists"
            continue
        }

        $refRow = $references | Where-Object { $_.system -eq $system -and [int]$_.T -eq $T } | Select-Object -First 1
        if ($null -eq $refRow) { throw "No matched IPOPT reference for $system T=$T" }
        $reference = $refRow.ipopt_objective

        $dataName = "network_data_${system}_T${T}_periodic.jls"
        $source = Join-Path $sourceData $dataName
        $destination = Join-Path $repo "ddp\results\network_filterddp\$dataName"
        if (-not (Test-Path -LiteralPath $source)) { throw "Missing exported data: $source" }
        Copy-Item -LiteralPath $source -Destination $destination -Force

        $env:REDUCED_PROFILE = "periodic"
        $env:REDUCED_CB = "1e-3"
        $env:TERMINAL_SOC_SOFT = "1"
        $env:FILTERDDP_SKIP_SOLUTION_WRITE = "1"
        $env:FILTERDDP_TIMING_DIAGNOSTIC = "1"
        $env:FILTERDDP_MAX_ITERATIONS = "400"
        $env:FILTERDDP_DIAG_HESSIAN = "1"
        $env:FILTERDDP_DIRECT_DIAG_HESSIAN = "1"
        $env:FILTERDDP_TRIPLET_SECOND_DERIVATIVES = "1"
        $env:FILTERDDP_CACHE_KKT_PATTERN = "1"
        $env:FILTERDDP_DIAG_HESSIAN_FLOOR = "1e-8"
        $env:FILTERDDP_FACTOR_BACKED_POLICY = $case.Factor
        $env:FILTERDDP_NEAR_OPT_REFERENCE = $reference
        $env:FILTERDDP_NEAR_OPT_GAP = "0.005"
        $env:FILTERDDP_NEAR_OPT_PRIMAL = $case.Primal
        $env:FILTERDDP_FEASIBILITY_DIAGNOSTIC = "1"

        Write-Status "START $system T=$T reference=$reference primal=$($case.Primal) factor=$($case.Factor)"
        try {
            & $julia --project=envs/ddp2026 ddp/examples/power_system/ieee123c_filterddp.jl $system $T solve 2>&1 |
                Tee-Object -LiteralPath $log
            if ($LASTEXITCODE -ne 0) { throw "Julia exited with code $LASTEXITCODE" }
        } catch {
            $_ | Out-File -LiteralPath $err -Append
            Write-Status "FAILED $system T=$T error=$($_.Exception.Message)"
            throw
        }

        & $julia --project=envs/ddp2026 ddp/examples/power_system/extract_filterddp_feasibility_trace.jl $log $trace $reference |
            Tee-Object -LiteralPath $err -Append
        if ($LASTEXITCODE -ne 0 -or -not (Test-Path -LiteralPath $trace)) {
            throw "Trace extraction failed for $system T=$T"
        }
        $last = Import-Csv -LiteralPath $trace | Select-Object -Last 1
        if ($null -eq $last) { throw "Empty trace for $system T=$T" }
        $completedAt = Get-Date -Format o
        Add-Content -LiteralPath $summaryCsv -Value "$system,$T,$reference,$($case.Primal),complete,$([IO.Path]::GetFileName($trace)),$completedAt"
        Write-Status "COMPLETE $system T=$T iteration=$($last.iteration) elapsed=$($last.elapsed_algorithm_s) objective=$($last.objective) gap=$($last.objective_rel_gap) primal=$($last.primal_inf)"

        git -c safe.directory="$repo" add -- $trace $summaryCsv $statusLog
        git -c safe.directory="$repo" commit -m "record combined Hessian rewrite $system T=$T" -- $trace $summaryCsv $statusLog
        git -c safe.directory="$repo" push origin ddp-hessian-rewrites-sep20
    }
    Write-Status "QUEUE COMPLETE"
} finally {
    Pop-Location
    Remove-Item -LiteralPath $lockFile -Force -ErrorAction SilentlyContinue
}
