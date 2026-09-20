param(
    [Parameter(Mandatory = $true)]
    [ValidateSet("ieee123C_1ph", "ieee2522C_1ph", "large10kC_1ph")]
    [string]$System,

    [Parameter(Mandatory = $true)]
    [ValidateSet(3, 6, 12, 24, 48, 96, 144, 192, 288, 384, 576)]
    [int]$Horizon
)

$ErrorActionPreference = "Stop"
$repo = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$runtimeRoot = $repo
$julia = $null
while ($runtimeRoot) {
    $candidate = Join-Path $runtimeRoot ".tools\julia-1.12.6\bin\julia.exe"
    if (Test-Path -LiteralPath $candidate) {
        $julia = $candidate
        break
    }
    $parent = Split-Path $runtimeRoot -Parent
    if ($parent -eq $runtimeRoot) { break }
    $runtimeRoot = $parent
}
$envPath = Join-Path $repo "envs\tadmm"
$runDir = Join-Path $envPath "root_level"
$generatedDir = Join-Path $envPath "processedData\${System}_T${Horizon}"
$resultsDir = Join-Path $repo "ddp\results\centralized_ipopt"
$rawDir = Join-Path $resultsDir "raw"
$localLogDir = Join-Path $repo "logs"
$tag = "${System}_T${Horizon}"
$stdout = Join-Path $localLogDir "centralized_ipopt_${tag}.out.log"
$stderr = Join-Path $localLogDir "centralized_ipopt_${tag}.err.log"

if (-not (Test-Path -LiteralPath $julia)) {
    throw "Julia not found at $julia. Install the repo-local runtime first."
}
New-Item -ItemType Directory -Force -Path $resultsDir, $rawDir, $localLogDir | Out-Null

$env:SYSTEM_OVERRIDE = $System
$env:T_OVERRIDE = [string]$Horizon
$env:USE_GUROBI_OVERRIDE = "false"
$env:USE_GUROBI_FOR_BF_OVERRIDE = "false"
$env:PRESERVE_IPOPT_FULL_LOG = "1"
$env:IPOPT_TIMING_STATISTICS = "1"

$startedUtc = [DateTime]::UtcNow
$process = Start-Process -FilePath $julia `
    -ArgumentList @("--startup-file=no", "--project=`"$envPath`"", "run_bf.jl") `
    -WorkingDirectory $runDir -RedirectStandardOutput $stdout `
    -RedirectStandardError $stderr -WindowStyle Hidden -PassThru

$peakBytes = 0L
while (-not $process.HasExited) {
    try {
        $process.Refresh()
        $peakBytes = [Math]::Max($peakBytes, [int64]$process.WorkingSet64)
    } catch { }
    Start-Sleep -Milliseconds 500
}
$process.WaitForExit()
$process.Refresh()
$exitCode = $process.ExitCode
if ($null -ne $exitCode -and $exitCode -ne 0) {
    throw "IPOPT run failed with exit code $($process.ExitCode); preserve $stdout and $stderr."
}

$summaryPath = Join-Path $generatedDir "results_socp_bf.txt"
$ipoptPath = Join-Path $generatedDir "ipopt_bf.log"
if (-not (Test-Path -LiteralPath $summaryPath) -or -not (Test-Path -LiteralPath $ipoptPath)) {
    throw "Run completed without the expected summary or IPOPT log."
}
$timingDeadline = (Get-Date).AddSeconds(30)
while ((Get-Date) -lt $timingDeadline) {
    $probe = Get-Content -LiteralPath $ipoptPath -Raw -ErrorAction SilentlyContinue
    if ($probe -match 'Timing Statistics:') { break }
    Start-Sleep -Milliseconds 500
}
$summary = Get-Content -LiteralPath $summaryPath -Raw
$ipopt = Get-Content -LiteralPath $ipoptPath -Raw

function Match-Value([string]$Text, [string]$Pattern, [string]$Name) {
    $match = [regex]::Match($Text, $Pattern, [Text.RegularExpressions.RegexOptions]::Multiline)
    if (-not $match.Success) { throw "Could not parse $Name" }
    return $match.Groups[1].Value.Trim()
}

function Match-Timing-Wall([string]$Text, [string]$Name) {
    $escaped = [regex]::Escape($Name)
    $match = [regex]::Match($Text, "(?m)^\s*$escaped\.+:\s+[0-9.]+\s+\(sys:\s+[0-9.]+\s+wall:\s+([0-9.]+)\)")
    if (-not $match.Success) { return [double]::NaN }
    return [double]$match.Groups[1].Value
}

$overallWall = Match-Timing-Wall $ipopt 'OverallAlgorithm'
$searchWall = Match-Timing-Wall $ipopt 'ComputeSearchDirection'
$pdWall = Match-Timing-Wall $ipopt 'PDSystemSolverTotal'
$augWall = Match-Timing-Wall $ipopt 'StdAugSystemSolverMultiSolve'
$factorWall = Match-Timing-Wall $ipopt 'LinearSystemFactorization'
$backsolveWall = Match-Timing-Wall $ipopt 'LinearSystemBackSolve'
$functionWall = Match-Timing-Wall $ipopt 'Function Evaluations'
$inferredFactorWall = if (-not [double]::IsNaN($augWall) -and -not [double]::IsNaN($backsolveWall)) {
    [Math]::Max(0.0, $augWall - $backsolveWall)
} else { [double]::NaN }

$status = Match-Value $summary '^Status:\s*(\S+)' 'status'
$validation = Match-Value $summary '^Status:\s*(FEASIBLE|INFEASIBLE|VALIDATION SKIPPED)' 'validation'
$row = [ordered]@{
    system = $System
    T = $Horizon
    delta_t_h = 24.0 / $Horizon
    status = $status
    validated = ($validation -eq "FEASIBLE")
    objective_usd = [double](Match-Value $ipopt '^Objective.*?\s([-+0-9.eE]+)\s*$' 'objective')
    iterations = [int](Match-Value $ipopt '^Number of Iterations\.*:\s*(\d+)' 'iterations')
    ipopt_reported_s = [double](Match-Value $ipopt '^Total seconds in IPOPT(?: \(w/o function evaluations\))?\s*=\s*([0-9.]+)' 'IPOPT time')
    jump_solve_time_s = [double](Match-Value $summary '^Solver time:\s*([0-9.]+)' 'JuMP solve time')
    solve_wall_s = [double](Match-Value $summary '^Wall-clock time:\s*([0-9.]+)' 'solve wall time')
    peak_working_set_mib = [Math]::Round($peakBytes / 1MB, 3)
    variables = [int](Match-Value $ipopt '^Total number of variables\.*:\s*(\d+)' 'variables')
    equality_constraints = [int](Match-Value $ipopt '^Total number of equality constraints\.*:\s*(\d+)' 'equalities')
    inequality_constraints = [int](Match-Value $ipopt '^Total number of inequality constraints\.*:\s*(\d+)' 'inequalities')
    equality_jacobian_nnz = [int](Match-Value $ipopt '^Number of nonzeros in equality constraint Jacobian\.*:\s*(\d+)' 'equality Jacobian nnz')
    inequality_jacobian_nnz = [int](Match-Value $ipopt '^Number of nonzeros in inequality constraint Jacobian\.*:\s*(\d+)' 'inequality Jacobian nnz')
    lagrangian_hessian_nnz = [int](Match-Value $ipopt '^Number of nonzeros in Lagrangian Hessian\.*:\s*(\d+)' 'Hessian nnz')
    linear_solver = Match-Value $ipopt '^This is Ipopt.*linear solver\s+([^\r\n]+)' 'linear solver'
    started_utc = $startedUtc.ToString("o")
    raw_log = "ddp/results/centralized_ipopt/raw/${tag}_ipopt.log"
    timing_overall_wall_s = $overallWall
    timing_search_direction_wall_s = $searchWall
    timing_pd_system_wall_s = $pdWall
    timing_augmented_system_wall_s = $augWall
    timing_factorization_reported_wall_s = $factorWall
    timing_backsolve_wall_s = $backsolveWall
    timing_factorization_inferred_wall_s = $inferredFactorWall
    timing_function_evaluations_wall_s = $functionWall
}

$rawTarget = Join-Path $rawDir "${tag}_ipopt.log"
Copy-Item -LiteralPath $ipoptPath -Destination $rawTarget -Force
Copy-Item -LiteralPath $summaryPath -Destination (Join-Path $rawDir "${tag}_summary.txt") -Force

$isKneeCase = ($System -eq "ieee2522C_1ph" -and $Horizon -gt 288) -or
              ($System -eq "large10kC_1ph" -and $Horizon -gt 48)
$csvName = if ($isKneeCase) { "centralized_ipopt_knee.csv" } else { "centralized_ipopt_timing.csv" }
$csvPath = Join-Path $resultsDir $csvName
$rows = if (Test-Path -LiteralPath $csvPath) { @(Import-Csv -LiteralPath $csvPath) } else { @() }
$rows = @($rows | Where-Object { $_.system -ne $System -or [int]$_.T -ne $Horizon })
$rows += [pscustomobject]$row
$systemOrder = @{ ieee123C_1ph = 1; ieee2522C_1ph = 2; large10kC_1ph = 3 }
$rows | Sort-Object @{ Expression = { $systemOrder[$_.system] } }, @{ Expression = { [int]$_.T } } |
    Export-Csv -LiteralPath $csvPath -NoTypeInformation

$row | Format-List
Write-Output "Updated $csvPath"
Write-Output "Next: validate the row, update both repository docs, then commit and push each repository before starting another case."
