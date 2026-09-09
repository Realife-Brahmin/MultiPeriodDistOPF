param(
    [Parameter(Mandatory = $true)]
    [ValidateSet("ieee123C_1ph", "ieee2522C_1ph", "large10kC_1ph")]
    [string]$System,

    [Parameter(Mandatory = $true)]
    [ValidateSet(6, 12, 24, 48, 96, 144)]
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
$summary = Get-Content -LiteralPath $summaryPath -Raw
$ipopt = Get-Content -LiteralPath $ipoptPath -Raw

function Match-Value([string]$Text, [string]$Pattern, [string]$Name) {
    $match = [regex]::Match($Text, $Pattern, [Text.RegularExpressions.RegexOptions]::Multiline)
    if (-not $match.Success) { throw "Could not parse $Name" }
    return $match.Groups[1].Value.Trim()
}

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
    ipopt_reported_s = [double](Match-Value $ipopt '^Total seconds in IPOPT\s*=\s*([0-9.]+)' 'IPOPT time')
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
}

$rawTarget = Join-Path $rawDir "${tag}_ipopt.log"
Copy-Item -LiteralPath $ipoptPath -Destination $rawTarget -Force
Copy-Item -LiteralPath $summaryPath -Destination (Join-Path $rawDir "${tag}_summary.txt") -Force

$csvPath = Join-Path $resultsDir "centralized_ipopt_timing.csv"
$rows = if (Test-Path -LiteralPath $csvPath) { @(Import-Csv -LiteralPath $csvPath) } else { @() }
$rows = @($rows | Where-Object { $_.system -ne $System -or [int]$_.T -ne $Horizon })
$rows += [pscustomobject]$row
$systemOrder = @{ ieee123C_1ph = 1; ieee2522C_1ph = 2; large10kC_1ph = 3 }
$rows | Sort-Object @{ Expression = { $systemOrder[$_.system] } }, @{ Expression = { [int]$_.T } } |
    Export-Csv -LiteralPath $csvPath -NoTypeInformation

$row | Format-List
Write-Output "Updated $csvPath"
Write-Output "Next: validate the row, update both repository docs, then commit and push each repository before starting another case."
