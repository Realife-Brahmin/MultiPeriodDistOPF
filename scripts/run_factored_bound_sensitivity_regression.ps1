param(
    [string[]]$Systems = @("ieee123C_1ph", "ieee2522C_1ph", "large10kC_1ph"),
    [int]$Horizon = 3,
    [double]$Tolerance = 1e-7,
    [int]$MaxIterations = 200
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
$julia = Join-Path $repo ".tools\julia-1.12.6\bin\julia.exe"
$depot = Join-Path $repo ".tools\julia_depot"
$project = Join-Path $repo "envs\ddp2026"
$driver = Join-Path $repo "ddp\examples\power_system\ieee123c_filterddp.jl"
$extractor = Join-Path $repo "scripts\extract_filterddp_trace.ps1"
$resultDir = Join-Path $repo "ddp\results\network_filterddp"
$logDir = Join-Path $repo "logs"
$statusLog = Join-Path $logDir "factored_bound_sensitivity_regression.status.log"

New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$env:JULIA_DEPOT_PATH = $depot
$env:FILTERDDP_MAX_ITERATIONS = $MaxIterations.ToString()
$env:FILTERDDP_OPTIMALITY_TOLERANCE = $Tolerance.ToString("R", [Globalization.CultureInfo]::InvariantCulture)
Remove-Item Env:FILTERDDP_MEMORY_DIAGNOSTIC -ErrorAction SilentlyContinue
Remove-Item Env:FILTERDDP_TIMING_DIAGNOSTIC -ErrorAction SilentlyContinue

foreach ($system in $Systems) {
    $tag = "factored_bounds_${system}_T${Horizon}"
    $stdout = Join-Path $logDir "${tag}.out.log"
    $stderr = Join-Path $logDir "${tag}.err.log"
    $trace = Join-Path $resultDir "${tag}_trace.csv"
    $started = Get-Date -Format "yyyy-MM-dd HH:mm:ss zzz"
    Set-Content -LiteralPath $stdout -Value ""
    Set-Content -LiteralPath $stderr -Value ""
    Add-Content -LiteralPath $statusLog -Value "$started SOLVE_START system=$system T=$Horizon tolerance=$Tolerance"

    $arguments = "--startup-file=no --project=`"$project`" `"$driver`" $system $Horizon solve"
    $process = Start-Process -FilePath $julia -ArgumentList $arguments `
        -WorkingDirectory $repo -WindowStyle Hidden -Wait -PassThru `
        -RedirectStandardOutput $stdout -RedirectStandardError $stderr

    & $extractor -InputLog $stdout -OutputCsv $trace | Out-File `
        -LiteralPath (Join-Path $logDir "${tag}_trace.log")
    $finished = Get-Date -Format "yyyy-MM-dd HH:mm:ss zzz"
    Add-Content -LiteralPath $statusLog -Value "$finished SOLVE_END system=$system T=$Horizon exit=$($process.ExitCode) trace=$trace"
    if ($process.ExitCode -ne 0) {
        throw "FilterDDP regression failed for $system with exit code $($process.ExitCode)"
    }
}
