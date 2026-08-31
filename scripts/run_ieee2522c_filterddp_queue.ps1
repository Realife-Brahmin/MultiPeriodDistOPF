param(
    [int[]]$Horizons = @(48, 96)
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
$julia = Join-Path $repo ".tools\julia-1.12.6\bin\julia.exe"
$depot = Join-Path $repo ".tools\julia_depot"
$exportProject = Join-Path $repo "envs\tadmm"
$solveProject = Join-Path $repo "envs\ddp2026"
$exportDriver = Join-Path $repo "ddp\examples\power_system\export_ieee123c_data.jl"
$solveDriver = Join-Path $repo "ddp\examples\power_system\ieee123c_filterddp.jl"
$extractor = Join-Path $repo "scripts\extract_filterddp_trace.ps1"
$resultDir = Join-Path $repo "ddp\results\network_filterddp"
$logDir = Join-Path $repo "logs"
$statusLog = Join-Path $logDir "ieee2522c_filterddp_queue.status.log"
$system = "ieee2522C_1ph"

New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$env:JULIA_DEPOT_PATH = $depot
$env:FILTERDDP_MAX_ITERATIONS = "200"

foreach ($horizon in $Horizons) {
    $dataFile = Join-Path $resultDir "network_data_${system}_T${horizon}.jls"
    if (-not (Test-Path -LiteralPath $dataFile)) {
        $exportOut = Join-Path $logDir "ieee2522c_export_T${horizon}.out.log"
        $exportErr = Join-Path $logDir "ieee2522c_export_T${horizon}.err.log"
        $started = Get-Date -Format "yyyy-MM-dd HH:mm:ss zzz"
        Add-Content -LiteralPath $statusLog -Value "$started EXPORT_START T=$horizon"
        $exportArgs = "--startup-file=no --project=`"$exportProject`" `"$exportDriver`" $system $horizon"
        $exportProcess = Start-Process -FilePath $julia -ArgumentList $exportArgs `
            -WorkingDirectory $repo -WindowStyle Hidden -Wait -PassThru `
            -RedirectStandardOutput $exportOut -RedirectStandardError $exportErr
        $finished = Get-Date -Format "yyyy-MM-dd HH:mm:ss zzz"
        Add-Content -LiteralPath $statusLog -Value "$finished EXPORT_END T=$horizon exit=$($exportProcess.ExitCode)"
        if ($exportProcess.ExitCode -ne 0) { break }
    }

    $stdout = Join-Path $logDir "ieee2522c_filterddp_T${horizon}.out.log"
    $stderr = Join-Path $logDir "ieee2522c_filterddp_T${horizon}.err.log"
    $trace = Join-Path $resultDir "ieee2522c_filterddp_T${horizon}_trace.csv"
    $started = Get-Date -Format "yyyy-MM-dd HH:mm:ss zzz"
    Add-Content -LiteralPath $statusLog -Value "$started SOLVE_START T=$horizon"
    $solveArgs = "--startup-file=no --project=`"$solveProject`" `"$solveDriver`" $system $horizon solve"
    $solveProcess = Start-Process -FilePath $julia -ArgumentList $solveArgs `
        -WorkingDirectory $repo -WindowStyle Hidden -Wait -PassThru `
        -RedirectStandardOutput $stdout -RedirectStandardError $stderr

    & $extractor -InputLog $stdout -OutputCsv $trace | Out-File `
        -LiteralPath (Join-Path $logDir "ieee2522c_trace_T${horizon}.log")
    $finished = Get-Date -Format "yyyy-MM-dd HH:mm:ss zzz"
    Add-Content -LiteralPath $statusLog -Value "$finished SOLVE_END T=$horizon exit=$($solveProcess.ExitCode) trace=$trace"
    if ($solveProcess.ExitCode -ne 0) { break }
}
