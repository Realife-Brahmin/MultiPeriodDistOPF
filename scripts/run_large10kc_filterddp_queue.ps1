param(
    [int[]]$Horizons = @(6, 12)
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
$julia = Join-Path $repo ".tools\julia-1.12.6\bin\julia.exe"
$project = Join-Path $repo "envs\ddp2026"
$driver = Join-Path $repo "ddp\examples\power_system\ieee123c_filterddp.jl"
$logDir = Join-Path $repo "logs"
$statusLog = Join-Path $logDir "large10kc_filterddp_queue.status.log"

$env:JULIA_DEPOT_PATH = Join-Path $repo ".tools\julia_depot"

foreach ($horizon in $Horizons) {
    $stdout = Join-Path $logDir "large10kc_filterddp_T${horizon}.out.log"
    $stderr = Join-Path $logDir "large10kc_filterddp_T${horizon}.err.log"
    $started = Get-Date -Format "yyyy-MM-dd HH:mm:ss zzz"
    Add-Content -LiteralPath $statusLog -Value "$started START T=$horizon"

    $arguments = @(
        "--startup-file=no",
        "--project=`"$project`"",
        "`"$driver`"",
        "large10kC_1ph",
        "$horizon",
        "solve"
    )
    $process = Start-Process -FilePath $julia -ArgumentList $arguments `
        -WorkingDirectory $repo -WindowStyle Hidden -Wait -PassThru `
        -RedirectStandardOutput $stdout -RedirectStandardError $stderr

    $finished = Get-Date -Format "yyyy-MM-dd HH:mm:ss zzz"
    Add-Content -LiteralPath $statusLog -Value "$finished END T=$horizon exit=$($process.ExitCode)"
    if ($process.ExitCode -ne 0) {
        break
    }
}
