param(
    [ValidateSet("large10kC_1ph")]
    [string]$System = "large10kC_1ph",
    [ValidateSet(96, 144, 192, 288)]
    [int]$Horizon = 96
)

$ErrorActionPreference = "Stop"
$repo = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$julia = Join-Path (Split-Path (Split-Path $repo -Parent) -Parent) ".tools\julia-1.12.6\bin\julia.exe"
if (-not (Test-Path -LiteralPath $julia)) {
    $julia = Join-Path (Split-Path (Split-Path (Split-Path $repo -Parent) -Parent) -Parent) ".tools\julia-1.12.6\bin\julia.exe"
}
$project = Join-Path $repo "envs\ddp2026"
$driver = Join-Path $repo "ddp\examples\power_system\centralized_ipopt_matched.jl"
$resultDir = Join-Path $repo "ddp\results\centralized_ipopt_matched_knee"
$logDir = Join-Path $repo "logs"
$tag = "${System}_T${Horizon}_periodic"
$stdout = Join-Path $logDir "matched_ipopt_knee_${tag}.out.log"
$stderr = Join-Path $logDir "matched_ipopt_knee_${tag}.err.log"
$ipoptLog = Join-Path $resultDir "${tag}_ipopt.log"
$statusLog = Join-Path $resultDir "${tag}_status.log"
New-Item -ItemType Directory -Force -Path $resultDir, $logDir | Out-Null

$env:REDUCED_PROFILE = "periodic"
$env:REDUCED_CB = "1e-3"
$env:TERMINAL_SOC_SOFT = "1"
$started = Get-Date
"STARTED=$($started.ToString('o'))" | Set-Content -LiteralPath $statusLog

$process = Start-Process -FilePath $julia `
    -ArgumentList @("--startup-file=no", "--project=`"$project`"", "`"$driver`"", $System, [string]$Horizon, "`"$ipoptLog`"") `
    -WorkingDirectory $repo -RedirectStandardOutput $stdout -RedirectStandardError $stderr `
    -WindowStyle Hidden -PassThru

$peakBytes = 0L
while (-not $process.HasExited) {
    try {
        $process.Refresh()
        $peakBytes = [Math]::Max($peakBytes, [int64]$process.WorkingSet64)
        "RUNNING pid=$($process.Id) peak_mib=$([Math]::Round($peakBytes / 1MB, 3)) elapsed_s=$([Math]::Round(((Get-Date)-$started).TotalSeconds, 1))" |
            Set-Content -LiteralPath $statusLog
    } catch { }
    Start-Sleep -Milliseconds 500
}
$process.WaitForExit()
$process.Refresh()
$elapsed = ((Get-Date) - $started).TotalSeconds
"FINISHED exit_code=$($process.ExitCode) peak_mib=$([Math]::Round($peakBytes / 1MB, 3)) elapsed_s=$([Math]::Round($elapsed, 3))" |
    Set-Content -LiteralPath $statusLog
exit $process.ExitCode
