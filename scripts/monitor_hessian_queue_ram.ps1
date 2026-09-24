param(
    [string]$Worktree = "C:\Users\Aryan Ritwajeet Jha\Documents\documents_general\MultiPeriodDistOPF\.codex\worktrees\hessian-rewrites-sep20",
    [int]$IntervalSeconds = 2
)

$ErrorActionPreference = "SilentlyContinue"
$resultDir = Join-Path $Worktree "ddp\results\hessian_rewrites\overnight"
$lockFile = Join-Path $resultDir "queue.lock"
$samples = Join-Path $resultDir "ram_samples.csv"
$peaks = Join-Path $resultDir "ram_peaks.csv"

if (-not (Test-Path -LiteralPath $samples)) {
    Set-Content -LiteralPath $samples -Value "timestamp,pid,start_time,working_set_gb,peak_working_set_gb"
}

$observed = @{}
while (Test-Path -LiteralPath $lockFile) {
    foreach ($process in Get-Process julia -ErrorAction SilentlyContinue) {
        $key = "$($process.Id)|$($process.StartTime.ToString('o'))"
        $working = $process.WorkingSet64 / 1GB
        $peak = $process.PeakWorkingSet64 / 1GB
        if (-not $observed.ContainsKey($key) -or $peak -gt $observed[$key].PeakGB) {
            $observed[$key] = [pscustomobject]@{
                PID = $process.Id
                StartTime = $process.StartTime.ToString('o')
                PeakGB = $peak
            }
        }
        Add-Content -LiteralPath $samples -Value ("{0},{1},{2},{3:F6},{4:F6}" -f (Get-Date -Format o), $process.Id, $process.StartTime.ToString('o'), $working, $peak)
    }
    Start-Sleep -Seconds $IntervalSeconds
}

Set-Content -LiteralPath $peaks -Value "pid,start_time,peak_working_set_gb"
foreach ($row in $observed.Values | Sort-Object StartTime) {
    Add-Content -LiteralPath $peaks -Value ("{0},{1},{2:F6}" -f $row.PID, $row.StartTime, $row.PeakGB)
}
