param(
    [ValidateSet("ieee2522C_1ph", "large10kC_1ph")]
    [string]$System = "ieee2522C_1ph",
    [ValidateSet(3, 12)]
    [int]$Horizon = 12,
    [int[]]$ThreadCounts = @(1, 2, 4, 8)
)

$ErrorActionPreference = "Stop"
$repo = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$runtimeRoot = $repo
$julia = $null
while ($runtimeRoot) {
    $candidate = Join-Path $runtimeRoot ".tools\julia-1.12.6\bin\julia.exe"
    if (Test-Path -LiteralPath $candidate) { $julia = $candidate; break }
    $parent = Split-Path $runtimeRoot -Parent
    if ($parent -eq $runtimeRoot) { break }
    $runtimeRoot = $parent
}
if (-not $julia) { throw "Repo-local Julia runtime not found" }

$probeDir = Join-Path $repo "tmp\ipopt_thread_probe"
New-Item -ItemType Directory -Force -Path $probeDir | Out-Null
$script = Join-Path $repo "ddp\examples\power_system\ieee123c_ipopt_matrix_diagnostic.jl"
$project = Join-Path $repo "envs\ddp2026"

foreach ($threads in $ThreadCounts) {
    $tag = "${System}_T${Horizon}_threads${threads}"
    $env:OMP_NUM_THREADS = [string]$threads
    $env:OPENBLAS_NUM_THREADS = [string]$threads
    $env:MKL_NUM_THREADS = [string]$threads
    $env:JULIA_NUM_THREADS = "1"
    $env:IPOPT_DIAGNOSTIC_OUTPUT = Join-Path $probeDir "${tag}_ipopt.log"
    $stdout = Join-Path $probeDir "${tag}.out.log"
    $stderr = Join-Path $probeDir "${tag}.err.log"
    Write-Output "START $tag"
    $process = Start-Process -FilePath $julia -ArgumentList @(
        "--startup-file=no", "--threads=1", "--project=`"$project`"", "`"$script`"",
        $System, [string]$Horizon
    ) -RedirectStandardOutput $stdout -RedirectStandardError $stderr `
      -WindowStyle Hidden -PassThru
    $process.WaitForExit()
    if ($process.ExitCode -ne 0) {
        throw "$tag failed; inspect $stderr"
    }
    Write-Output "DONE $tag"
}
