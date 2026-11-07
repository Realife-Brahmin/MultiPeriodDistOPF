param([string]$SystemName='ieee123C_1ph',[int]$Horizon=6,[string]$RunLabel='pilot')
$ErrorActionPreference='Stop'
$repoRoot=(Resolve-Path (Join-Path $PSScriptRoot '../../../..')).Path
Set-Location -LiteralPath $repoRoot
$juliaExe='C:\Users\Aryan Ritwajeet Jha\Documents\documents_general\MultiPeriodDistOPF\.tools\julia-1.12.6\bin\julia.exe'
$env:JULIA_DEPOT_PATH='C:\Users\Aryan Ritwajeet Jha\Documents\documents_general\MultiPeriodDistOPF\.tools\julia_depot'
$env:JULIA_NUM_THREADS='1'
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PROFILE_PERIODIC='1'
$env:REDUCED_PROFILE='periodic'
$env:REDUCED_CB='system'
$env:TERMINAL_SOC_SOFT='0'
$env:OPENDSS_ZERO_LINE_C='1'
$outDir=Join-Path $repoRoot "problems/mpopf/tadmm/review_pilot/results/${SystemName}_T${Horizon}_${RunLabel}"
if(Test-Path -LiteralPath $outDir){throw "Refusing to overwrite $outDir"}
New-Item -ItemType Directory -Path $outDir | Out-Null
python (Join-Path $PSScriptRoot 'capture.py') $outDir $SystemName
if($LASTEXITCODE -ne 0){throw 'Source capture failed'}
$projectArg='--project=problems/mpopf/tadmm/review_pilot'
& $juliaExe --startup-file=no $projectArg problems/mpopf/spar_ddp/examples/power_system/export_ieee123c_data.jl $SystemName $Horizon *> (Join-Path $outDir 'export.log')
if($LASTEXITCODE -ne 0){throw "Export failed: $outDir"}
$argLine="--startup-file=no $projectArg problems/mpopf/tadmm/review_pilot/run.jl $SystemName $Horizon `"$outDir`""
$child=Start-Process -FilePath $juliaExe -ArgumentList $argLine -WorkingDirectory $repoRoot -WindowStyle Hidden -PassThru -RedirectStandardOutput (Join-Path $outDir 'solve.log') -RedirectStandardError (Join-Path $outDir 'solve.err')
$memPath=Join-Path $outDir 'process_samples.csv'
'utc,working_set_bytes,private_bytes,cpu_seconds' | Set-Content -LiteralPath $memPath
do {
    $child.Refresh()
    if(-not $child.HasExited){
        '{0},{1},{2},{3}' -f ([DateTime]::UtcNow.ToString('o')),$child.WorkingSet64,$child.PrivateMemorySize64,$child.TotalProcessorTime.TotalSeconds | Add-Content -LiteralPath $memPath
    }
    Start-Sleep -Milliseconds 500
} while(-not $child.HasExited)
$child.WaitForExit()
if($child.ExitCode -ne 0){throw "Optimization failed (exit $($child.ExitCode)): $outDir"}
foreach($method in @('central','tadmm')){
    & $juliaExe --startup-file=no $projectArg problems/mpopf/spar_ddp/examples/power_system/opendss_check.jl $SystemName $Horizon (Join-Path $outDir "${method}_dispatch.jls") $method *> (Join-Path $outDir "${method}_opendss.log")
    if($LASTEXITCODE -ne 0){throw "OpenDSS check failed: $method $outDir"}
}
python (Join-Path $PSScriptRoot 'validate.py') $outDir
if($LASTEXITCODE -ne 0){throw "Evidence validation failed: $outDir"}
Write-Output "CASE_FINISHED $outDir"
