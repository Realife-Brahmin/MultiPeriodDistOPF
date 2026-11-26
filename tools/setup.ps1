param(
    [ValidateSet("spar-ddp", "tadmm", "msopf")]
    [string]$Project = "spar-ddp",
    [string]$JuliaExecutable = "julia",
    [string]$JuliaChannel = "",
    [switch]$CheckOnly,
    [switch]$Precompile
)
$ErrorActionPreference = "Stop"
$setupRoot = Split-Path -Parent $PSScriptRoot
if (-not (Get-Command $JuliaExecutable -ErrorAction SilentlyContinue)) {
    throw "Julia not found. Install Julia/Juliaup: https://julialang.org/downloads/"
}
$setupArgs = @()
if ($JuliaChannel) { $setupArgs += "+$JuliaChannel" }
$setupArgs += @("--startup-file=no", (Join-Path $setupRoot "tools/setup.jl"), $Project)
if ($CheckOnly) { $setupArgs += "--check-only" }
if ($Precompile) { $setupArgs += "--precompile" }
# A script file plus an argument array avoids PowerShell's Julia -e quoting pitfalls.
& $JuliaExecutable @setupArgs
if ($LASTEXITCODE -ne 0) { throw "Setup failed (Julia exit code $LASTEXITCODE)" }
