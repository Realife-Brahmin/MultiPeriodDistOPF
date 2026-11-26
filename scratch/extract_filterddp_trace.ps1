param(
    [Parameter(Mandatory = $true)]
    [string]$InputLog,
    [Parameter(Mandatory = $true)]
    [string]$OutputCsv
)

$pattern = '^\s*(\d+)\s+([+-]?[\d.]+e[+-]\d+)\s+([+-]?[\d.]+e[+-]\d+)\s+([+-]?[\d.]+e[+-]\d+)\s+([+-]?[\d.]+e[+-]\d+)\s+([+-]?[\d.]+)\s+(\S+)\s+([+-]?[\d.]+e[+-]\d+)\s+(\d+)\s*$'
$rows = foreach ($line in Get-Content -LiteralPath $InputLog) {
    if ($line -match $pattern) {
        [pscustomobject]@{
            iteration = [int]$Matches[1]
            objective = [double]$Matches[2]
            primal_inf = [double]$Matches[3]
            dual_inf = [double]$Matches[4]
            complementarity_inf = [double]$Matches[5]
            log10_barrier = [double]$Matches[6]
            log10_regularization = if ($Matches[7] -eq '-') { $null } else { [double]$Matches[7] }
            step_size = [double]$Matches[8]
            line_search_backtracks = [int]$Matches[9]
        }
    }
}

$rows | Export-Csv -LiteralPath $OutputCsv -NoTypeInformation
Write-Output "exported=$OutputCsv rows=$($rows.Count)"
