param(
    [Parameter(Mandatory = $true)][string]$InputLog,
    [Parameter(Mandatory = $true)][string]$OutputCsv
)

$ErrorActionPreference = "Stop"
$culture = [Globalization.CultureInfo]::InvariantCulture
$number = '[-+0-9.eE]+'
$stagePattern = "FILTERDDP_TIMING iteration=(\d+) barrier_iteration=(\d+) stage=(\d+) derivative_s=($number) first_order_s=($number) second_order_s=($number) algebra_s=($number) kkt_assembly_s=($number) factor_s=($number) solve_s=($number) update_s=($number) K_nnz=(\d+) rhs_cols=(\d+) Vxx_nnz=(\d+)"
$iterPattern = "FILTERDDP_ITER_TIMING iteration=(\d+) barrier_iteration=(\d+) backward_s=($number) forward_s=($number) total_s=($number) outcome=([a-z_]+) step_size=($number) backtracks=(\d+) mu=($number)"

$groups = @{}
$events = @{}
foreach ($line in Get-Content -LiteralPath $InputLog) {
    if ($line -match $stagePattern) {
        $key = "$($Matches[1]):$($Matches[2])"
        if (-not $groups.ContainsKey($key)) {
            $groups[$key] = [ordered]@{
                iteration=[int]$Matches[1]; barrier_iteration=[int]$Matches[2]; stages=0
                derivative_s=0.0; first_order_s=0.0; second_order_s=0.0; algebra_s=0.0
                kkt_assembly_s=0.0; factor_s=0.0; solve_s=0.0; update_s=0.0
            }
        }
        $g = $groups[$key]; $g.stages++
        $names = @('derivative_s','first_order_s','second_order_s','algebra_s','kkt_assembly_s','factor_s','solve_s','update_s')
        for ($i=0; $i -lt $names.Count; $i++) { $g[$names[$i]] += [double]::Parse($Matches[4+$i], $culture) }
    } elseif ($line -match $iterPattern) {
        $key = "$($Matches[1]):$($Matches[2])"
        $events[$key] = [ordered]@{
            backward_s=[double]::Parse($Matches[3],$culture); forward_s=[double]::Parse($Matches[4],$culture)
            measured_total_s=[double]::Parse($Matches[5],$culture); outcome=$Matches[6]
            step_size=[double]::Parse($Matches[7],$culture); backtracks=[int]$Matches[8]
            barrier_mu=[double]::Parse($Matches[9],$culture)
        }
    }
}

$rows = foreach ($key in $groups.Keys) {
    $g = $groups[$key]; $e = $events[$key]
    $componentSum = $g.derivative_s + $g.algebra_s + $g.kkt_assembly_s + $g.factor_s + $g.solve_s + $g.update_s
    [pscustomobject][ordered]@{
        iteration=$g.iteration; barrier_iteration=$g.barrier_iteration; outcome=$e.outcome
        stages=$g.stages; backward_s=$e.backward_s; forward_s=$e.forward_s
        total_s=$e.measured_total_s; derivative_s=$g.derivative_s
        first_order_s=$g.first_order_s; second_order_s=$g.second_order_s
        algebra_s=$g.algebra_s; kkt_assembly_s=$g.kkt_assembly_s
        factor_s=$g.factor_s; solve_s=$g.solve_s; update_s=$g.update_s
        unclassified_backward_s=$e.backward_s-$componentSum; step_size=$e.step_size
        backtracks=$e.backtracks; barrier_mu=$e.barrier_mu
    }
}
$rows | Sort-Object iteration,barrier_iteration | Export-Csv -NoTypeInformation -LiteralPath $OutputCsv
Write-Output "wrote $($rows.Count) timed backward passes to $OutputCsv"
