# Revision experiment criteria — November 9

Use current corrected OPF models first. Preserve every run, failure and source
manifest. Centralized remains fully converged with its existing tolerances;
only tADMM local subproblems may use persistent iteration caps. Timing includes
coordination/communication; physical and OpenDSS audits remain outside it.

## Required evidence

For every centralized and tADMM candidate, explicitly evaluate branch-by-period
SOCP slack `v_i*ell_ij - P_ij^2 - Q_ij^2` (v is squared voltage), not merely
inequality violation. Retain signed slack, absolute maximum/p95 and relative
maximum/p95 with denominator and zero-flow floor stated. The current exploratory
absolute target is1e-5 per-unit squared; report magnitudes and physical replay,
not a blanket exactness theorem. Small absolute gaps near zero flow can have
large relative errors. Existing physical metrics also retain implied excess
line-loss energy; report this alongside objective and constraint agreement.

`audit-socp-storage-nov09.py` independently recomputes these gaps from exported
branch/voltage tables after checking table SHA256. Per-branch records are saved
as `socp-tightness.csv.gz`. The audit is passive and cannot change a solution.

Storage must execute chronologically in OpenDSS with100% efficiencies, zero
idle/leakage losses, state initialized once and no resets/clamping used to
repair the schedule. Report OPF/OpenDSS SOC differences in percentage points
of each battery's rated capacity, energy-accounting and bound errors as
percent of that capacity, and dispatch error as percent of rated power.
Rated-power normalization avoids division by nearly zero commanded power.
The initial engineering target is0.05 for these percentages;0.5 is concerning.
Keep raw kW/kWh and unchanged strict feasibility flags separately: an engineering
agreement target does not redefine strict mathematical convergence.

Continue recording per-process/node RSS, PSS, private memory and CPU, with
Slurm memory/accounting and sampled-peak limitations. Keep centralized and
tADMM phases distinct using recorded phase events. A centralized failure or
memory exhaustion is evidence only when logs, allocation limits and solver
status establish its cause; do not assume an HPC failure from PC experience.

## Tuning and horizon progression

Measure local-solve median/p95/max and per-iteration slowest-versus-median
time, plus truncation counts, outer iterations, end-to-end time and final
physical errors. Select caps empirically per system/horizon, with equal
models/inputs and a fully converged centralized comparator. Cap10/30 are
initial trials, not a universal cap heuristic. Do not declare a speedup from
iteration counts or objective agreement alone.

Review fresh med2522 T96 cap10/30 outcomes before expanding. If they are
competitive and physically sound, extend large10k T24 toT48 thenT96, preserving
one worker per time step; extend med2522 beyondT96 only with a declared
hardware allocation capable of the worker count. Hold/retune when physical
errors, SOCP slack or runtime is unacceptable. First continue trajectories
to stricter thresholds and record milestone times where justified; largerT
does not itself guarantee better convergence. Use matched repeats for paper
timing claims and retain unsuccessful configurations in the tuning record.

## Advertised July configuration

The immutable advertised tag `260731-ias-trans-tadmm` resolves locally to
commit50b8901c324ba9eefd0d672c6e484e0ededc7a57. Its test inputs are under
`rawData/`, unlike current `test_networks/`. Its release annotation explicitly
calls out the large10k T48 memory baseline. That historical observation is not
an HPC memory result.

If needed, create a separately named release-input comparison using those
network/resource/profile inputs. Retain corrected line-loss equations and
the same explicit SOCP/chronological OpenDSS validation. Pin all input hashes,
document differences in resource sizes/load levels/profiles and apply the
same configuration to centralized and tADMM. Do not silently replace the
current-input benchmark or restore a known modeling error. No release-input
optimization campaign has yet been launched.

## First saved med2522 replay

Gate31667879 passed native tests and replayed the saved T96 MA57 case31660059.
Centralized passes chronological execution. The old tADMM endpoint had hit
1000 iterations, did not meet the final RMS criterion, and had38.97% objective
disagreement before replay. Its96.61055kWh requested SOC-bound violation and
up to6.385kW dispatch refusal are failed-run diagnostics, not converged results.
Independent maximum absolute cone gaps are7.79199e-7 (centralized) and
8.19689e-7 (failed endpoint); p95 relative gaps are0.23945% and0.38098%, with
maximum relative gaps near unity at tiny flows. Full chronological and
original-model validation are indispensable even when absolute cone slack
is small. Fresh single-node cap10/30 results remain pending.
