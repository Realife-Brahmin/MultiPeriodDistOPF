# Nested large10k battery-count study (authorized October21)

10k_X retains250 original battery sites;10k_Y retains100, nested in X.
Original network1020 batteries is context, NOT a fresh matched timing arm.
T24/96/144/288 x2 variants x MA57/MA97/DDP x2 repeats =48 measured solves.
No objective/cost retuning or total-capacity rescaling. C_B=2.4e-5 gamma420.
PV and every non-battery dictionary field are unchanged. Internal system ID
stays large10kC_1ph in isolated variant snapshots so coefficient lookup is safe.

Construction orders original batteries by feeder-root branch (integer bus ID),
then bus ID within branch, interleaving one rank across all sorted branches.
Y is first100, X first250, same sites all horizons. Reports all1020 ranked sites,
feeder coverage and original retained kW/kWh. Each battery parameter dictionary
is filtered consistently. Matching private Storage.dss keeps exact original
declarations; PVSystem.dss must be byte-identical. No shared input/source edits.
All non-battery fields and retained parameter values are asserted equal.

Original periodic T24/96/144/288 exports are copied/hashed from frozen Oct16
inputs. Only correctness T3 is generated from unchanged original DSS first.
Preparation runs T3 on X and Y with all3 methods, original-model residual
audit, same objective and independent OpenDSS before releasing8 case jobs.
Each case allocation runs both repeats, reverses method order on repeat2.
Dependencies afterok serialize all allocations and enforce smoke gates.
Matching Gold6438Y+64CPU600GiB/full-node occupancy check, Julia1.12.4,
OMP/BLAS/Julia1; DDP W=min(T,63), full warmup every optimizer.
Four-hour allocation cap per job,50min cap per individual solve (including
build and warmup). Failures/preemptions preserved, no automatic scope expansion.

Current central source had lost Oct16 uncommitted benchmark warmup/matched-start
controls. Restore opt-in full IPOPT_WARMUP and IPOPT_MATCH_DDP_START, matching
the frozen Oct16 implementation/options:tol1e-6,monotone,original inf_pr,
linear_system_scaling=none,MA97AMD or MA57pivot_order2. One actual numerical
case per method warmed then restarted at original physical values, checked by
initial objective. No nested multithreading or private HSL binaries in Git.
Shared initial ell incorporates the same positive0.01 cone slack as DDP.
Existing DDP process original-feasibility audit enabled, objective added to its
log for exact pairing with phase/timing records. No expensive per-merge tracing.

Log full optimizer, setup/build, initialization(first-call/JIT included), warmup,
validation, total process wall, Ipopt timing statistics, iterations/objective/
original primal/dual/complementarity and existing DDP scan/pre/post/directions/
trials/GC. Parse per-pass CSVs after each solve; print no new per-merge events.
Record shared buffer 3nB^2+2nB Float64 entries per time-step and exact all-period
bytes. Worker RSS lifetime maxima include warmup/shared mappings; their sum is
NOT simultaneous RAM. Batch cgroup/Slurm high water is separate metadata.

Posthoc reporting target: objective gap<=0.001 AND original primal<=1e-5,
using each variant/horizon's own mutually agreeing central references.
Use2repeat medians and faster central backend median, never fastest individual.
Ipopt original-primal printed precision retained. No main-paper edits until
collected evidence. The central callback does not report complementarity;
retain full internal logs/table/final complementarity without inventing a
callback field. Analyze full/quality clocks separately. All raw dispatches saved.
If afterok blocks, inspect parentfailure; never bypass numerical/smoke gates.
