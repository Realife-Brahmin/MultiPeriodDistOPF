# Nested battery-count campaign: final collected results

Completed and collected October22,2026. Forty-eight measured solves and six
T3 correctness smokes validated. This changes the battery population, not the
algorithm: a physical sensitivity study, NOT a pure computational ablation.
No production configuration or paper was changed on the basis of these results.

## Matched comparison

10k_X retains250 original battery sites (91414.8kW/365659.2kWh), covering102
root branches. 10k_Y retains100 (33330kW/133320kWh), covering100 branches,
and is nested inside X. No rating rescaling; network, loads, PV, retained
battery ratings and cost coefficients unchanged. Same selected sites at every
horizon. C_B2.4e-5, terminal soft penalty420. Each variant/horizon has its OWN
four mutually agreeing centralized references. Original1020 is historical
context, not a freshly matched third arm.

Same Gold6438Y+ node class,64CPUs,300GiB,Julia1.12.4,OMP/BLAS/Julia1.
Every measured case ran on cn164 except the T24 cases, which used cn167;
hardware records identify the same CPU model. DDP W=min(T,63), exact curvature,
tree/half-condensation, balanced scan, pass-GC, fixed production screening.
Full solve warm-up and identical physical initialization across methods.
Two repeats, order MA97/MA57/DDP then DDP/MA57/MA97 within each case.
No expensive overlap. Setup and validation are not optimizer times.

## Two-repeat median optimizer results

Full-convergence DDP versus fastest centralized backend (MA97 in all cases):

- X250,T24:2.3265 versus12.4515s;5.35x.
- X250,T96:6.8015 versus59.5335s;8.75x.
- X250,T144:10.5530 versus86.5430s;8.20x.
- X250,T288:18.2290 versus172.5550s;9.47x.
- Y100,T24:1.7045 versus12.2095s;7.16x.
- Y100,T96:5.0905 versus54.4625s;10.70x.
- Y100,T144:8.3005 versus81.9390s;9.87x.
- Y100,T288:16.8370 versus166.0100s;9.86x.

Shared-quality means the first observed point satisfying BOTH objective gap
<=0.001 relative to its own central reference and original-model primal<=1e-5.
No interpolation; fastest backend selected by its two-repeat quality median:

- X250,T24:DDP1.3527 versusMA976.8578s;5.07x.
- X250,T96:DDP3.4794 versusMA9736.0186s;10.35x.
- X250,T144:DDP4.7914 versusMA9750.9776s;10.64x.
- X250,T288:DDP11.3546 versusMA9760.8192s;5.36x.
- Y100,T24:DDP0.9110 versusMA974.4648s;4.90x.
- Y100,T96:DDP2.7948 versusMA9731.7231s;11.35x.
- Y100,T144:DDP3.9067 versusMA9723.6803s;6.06x.
- Y100,T288:DDP9.8446 versusMA9748.4010s;4.92x.

Thus DDP wins all tested cases on both clocks. Shared-quality speedup does
NOT always improve on full-convergence speedup, nor increase monotonically
with horizon. Feasibility and objective crossings occur at different iterates
in the two methods. Preserve these outcomes rather than extrapolating the
initial T24/T96 trend. Two repeats are descriptive, not strong statistical
evidence. DDP full-time repeat spreads are about0.2--11.4% of pair median;
the largest is the short Y100T24 solve. Per-run values remain available.

## Phase B interpretation

Median accumulated scan wall times for X250 at T24/96/144/288 are
0.426/1.037/1.988/3.172s (about15--19% of full solve); for Y100 they are
0.067/0.221/0.415/0.883s (about4--5%). This is mixed computation, memory
access, coordination and waiting, NOT communication-only time. At Y100T288,
Phase A is7.83s versus scan0.88s: after shrinking the battery state, Phase B
is no longer the dominant component in this variant.

Accepted iteration counts for X250 are12/15/15/14, while backward-pass counts
are20/24/24/22. Y100 counts are12/14/15/15 iterations and20/23/24/24 passes.
Barrier passes must not be mistaken for extra accepted iterations.
medians.csv/per_run.csv retain per-pass means/medians/min/max as well as totals.
For example X250T288's per-pass median is0.125s and max0.257s, not constant;
Y100T288 has median0.036s and max0.045s. These are medians over two repeats
of each per-run statistic. No causal communication decomposition follows.

The shared scan payload is8T(3nB^2+2nB) bytes: at T288, X250433152000bytes
and Y10069580800bytes. This excludes all private solver/model memory.
The Phase B improvement is motivation for structured value representations,
not proof of an exact low-rank model or safe physical battery aggregation.

## RAM and other clocks

Whole six-arm job cgroup high-water marks for X at T24/96/144/288 are
22.82/90.19/125.01/180.43GiB; Y22.54/88.92/127.58/180.93GiB.
All stayed below300GiB. These are allocation-wide peaks, not backend-specific
simultaneous peaks and not optimizer-only RAM. DDP coordinator/worker lifetime
RSS maxima and central process RSS are separately retained. Summing worker
maxima double-counts shared mappings and is NOT node peak RAM. Similar total
memory for X/Y shows that shrinking the scan payload does not eliminate the
dominant network/model/process storage.

per_run.csv also retains total process wall including build,warm-up,OpenDSS
and independent physical audit. These totals are deliberately not solve times.
Central build includes separately logged screening and power-flow initialization;
the first call can include compilation. Full internal Ipopt timing categories,
callback traces, printed original primal/dual/complementarity data and DDP
phase/direction/trial/GC traces remain in evidence.tgz. Do not relabel callback
mu as complementarity; that quantity is available in the internal raw table.

## Validation, provenance and reproduction

Source snapshot48099dbf, immutable location:
`/home/aryan.jha/kamiak_jobs/battery_count_oct21_campaign`.
Code branch later received history-only master merge9fd54956 and separate
representation note6e1f743c; neither changed the executed snapshot.
All8 case groups passed independently recomputed original-model residual<=1e-5,
objective agreement and OpenDSS. Both repeats have matching iteration/objective/
original-primal trajectories within the collector's explicit tolerance.
Initial objective agrees across all six arms in every group. Source and input
checksum checks passed at every case start. Hardware/affinity and private HSL
library hashes are archived; no licensed sources/binaries are included.

After extraction, the existing OpenDSS and DDP/central log validators were
rerun locally:108 checks across54 solves passed. Re-running the collector
from its repository location reproduced all five summary files byte-for-byte.

Preparation31345168 was preempted during setup after36seconds, automatically
requeued, then completed. sacct duplicate-attempt history is preserved.
The launcher reuses job-ID log paths; the first attempt's transient setup logs
were overwritten by Slurm/the restart, so do not claim all attempt stdout is
separately preserved. No interrupted warm-up is reported as a benchmark.
All measured jobs completed without preemption/failure.

Actual allocations were reduced from600 to300GiB while pending on user
instruction; frozen submit.sh still records600G. Hardware/cgroup logs establish
the actual300GiB allocation. No CPU,source,input or numerical policy changed.

evidence.tgz SHA256:
83affd96ee63893c5fc39056b1af874ef10762c67d95c2b50d7b3be6ad0d2928

The archive contains all raw text/CSV logs, manifests, construction/selection,
checks and completion records, excluding solution.jls binaries. All54 raw
dispatches remain at the immutable cluster path (about2.7GB total results);
solution.sha256 records their hashes. This is not a second off-cluster backup
of the dispatch binaries. Slurm attempt/accounting is in slurm_attempts.psv.

Extract evidence.tgz to an empty directory; run
`python ../collect.py EXTRACTED OUTPUT`. The collector checks completion,
hardware/source/input gates, physical residuals, paired trajectories, clock
alignment, own-variant references and finite monotone quality traces, then
regenerates per_run.csv,medians.csv,comparison.csv,job_memory.csv,validation.json.
DDP time is reconstructed from printed full solve minus summed per-iteration
times plus cumulative prior iterations and current backward endpoint. This
places residuals at the end of their corresponding backward pass and conservatively
assigns the small unaccounted timer remainder before the first point. Central
time is callback elapsed; original primal comes from the internal printed table
at its retained precision. These clocks are not interpolated exact crossings.
