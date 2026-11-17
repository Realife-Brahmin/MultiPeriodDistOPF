# DDP scaling effort, November 16-17, 2026

Question (user): DDP's advantage over Ipopt should grow with the horizon, but on
large10k the ratio falls after `T=96` (October campaign, Ipopt / DDP at the 0.1%
rule: 2.47 at `T=96`, 1.67 at `T=144`, 1.29 at `T=288`). Find out why and fix it
with methods that help every system, and add the PESGM systems med3980 and
sf6937 as test cases. Branch `ddp-scaling-effort-nov16`. Exploratory: no
campaign row, solver default or paper is changed.

## Findings

1. **Why the ratio falls: the battery sweep runs out of workers.** The pairwise
   rounds do about (T/2) log2(T) full combinations a pass. On 63 workers they are
   7 combinations deep at `T=96` but 15 at `T=144` and 25 at `T=288`, each a dense
   n_B x n_B combination (n_B = 1020 on large10k).
2. **Chunked scan (`FILTERDDP_PROCS_CHUNKED=1`, commit 2c0fa98c), general.**
   Each worker reduces its run of consecutive periods, the rounds run over the 63
   summaries, each worker closes its periods: about 2T/P + log2(P) deep. Used only
   where a depth model says it is shallower (not at `T <= 96` on 63 workers).
   Same iterations and objective to every printed digit wherever compared.
   - Kamiak, large10k `T=144`, 63 workers, both arms in job 31842826 on one node
     (cn150): solve 132.4 s (rounds) -> **76.1 s** (chunked), 18 iterations each,
     objective 2707473.330909072887 each, OpenDSS clean (0 of 1,486,080
     bus-periods). The battery sweep 3.88 -> 1.88 s a pass (depth model 15 -> 9.2).
     For scale only (other node, other code): the October campaign's Ipopt took
     127 s (MA97) and 145 s (MA57) to full convergence at `T=144`.
   - Same job, large10k `T=288`: 237.0 s (rounds) -> **116.4 s** (chunked), 18
     iterations each, objective 2710205.787244017702 each, OpenDSS clean (0 of
     2,972,160 bus-periods). The battery sweep 7.27 -> 2.78 s a pass (depth model
     25 -> 12). The October campaign's Ipopt: 243 s (MA97), 381 s (MA57).
   - 309 PC, med2522 `T=24` on 4 workers (one feeder): battery sweep 103 -> 58 ms
     a pass, solve 32.7 -> 29.4 s, 52 iterations, same objective.
3. **Step 1 allocations (commit 83601e98), general, bit-identical.** Dense zero
   second derivatives, copies of sparse matrices, a rebuilt fixed-pattern
   Hessian: 0.54 -> 0.40 GiB a pass on large10k `T=6`, no full collection.
4. **Feeder blocks (`FILTERDDP_VALUE_BLOCKS`, commit dd58034b): large10k only.**
   The battery value function splits exactly by feeder where the feeders hang
   off the slack bus. large10k has 102 such feeders (193 -> 55 s sequential at
   `T=24`); ieee123, med2522 have one, and med3980 and sf6937 feed every
   substation through one line (`feeder_blocks.jl`, `feeder_blocks_dss.jl`). The
   user rejected topology-specific speedups (2026-10-06, again 2026-11-17): kept
   off by default as an observation and a lead for spatial decomposition.
5. **med3980 and sf6937 with one substation (commits bb23b311, 8cd23733,
   9f87b11e).** The shared parser reads named buses; fixed capacitors are modeled
   exactly (q = qC v, linear in the branch-flow voltage variable); decks
   `test_networks/ieee9500xf_sub1_1ph` and `smartds_p1u_xf_sub1_1ph`:
   substation 1 feeds the network radially through one tie per other area
   (`single_source_scan.jl`). Without capacitors neither solves at full load;
   with them med3980 does (load scale 1), sf6937 needs 0.6 (loads, capacitors,
   PV and batteries scaled together, `scale_deck.py`, `lambda_search.sh`).
   First runs (`new_systems_pipeline.sh`, 309 PC): Ipopt `T=24` solves both (103
   and 80 iterations), DDP `T=6` converges on both (99 and 88 iterations, within
   1.1e-5 and 1.3e-6 of Ipopt), every dispatch clean in OpenDSS (6e-6 pu or
   better). DDP needs many more iterations here than on large10k (13-18).

## med3980 against Ipopt (Kamiak job 31842910, one 64-core node)

Campaign protocol (`kamiak/new_systems_job.sh`): DDP on 63 worker processes,
Ipopt MA57/MA97 on one thread, shared load-flow start, full warm-up, every
dispatch clean in OpenDSS. Time to the 0.1% gap with primal infeasibility at
most 1e-5 (internal residuals for both, `near_opt_posthoc.py`,
`ipopt_trace_posthoc.py`), full solve in brackets:

| T | DDP | Ipopt MA57 | Ipopt MA97 | best Ipopt / DDP |
|---|---|---|---|---|
| 24 | 10.9 s (11.3) | 41.6 s (42.8) | 37.9 s (39.1) | 3.5 |
| 96 | 33.1 s (34.5) | 187.1 s (196.2) | 165.9 s (175.1) | 5.0 |
| 144 | 41.4 s (41.7) | | 235.3 s (248.5) | 5.7 |
| 288 | 81.8 s (82.3) | | 487.0 s (530.5) | 6.0 |

105 and 132 DDP iterations against 110 and 112 for Ipopt; DDP's final objective
within 1.8e-5 and 4.0e-5 of Ipopt's. On this single-feeder system the ratio
grows with the horizon. Logs in `kamiak/results_31842910_med3980/`.

## sf6937 against Ipopt (Kamiak job 31842911, same protocol)

Load scale 0.6 (one substation). Time to the 0.1% rule, full solve in brackets:

| T | DDP | Ipopt MA57 | Ipopt MA97 | best Ipopt / DDP |
|---|---|---|---|---|
| 24 | 33.5 s (36.9) | 66.2 s (68.3) | 62.4 s (64.8) | 1.9 |
| 96 | 134.3 s (141.7) | 365.4 s (377.6) | 333.2 s (346.2) | 2.5 |
| 144 | 165.9 s (176.8) | | 530.0 s (563.8) | 3.2 |
| 288 | 323.3 s (323.3) | | 1314.4 s (1439.1) | 4.1 |

118 and 184 DDP iterations against 102 and 127 for Ipopt; DDP's final objective
within 4.6e-6 and 1.6e-5 of Ipopt's; every dispatch clean in OpenDSS. At `T=96`
a pass takes 637 ms: battery sweep 341 ms (rounds, 7 deep on 63 workers; 506
batteries), step 1 212 ms, step 3 63 ms. The sweep and the iteration count are
the levers left here. Logs in `kamiak/results_31842911_sf6937/`.

## Separate primal and dual step lengths (`FILTERDDP_SPLIT_STEP=1`)

Most of sf6937's DDP steps were cut short by the fraction-to-the-boundary rule
(106 of 193 steps below 0.1 at `T=96`, no backtracks). The solver's existing
option gives the bound multipliers their own step, as Ipopt does. Iterations at
`T=24` (309 PC, `run_paired.sh count`, same objective to 5e-7 in every pair):

| | ieee123 | med2522 | med3980 | sf6937 | large10k |
|---|---|---|---|---|---|
| one step | 22 | 52 | 105 | 118 | 14 |
| separate steps | 25 | 48 | 86 | 83 | 14 |

At `T=48`: 24 -> 25, 60 -> 52, 115 -> 94, 153 -> 100, 14 -> 14 (sf6937's gain
grows with the horizon). Objectives agree to 2e-6, med3980 `T=48` to 1.3e-4.
Not yet in our solver: Ipopt's safeguard that pulls a bound multiplier back
towards mu / slack when it drifts too far (Waechter and Biegler 2006, the
source of the separate multiplier step).

(October's test, 60 against 31 on ieee123, was under the diagonal Hessian.)
Used for the two-repeat runs (`kamiak/campaign_job.sh`, `DDP_SPLIT=1`).
Stopping early and repairing (tADMM's step) was considered and not pursued:
DDP's objective enters the 0.1% band only a few iterations before its residual
reaches 1e-5 (sf6937 `T=96`: 159 against 175 of 184; med3980: 122 against 127).

## Files

- `kamiak/job.sh`: large10k rounds vs chunked (jobs 31842814, failed at launch on
  CRLF line endings; 31842826). `kamiak/results_31842826/`: summary, hardware and
  logs (gzipped). `kamiak/new_systems_job.sh`: med3980 and sf6937 against Ipopt.
- `local_309/`: logs of the 309 PC runs above (gzipped).
- `lambda_search_logs.tgz`, `first_runs_logs.tgz`: the load-scale search and the
  first runs on the new systems.
- `check_jacobian.jl`: stage Jacobian against central differences.

`T=144` and `T=288` (jobs 31843080, 31843081; MA97 only, the faster Ipopt at
`T<=96`) use the chunked scan (depth 15 -> 9.2 and 25 -> 12). DDP iterations:
med3980 129 and 155, sf6937 195 and 222 (Ipopt 112, 114, 134, 162). med3980
`T=144` ends 4.5e-4 above Ipopt's objective, inside the 0.1% rule but larger
than elsewhere (1e-5 to 1e-4). Logs in `kamiak/results_31843080_med3980/` and
`kamiak/results_31843081_sf6937/`.

## Not done

Two repeats per case for any reportable number; the paper.
