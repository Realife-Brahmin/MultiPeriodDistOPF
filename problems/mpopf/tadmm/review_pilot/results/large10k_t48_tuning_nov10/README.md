# Large10k T48 tuning campaign

Submitted 24 predeclared candidates in 12 jobs, two starting penalties per
job, with independent physical/timing audits. The active campaign is
`/home/aryan.jha/kamiak_jobs/tadmm_review_nov10_large10k_t48_tuning_v3`.
`submissions.tsv` records every configuration and job. The frozen source is
ad55d154. The native gate reproduced the historical default penalty sequence
and checked the new residual gate, decrease factor, floor, and Julia parsing.
The initial launcher attempts stopped before submission due to the symlinked
Julia path and soft scope in the gate. Their directories are retained.

Each job allocates all 64 physical cores on one Xeon 6438Y+ node, with 48
workers and 128 GiB requested memory. Each compares against a full MA57
centralized solve on its own allocation, with the same T48 input, objective
and initial state. tADMM uses full warmup and reset. The candidates vary
initial rho, decrease factor, floor, an optional all-copy RMS reduction
gate and exact versus 10/15/20/30-iteration local caps. All schedules in
this grid decrease or hold rho. Original stopping, strict feasibility,
and independent replay thresholds are retained. Save unsuccessful trials
and tuning costs. A 2x claim requires a validated candidate at half its
matched central solve time or less, then a matched confirmation run.

## Completed results

All 24 candidates and 12 independent audits completed. The per-candidate
results are in `completed-candidates.csv`, ordered by the ratio of matched
centralized time to tADMM core time. Full compact evidence is archived per
job. The complete replay tables and saved dispatches remain on Kamiak.

The fastest candidate used cap20, initial rho 15000, decrease factor 5 and
floor 24. It took 123 rounds, 125.85 s core and 140.70 s loop versus 87.85 s
centralized on cn168. Its objective gap was 3.93e-6 relative, maximum SOC
agreement error 0.3970 kWh or 0.02978 percentage points of rated capacity,
and maximum requested bound error 0.5419 kWh or 0.04064 percentage points.
The strict original-model maximum violation was 5.42e-4. The absolute cone
gap matched the central value of 1.55854e-5 pu squared to displayed precision.

None reached the physical near-optimal milestone or passed strict
chronological replay. None reached the combined objective and
original-model thresholds by half of matched centralized time. The
engineering storage thresholds were met at the endpoint by 22 of 24.
Peak sampled process-tree PSS ranged from 53.24 to 54.88 GiB.
The gated variants had a median core time of 261.65 s versus 157.32 s for
the ungated set. The grid also varied other settings, so this is a
descriptive comparison. No candidate qualifies for a speedup confirmation
or horizon extension under the declared criteria.

The fresh exact comparison 31679947 on cn168 took 129 rounds, 154.9708 s
core and 171.3019 s loop versus 87.7364 s centralized. Cap10 on the same
node took 176 rounds and 154.6995 s core versus 87.8136 s centralized.
The exact endpoint's maximum raw cone slack was 1.55854e-5 pu squared and
maximum SOC disagreement 0.03615 percentage points of rated capacity.
Its strict physical and chronological replay flags were false. These
single trajectories show that cap10 gave essentially the same core time
and higher raw cone slack than exact, rather than a validated speedup.

`exact-T48-31679947-evidence.tgz` preserves the exact timing traces,
source/input manifests, logs, memory and physical reports.
SHA256 `29d208d690813322422d678f1e33c30f9453db18e5300c203c845f198b700a4d`.
