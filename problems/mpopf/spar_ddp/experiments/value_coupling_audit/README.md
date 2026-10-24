# Value-function coupling and compact-computation investigation

**Status: bounded campaign complete; revise consumers before integration.**
Read [the final decision](results_campaign_oct23/FINDINGS.md) first.
This directory preserves useful research code and evidence, not a production
compressed solver. No script here runs automatically from the normal driver.
All physical batteries and constraints were retained in these tests.

## Read in this order

1. [Final all-four-stage campaign](results_campaign_oct23/FINDINGS.md):
   matching elimination records, direct initial generation, compact operators,
   all-stage directions and honest fresh-input costs. Accuracy passed;
   the total-cost gate failed. No optimizer integration was justified.
2. [Earlier direct-generation coverage](results_direct_initial_oct23/FINDINGS.md):
   only stage1 records were then available. The final campaign closes this gap.
3. [Compact chaining](results_chain_oct23/FINDINGS.md) and
   [component kernel](results_compact_oct23/FINDINGS.md): component progress,
   not interchangeable whole-workload timings.
4. [Block-scaled compression](results_oct23/FINDINGS.md) and
   [frozen audit](results_oct22/FINDINGS.md): error amplification, reference-path
   diagnosis, rank growth and adverse medium-late findings.

## Source and evidence map

- `campaign_*`: final capture, complete direct-input/consumption workload,
  synthetic checks and immutable-job harnesses.
- `direct_initial*`: factors generated from existing stage/network elimination.
- `compact_kernel*`, `compact_chain*`: compact operators and retained-interior
  compositions, including raw asymmetry/original-triangle handling.
- `followup*`, `hierarchy*`, `replay*`: earlier representation and error studies.
- `analyze_*.py`: strict posthoc collectors; raw evidence is under dated
  `results_*/` directories, including failed/interrupted attempts.
- `CAMPAIGN.md`, `DIRECT_INITIAL.md`, `COMPACT_CHAIN.md`, `COMPACT.md`,
  `FOLLOWUP.md`, `PLAN.md`, `REPLAY.md`: dated experimental specifications.
  Their future-tense tasks and absolute cluster paths are historical context,
  not live queues or portable user entry points.

Binary Newton captures remain at the immutable locations recorded in their
manifests, not in a fresh clone. Some harnesses deliberately depend on pinned
older snapshots. Reproduction requires those inputs and the recorded commit;
do not run an old launcher against today's working tree without review.
Coefficient storage, object size and resident memory are different measures;
component timings do not establish parallel Phase B or optimizer performance.
