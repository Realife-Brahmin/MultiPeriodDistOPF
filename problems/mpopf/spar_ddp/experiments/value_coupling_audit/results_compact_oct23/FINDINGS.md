# Compact Schur component: bounded assessment

Implemented and verified a genuinely compact factor/operator component. With
input factors available, its merge/factor-update, condition-probe and three
Schur-action costs total a median0.02595s versus0.16238s for dense merge/actions
(about6.26x faster). However, fresh factor discovery plus the dense consumer
boundary remove that benefit: the charged hybrid pipeline is2.89452s versus
1.29224s dense (2.24x slower). No medium or optimization tests were triggered.
This supports compact-kernel engineering, not a claim that compression failed
or that a compressed SPAR-DDP optimizer is already faster.

## Implemented operation and accuracy

One real periods3/4 leaf merge in the SAME original1020 late frozen four-period
window, true dense future from the24-period capture. Sourcecede3fc6,
job31368876 COMPLETED0:0 in1m18; Gold6438Y+,64CPU300GiB,Julia1.12.4/BLAS1/OMP1,
occupancy guard. All source/input/capture hashes, synthetic unit tests and
saved dense policies pass; saved probe discrepancy0. Strict collection passes
three alternating identical-input warmed pairs in one process.

The kernel uses hierarchical diagonal/local blocks plus off-diagonal factors,
recursive Woodbury elimination and reused child factorizations. Factor updates
use QR and small-core SVD, rather than SVD of newly reconstructed full outputs.
All physical batteries/constraints, full vectors, strong local blocks and raw
input asymmetry remain. Original-coordinate upper completion explicitly matches
production's Symmetric(K,:U)/copytri! convention, even for interleaved clusters;
it is not arbitrary symmetrization of Mab or raw input storage.

Merged Maa/Mab/Mbb/ma/mb and hybrid suffix values have maximum scaled discrepancy
1.031e-12; Schur probe actions7.633e-13 and compact-solve residual1.451e-16.
Recovered battery paths and selected stage1 full control/equality/bound directions
pass the unchanged1e-7 scaled and1e-5 relative-with1e-10-floor gates. Maximum
direction discrepancy6.835e-21 is essentially the dense-window association
baseline, NOT proof of bit-identical compact arithmetic. A deliberately corrupted
merged gradient produces selected full-direction discrepancy15.01 and is rejected:
this proves the saved-policy audit actually depends on the tested merge.

No full nB-by-nB matrix is assembled inside the compact merge/actions. Measured
local dense factor size<=60 (cap64), Woodbury core size<=4, QR/SVD core size<=30;
563 individual cross-block dense fallbacks are charged. Maximum retained rank3
refers to updated/recompressed interfaces in this component, NOT all original
input factors or a temporal/horizon rank bound. General rank growth is not
resolved, and a core reaching its parent dimension would explicitly reject the
compact path. The condition probe is an uncertified heuristic, not a condition
bound or feasibility certificate.

## What is not yet compact

Outputs are factored Schur operators, not a closed HODLR representation ready
for the next arbitrary scan combination. Deferred actions do real compact
computation, but today's downstream scan/policy code requires dense blocks.
An explicit export applies operators in batches of16 columns and constructs
all three dense matrices OUTSIDE the compact kernel. Its full cost is charged;
it is not hidden as validation-only free work. A hybrid four-period replay
replaces ONLY this merge, with remaining operations unchanged and dense.
This is neither a fully compact scan nor PhaseB parallel wall timing.

## Warmed cost accounting and reuse assumptions

Medians over3paired repeats: dense merge0.16000s, actions0.00236s. Compact
factor/update0.01927s, condition probes0.00139s, actions0.00530s. Input cross-block
factor discovery (six raw blocks) costs0.83599s each fresh use, dense export
0.88723s, and unchanged downstream replay0.74303s (dense0.73940s). Explicit
cleanup GC costs0.39439s compact versus0.39083s dense. Component medians need
not sum to median phase totals. All phases record allocation and within-phase GC.

Total measured phase sums: compact2.86961--2.89489s, median2.89452s;
dense1.28980--1.29556s, median1.29224s. These include matching probe consumption
and explicit cleanup, unlike the prior0.89s four-period dense-only clock.
Without explicit cleanup, execution medians are2.49077s versus0.90205s; the
negative pipeline conclusion does not depend on common cleanup overhead.
Even with already-cached inputs (dropping fresh discovery from the calculation),
the dense-export pipeline median2.05630s exceeds dense1.29224s. Do not use
unbounded cache reuse to claim an end-to-end win.

Caches are valid only while this frozen iterate's matrices are unchanged.
Initial dense cross-block SVD discovery is INCLUDED in fresh-input totals;
it has not been replaced by a generator using original stage tree-KKT factors.
The existing topology and boundary/interface factors are reused within the
component; this is not evidence of reusing numerical factors across Newton
iterations. Cold initial topology setup0.86754s includes compilation and is
reported separately, NOT mixed into a warmed benchmark. Fixed topology may
be reused for a fixed network/battery set, but that does not make fresh numerical
factor discovery free. Raw phase data and individual repeats are preserved.

## Storage and memory limits

Two compact input elements store1,321,104 coefficient bytes versus49,971,840
for their six raw blocks/vectors (about37.83x fewer). This denominator differs
from the previous29x mean across a whole window's represented elements. The
compact operator/factor object graph occupies3,785,326 bytes by summarysize;
tree/index/object overhead is included there, not in coefficient bytes.
Neither number includes a dense legacy export or proves resident-RAM savings.

Allocation is cumulative traffic, not resident memory: export alone allocates
about2.687GB per repeat; new input factor discovery about0.415GB. Full reference
and compact data coexist in this diagnostic. Whole-job RSS includes loading,
validation, both paths and dense exports, so no isolated candidate/dense peak
resident-memory comparison is available (Slurm sampled batch MaxRSS2,560,292KiB).
No memory saving for an optimizer,
whole PhaseB or full horizon is claimed.

## Decision and specific next obstacle

Large-window value/action/direction accuracy passes for this component; total
cost does NOT justify the user-gated extension to medium or full optimization.
Published benchmark rows, production defaults and paper remain untouched.
No compressed optimizer has been tested: original-model feasibility/OpenDSS,
objective quality, iteration counts and near-optimal/full convergence times
remain unmeasured, not failed.

The specific next obstacle is CLOSED compact-to-compact Schur composition and
compatible factor-aware consumers, avoiding dense export and repeated initial
factor discovery. A future kernel must propagate/recompress updated boundary
factors, preserve original-coordinate triangle semantics and full-stage
direction sensitivity, and charge interface-rank growth plus fallback costs.
Direct generation from stage network separator factors is still unimplemented;
small per-period interfaces alone do not prove compact temporal closure.
This report does not automatically launch that larger implementation/campaign.

## Evidence and excluded pilot

v2_evidence.tgz contains all raw logs/CSV/source/hardware/accounting/manifests,
SHA25680f6d9f9b3d16dbb559a8c16433ab29e8dedbe7ec3dbec9dd5efceda328e309f.
Reproduce with analyze_compact.py EXTRACTED_ROOT OUTPUT. analysis/validation.json,
phases.csv,totals.csv and evidence.sha256 preserve strict completion and costs.
Initial pilot31368859 is preserved separately in initial_pilot_evidence.tgz:
its periods1/2 merge did not feed selected policy probes and its first dense
action had compilation. It is NOT used as the validated warmed comparison.
