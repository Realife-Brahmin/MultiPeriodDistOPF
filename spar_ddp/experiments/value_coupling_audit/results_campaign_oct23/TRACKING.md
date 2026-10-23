# Autonomous compact campaign tracking

Authorized scope/gates: ../CAMPAIGN.md (takes precedence over historical narrow
no-capture limits). Dedicated heartbeat spar-ddp-compact-campaign-to-decision.
Base52d1183c; latest branch must be verified on every continuation.

## Input coverage

Original capture records only1/12/24. New minimal private runtime copies1.1MB
driver code, links immutable package/env/raw input/validators; late-only packets
1:4 and all24leaves are expected to add about1GiB. Home initially5.3GiB free,
disk guard>=2GiB before starting; no existing-work deletion.

Initial capture31369406/70442034 failed BEFORE any Julia simulation: inherited
manifest used paths relative to the ORIGINAL root but was checked in the new
sparse runtime root. Thus59 intentionally uncopied files appeared missing and
the deliberately replaced helper differed. Original root independently rechecked
and EVERY inherited checksum passed. Fix is cd to original root for INHERITED
check; new runtime/overlay/input checks remain separate and unchanged. No hash
filtered, no scientific gate loosened. Old snapshot preserved immutable at
/home/aryan.jha/kamiak_jobs/compact_campaign_capture_oct23.
capture_setup_failed.tgz SHA256
69ee01a5e70568cd21da0b6fa8e4163b63d37317203953d5d706b301d7eb5c6b
contains original failure stdout/stderr/check logs/overlay/manifests/accounting;
the1.1MB runtime stays on cluster pinned by runtime.sha256. No trajectory data
was produced. Safe replacement must be a NEW snapshot.

Next: replacement smallT3 gate then matching largeT24 capture; inspect full
logs and numerical/physical/OpenDSS/capture-equivalence checks. Once verified,
apply direct generator to ALL1:4 and implement/measure all four-stage suffix
and policy consumers with honest boundary cost. Do not stop at capture.

Replacement31369418/a2e16a52 completed the T3 dense solve, numerical validation,
independent physical original_max2.10e-16 and OpenDSS(max voltage1.59e-9pu),
then FAILED in posthoc checker because its top-level probe accumulator used
Julia soft scope. LargeT24 never started. All data/logs remain in
/home/aryan.jha/kamiak_jobs/compact_campaign_capture_oct23_v2.
capture_checker_failed.tgz SHA256
aa4a4d69a232a0062af35ff1194456b9458f5f75c039cc0b31e57a3faf8edb54
preserves full failed attempt excluding reused runtime/depot. Checker is now
wrapped in a function; NO numerical change. Next snapshot REUSES this already
completed smoke data read-only, revalidates with corrected checker into NEW
report location, and ONLY then runs largeT24. Old failed job stays failed;
no fake COMPLETE or repeated smoke simulation. Source/runtime/log hashes pin
the reused evidence. All four-stage compact value/operator core and synthetic
unit have also been implemented/passed; full fresh-cost harness remains next.
