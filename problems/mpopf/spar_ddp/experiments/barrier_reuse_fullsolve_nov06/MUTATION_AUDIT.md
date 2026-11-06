# Eligibility and actual mutation ownership

The immutable runtime is copied from dense_merge_teams96_nov05, the validated
current grouped-b2/right implementation. This overlay does not alter the
repository production driver or enable any default. The scoped run.sh clears
FILTERDDP/IPOPT settings then installs the known exact/tree configuration;
STUDY_BARRIER_REUSE survives that reset and is the sole baseline/cache toggle.

The private solver loop emits solve_start after reset, barrier only in the
actual converged-barrier branch, and accepted immediately before nominal
update. The gate consumes the event once, requires same epoch/k, exactly next
j, lower mu, identical configuration and complete worker ownership. Every
other/unknown event bumps the epoch. There is no positive permission inferred
from a repeated iteration number alone. The snapshot fixes model, source,
input, terminal penalty and path for each process, with before/after manifests.

Worker mutation sites in procs_pass:
- procs_worker_setup: deserialize data, set Cb, construct model/terminal gamma,
  initialize solver and shared slots. All are before any population.
- pw_step1(elem): install initial nominal only on new solve, which clears cache.
- pw_pass_start(advance): swap nominal/current only after accepted event. Fresh
  population rewrites every cached matrix even when storage buffers persist.
- pw_rollout: write trial current, leaving nominal x/u/phi/z/slacks untouched.
  A rejected trial cannot grant barrier permission. Accepted swap invalidates.
- pw_nominal at solve exit: return/swap final accepted state; no subsequent pass
  within that solve. New solve clears validity and payloads.
- _tp_pre: fresh derivative/Hessian/network/KKT preparation and RHS at every
  pass, including hits. Global KKT pattern/layout caches are scratch, not saved
  Phase B factors. Its supported route has no diagonal shift/regularized retry;
  a failure returns nothing/stops. Applied regularization is established from
  this route, never from the downstream reset of d.reg_last.
- _tp_post: write current policy/update objects from current S and suffix V/v.
  They are deliberately not reused. Pass-GC drops S/policies, keeping only the
  separate minimal Phase B cache. Every eventual closed V is restored before
  current-stage consumers read it; hit intermediate matrices are unmaterialized.

No runtime mutation callback writes model data, bounds or terminal coefficients
in this pinned full-solve driver. In-place external mutation of those objects
is outside its supported contract. barrier_cache_invalidate!(reason) invalidates
coordinator eligibility and, if called locally on a worker, its own validity.
During a hit, a late local invalidation fails explicitly. A future API that
mutates a distributed model must call invalidation on every affected owner and
coordinator BEFORE the pass; omitting that is unsupported. Debug mode hashes
actual nominal state, duals, bounds, matrix-driving derivatives, factors and
terminal stage to audit the ownership proof. It also deliberately mutates seven
categories of actual stage data and rejects their fingerprints. These expensive
checks are absent from clean clocks. This is a scoped exact driver optimization,
not a claim that arbitrary unannounced mutation can be detected in constant time.

The complete source-level barrier invariance argument and global original-
balanced capture verification remain in phase_b_barrier_cache_nov05. Here the
right-sided grouped kernel is separately checked at T7/T24/T96. Full runtime
validation compares against a fresh shadow using the same BLAS thread count,
raw UPPER convention and actual full horizon; matrix digests require exactness.
Current policies/actions/directions retain the unchanged scientific error gates.

Diagnostic-only repair after live budget diagnosis: fresh shadow computation
moves to an existing worker already pinned to two physical cores. All Phase A
owners are synchronized first; the reference uses private leaf copies and leaves
shared operands untouched. The coordinator retains its single-core affinity.
Clean baseline/cache scheduling and kernel remain identical. Both warm and
measured exits now execute the same retirement method, and debug warmup logs
are retained separately. Negative-control flags reset at new solve so validation
evidence is present in the retained measured debug trace as well as warmup.
