# Changes made during the FilterDDP experiment

## Changes to the authors' repository: one patch, added 2026-08-12

`ddp/external/FilterDDP.jl` is at upstream commit `513a104` plus exactly one
patch, saved as [`ddp/patches/per_stage_data.patch`](../patches/per_stage_data.patch)
(108 lines, touching `src/ocp/ocp.jl`, `src/backward_pass.jl`,
`src/forward_pass.jl` and nothing else).

**What it does.** `OCP` held one `stage_objective` and one `constraints` for the
whole horizon, so per-period data could not be indexed by `t`. The patch adds two
optional fields, `stage_objectives` and `stage_constraints`, each a vector with
one entry per stage, and two accessors `stage_obj(ocp, t)` / `stage_con(ocp, t)`
that the two passes now call. Both fields default to `nothing`, in which case
every call site falls back to the single shared function.

**Why it is justified.** This is not an extension of the algorithm. Remark 1 of
the global-convergence paper states plainly that the objective, dynamics and
constraint functions *"can in general be time-varying but we avoid specifying
this explicitly for notational simplicity."* The capability was in the method and
merely absent from the implementation, because every example the authors ship is
time-invariant (grep for a time index across `experiments/` returns nothing).

**Regression evidence that the default path is unchanged.** The authors' own
`experiments/filterddp/double_integrator.jl` was run on the patched clone:

| | iterations | objective | primal infeasibility |
| --- | --- | --- | --- |
| Stage 3 record, unpatched | 51 | `1.26574863e+00` | `8.09e-08` |
| same script, patched | 51 | `1.26574863e+00` | `8.0917e-08` |

Identical. Their `results/` directory was restored with `git checkout` afterwards,
so the only difference from upstream is the three patched source files.

**What it bought.** The copper-plate model no longer needs the Lagrange
interpolant or the time-index state: `nx` drops from 2 to 1, and the horizon
ceiling disappears. With the interpolant, `T = 48` died with a
`StackOverflowError` inside Symbolics while still *building* the constraint
function. With per-stage data, `T = 192` builds and solves in 12 s to
`|ΔJ| = 7.1e-15` against the closed form.

## Previous status (Stages 1-10): no changes to the authors' repository


`ddp/external/FilterDDP.jl` is at upstream commit `513a104` and
`git status --porcelain` on it is **empty**. No source file, example or result file was
modified. Concretely:

- **No code fixes were needed.** The repository installed, precompiled and ran its example
  correctly on the first attempt under Julia 1.12.6. Nothing was patched, no warnings were
  suppressed and no tolerance was loosened to obtain a passing run.
- `Pkg.instantiate()` created `Manifest.toml` files inside the clone. These are
  **gitignored upstream** (`.gitignore:24`), so they do not dirty the tree.
- `Pkg.instantiate()` rewrote `experiments/filterddp/Project.toml` with CRLF line endings
  on Windows. `git diff` showed **no content change**, only line endings. Restored with
  `git checkout -- experiments/filterddp/Project.toml`.
- Running `experiments/filterddp/double_integrator.jl` unmodified overwrites
  `experiments/filterddp/results/double_integrator.txt` inside the authors' tree, because
  the script writes there by design. The shipped file was copied to
  `results/official_example/double_integrator_AUTHORS_shipped.txt` **before** the run, our
  output was copied to `..._OUR_run.txt`, and the original was then restored. Both are kept
  so the reproduction can be checked (they agree — see README).

## Environment decisions

| Decision | Reason |
| --- | --- |
| **Julia 1.12.6**, not 1.13.0-beta3 | `Project.toml` declares `julia = "1.12.4"`; 1.12.6 satisfies it. The paper's beta is a benchmarking detail, not a requirement. Everything ran clean, so no pin was needed. |
| `--startup-file=no` on every invocation | The user's global `startup.jl` loads Revise and OhMyREPL, which emit precompile noise into the logs and are irrelevant to reproducibility. Discovered during the Stage 2 smoke test. |
| Separate project at `ddp/env/` for our own examples | Keeps the Stage 5/6 work from depending on, or writing into, the authors' `experiments/filterddp` environment. |
| Papers downloaded from arXiv | The brief said they would be at `papers/`; they were not present anywhere on disk. Fetched `2504.08278` and `2606.01487` and extracted text with `pdftotext -layout` (the PDF reader in use misreports these files as password-protected; they are not, and have no `/Encrypt` dictionary). |

## Corrections to my own work

- **`copper_plate_battery.jl`, Note 3.** I originally asserted that `w = 0` would make the
  reduced Hessian singular and the solution non-unique, and wrote a probe to demonstrate it.
  The probe **refuted the claim**: with `w = 0` the solve still converged in one iteration
  with zero regularisation and returned the same trajectory from two very different
  starting guesses. The reason is that the balance equality ties `pb` to `pg`, so the null
  space of `∇_u c = [1,1]` is spanned by `(1,-1)`, along which the cost has curvature
  `2a(τ) > 0`. The note and the probe's commentary were rewritten to state the negative
  result. `w > 0` is about inter-period coupling and physical sense, not well-posedness.
- **Stage 6, `T = 6` all-bounds configuration.** My first choice of bounds
  (`pg ≤ 1.25`, `pb ≤ 0.4`) was **infeasible**: period 4 has `d = 1.8` but
  `pg + pb ≤ 1.65`. The reference solver correctly reported no KKT point; I had it error
  out rather than return `nothing`. Fixed both: the reference now returns `nothing` for an
  empty feasible set, the `T = 6` case uses feasible bounds, and the infeasible
  configuration was kept deliberately as probe **6g** because FilterDDP's behaviour on it
  is informative.
- **Objective changed to the tADMM paper's form (2026-08-04).** Stages 5-6 originally used
  `Σ a^t(P_Subs^t)² + b^t P_Subs^t`, with `a^t > 0` chosen purely to get strict convexity
  and a closed-form reference. That is not the objective the MPOPF work actually uses, and
  it was not needed: the paper's `Σ c^t P_Subs^t Δt + C_B Σ (P_B^t)² Δt` is strictly convex
  on the feasible set through `C_B` alone. Both example files, both reference solutions and
  the model section were switched over. **This changed a reported result:** under the old
  objective Stage 6d exited status 7 and the tolerance sweep failed below `~1e-9`; under
  the new one it converges in 16 iterations and certifies to `1e-12`. The earlier "accuracy
  floor" finding was an artefact of my objective choice and has been corrected in
  `README_FILTERDDP_EXPERIMENT.md` and `MPOPF_APPLICABILITY.md` §4 and §7.
- **Round-trip efficiency `η` removed entirely (2026-08-04).** The storage dynamic was
  written `B^{t+1} = B^t − η·P_B^t·Δt` with `η = 1`. A single `η` applied to a *signed*
  `P_B` cannot be right in both directions: charging stores `η_c·|P_B|`, but discharging
  `P_B` into the bus must drain `P_B/η_d`, and `η` and `1/η` coincide only at `η = 1`.
  Carrying the symbol at unity was an invitation for someone to later set it to `0.95` and
  get a quietly wrong model, so it was dropped and the battery is now explicitly lossless —
  matching the tADMM paper, which is also lossless. **No numerical result changed**, since
  `η` was `1.0` throughout. The paper gained a paragraph on why, and a remark on what a
  proper two-efficiency model would cost in the DDP setting: `P_c, P_d ≥ 0` map onto native
  control bounds and keep the rank condition comfortably (`n_u = 3` vs `n_c = 1`), but
  non-simultaneity `P_c·P_d = 0`, if ever imposed explicitly rather than left implicit,
  becomes a complementarity constraint whose Jacobian loses rank exactly at the solutions
  of interest.
- **`Δt` is now carried symbolically in the code.** It previously existed only in the
  model section, with the code hardcoding `Δt = 1`; the two agreed only because the
  instance sets `Δt = 1`. Closed.
- **Bug in my own reference solver, found by a nonsense result.** The active-set
  enumeration accepted a subset whose KKT matrix is *singular* — e.g. `Psub[1] ≤ 1.35` and
  `Psub[1] ≥ 0` both active. Julia's `\` does not always throw on those; it returned a
  least-squares vector that happened to be primal and dual feasible, so it was accepted as
  the optimum. This produced a "reference" objective of `317.555` against FilterDDP's
  correct `3.676` at T = 6 — i.e. **the solver was right and my reference was wrong.**
  Fixed by checking that the KKT system was actually solved (`‖K·sol − rhs‖` small), not
  merely that its output looks feasible.
- **Stage 6 bound sets corrected.** The T = 3 "all bounds" case combined the 6a and 6b caps,
  which is infeasible: `Psub ≤ 1.35` with `P_B ≤ 0.10` caps period 2 at `1.45 < p_L = 1.6`.
  Each bound worked alone only because the other variable was free to cover the peak.
  Loosened so the combined set is feasible while still binding. The T = 6 all-bounds case
  had no active constraints at all and was tightened so it tests something.
- Two PowerShell quoting bugs of mine (not repository issues) produced spurious failures
  that are visible in the logs and then corrected in place: the Stage 2 smoke test lost its
  string quotes, and `Pkg.test("FilterDDP")` lost its argument to native-command argument
  parsing. Both were rerun correctly; `logs/tests.log` contains the failed attempt followed
  by the corrected one.

## Things deliberately not done

- No feasibility restoration phase was added, no penalty was substituted for FilterDDP's
  method, and no convergence tolerance was weakened to turn a failure into a pass. Stage 6d
  and 6g are reported as failures because they are failures.
- No distribution-network architecture was built. Stopped at the copper-plate model, as
  instructed.

## Network-scale local fork changes (2026-08-21--22)

The later network evaluation necessarily changed the local ignored clone. The complete
diff against upstream commit `513a104` is tracked as
`ddp/patches/dynamic_network_scaling.patch`; apply it from inside the clone with
`git apply ../../patches/dynamic_network_scaling.patch`.

The patch carries the previously recorded per-stage objective/constraint support, replaces
large statically sized trajectory and solver fields with ordinary dynamic vectors and
matrices, and erases generated-function types from `Solver` to avoid compiler recursion on
network models. It also contains the later sparse-KKT follow-up described below, so applying
this one patch reproduces the final working fork rather than only the intermediate failure.

The initial dynamic-storage subset let the IEEE2522C `T=3` model and solver construct, but
did not make FilterDDP sparse. At that intermediate revision the backward pass still built
a 13,358-by-13,358 dense Hessian, performed dense QR on a 13,358-by-10,337 equality
Jacobian, and materialised the full orthogonal matrix. The guarded run therefore completed
zero iterations. This remains reported as an intermediate architectural scalability limit,
not retroactively presented as a successful solve.

## Sparse network KKT follow-up (2026-08-22)

The dense limitation was subsequently removed with sparse derivative callbacks and a
direct sparse saddle-point solve in the backward pass. IEEE2522C at `T = 3` now completes
56 FilterDDP iterations in 233.098 s, matches the centralized objective within `7.342e-4`,
and satisfies the equality equations to `5.631e-10`. The implementation and successful
result are recorded separately so the original dense failure remains visible in the
research chronology.

## Factored bound sensitivities (2026-08-31)

The sparse path originally retained two dense `nu x nx` bound-dual maps,
`ζl = -Σ_L .* beta` and `ζu = Σ_U .* beta`, at every horizon stage. These are
exact row-scaled copies of the already-retained feedback matrix `beta`. The
follow-up patch stores only the two length-`nu` scale vectors and evaluates
their action using the `beta*δx` product already required by rollout.

Apply `ddp/patches/factor_bound_sensitivities.patch` after
`dynamic_network_scaling.patch`. It reduces the measured large10k update-rule
storage from 1606.982 MiB to 757.012 MiB per stage while reproducing the
IEEE123C `T = 3` iteration count, objective, and residual exactly. Full details
are in `ddp/notes/FILTERDDP_FACTORED_BOUND_SENSITIVITIES.md`.
