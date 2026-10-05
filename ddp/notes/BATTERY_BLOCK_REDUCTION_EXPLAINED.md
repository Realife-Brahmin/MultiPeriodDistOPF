# The battery-block (Kron / Schur-complement) reduction, explained

Agenda item of 2026-10-07: *which variables are eliminated, which are kept,
how the DDP updates are recovered, and is it exact?* Short answers first.

- **It is exact for the stage matrix it is given.** It is the same linear
  solve, done in a different elimination order. Nothing is dropped.
- **It is a separate thing from the Hessian approximation.** The diagonal
  Hessian decides what the stage matrix *is*; the reduction decides how that
  matrix is *solved*. Either can be changed without the other.
- **The idea needs no special topology; the fast way of doing it does.** The
  leaves-to-substation sweep is for radial networks.

## 1. What one stage has to solve

At stage `t` of a backward sweep, FilterDDP solves

```
K [ α  β ]  =  - [ Q̂u   B̃  ]        K = [ Ĥ    c_uᵀ ]
  [ ψ  ω ]       [ c    c_x ]            [ c_u   0   ]
```

- Unknowns: one column `[α; ψ]` (the step in every control `u^t` and the
  change in the constraint multipliers) and `n_x` columns `[β; ω]` (how
  every control and multiplier responds to each battery's energy `B^t`).
- `u^t` holds every network variable of the period: `P_Subs`, `Q_Subs`,
  line `P`, `Q`, `ℓ`, voltages `v`, battery powers `P_B`, DER reactive
  controls, and the slacks. At large10k that is 54,665 controls and 42,303
  multipliers: `K` is 96,968 x 96,968, and there are 1 + 1,020 columns.

## 2. Which unknowns are kept, which are eliminated

Split the rows/columns of `K` into

- **kept, `E`** (2 `n_B` of them): the battery powers `P_B` and the
  multipliers of the battery energy rows;
- **eliminated, `F`** (everything else): line flows, currents, voltages,
  substation and DER variables, slacks, and the multipliers of the balance,
  voltage-drop and SOC rows.

Why this split:

1. The dynamics are `B^{t+1} = B^t - Δt P_B`, so `f_u` touches only `P_B`.
   Hence `B̃ = f_uᵀ V_xx f_x` (plus the terminal `ℓ_ux`) is nonzero only in the
   `P_B` rows, and `c_x` only in the energy rows: **the `n_x` feedback
   right-hand sides are supported on `E`.**
2. The value update reads only `E` rows of the solution:
   `V_xx^t = C + β_Bᵀ B̃_B + ω_Eᵀ c_x,E` uses the `P_B` rows of `β` and the
   energy rows of `ω`; and `βᵀQ̂u + ωᵀc`, needed for `V_x`, equals
   `B̃ᵀα + c_xᵀψ` because `K` is symmetric.

So of the 96,968 x 1,020 feedback block, 2,040 x 1,020 (2%) is used.

## 3. The elimination

Write `K x = b` in the two blocks:

```
K_FF x_F + K_FE x_E = b_F
K_EF x_F + K_EE x_E = b_E
```

Solve the first for `x_F = K_FF⁻¹ (b_F - K_FE x_E)` and substitute:

```
S x_E = b_E - K_EF K_FF⁻¹ b_F,        S = K_EE - K_EF K_FF⁻¹ K_FE.
```

`S` (2 `n_B` x 2 `n_B`) is the Schur complement of `K` onto `E`: the stage
network *Kron-reduced onto the batteries*. It says how each battery's power
and energy-row multiplier respond to every other battery's, with the whole
network already accounted for. This is ordinary block Gaussian elimination;
the only requirement is that `K_FF` is nonsingular (it is: the barrier terms
make it so, and every local factorization would report otherwise).

## 4. How the DDP quantities come back

| quantity | how |
|---|---|
| battery rows of the feedback, `[β_B; ω_E]` | `S⁻¹ (-[B̃_B; c_x,E])` — the right-hand side has `b_F = 0` |
| `V_xx^t` | `C + β_Bᵀ B̃_B + ω_Eᵀ c_x,E` |
| step `[α; ψ]` (all controls) | `x_E` from the reduced system, then `x_F = K_FF⁻¹ (b_F - K_FE x_E)` |
| `V_x^t` | `ℓ_x + f_xᵀ V_x + B̃ᵀ α + c_xᵀ ψ + c_xᵀ φ` |
| forward pass, `β δx` | one full solve with right-hand side `-[B̃ δx; c_x δx]` |

Network quantities are **not** approximated or removed: every flow and
voltage of the step is recovered by the back-substitution for `x_F`.

## 5. Is it exact? Evidence

- Algebraically: yes, for the `K` given (Section 3).
- Numerically: the battery rows agree with a full UMFPACK solve to 5e-16
  (ieee123), 3e-12 (med2522), 8e-16 (large10k); `S` agrees with MUMPS's Schur
  complement to 1e-16 .. 3e-11; one-column solves have residuals no larger
  than UMFPACK's.
- In the algorithm: FilterDDP reaches near-optimality at the same iteration
  with the same objective on every tested case (one borderline case moves by
  one iteration), and run to strict tolerance it reproduces UMFPACK's result.

## 6. How this differs from the Hessian approximation

Two independent choices:

| | full solve of `K` | reduced solve of `K` |
|---|---|---|
| **exact stage Hessian** | original FilterDDP | not built (see below) |
| **diagonal stage Hessian** | Table II baseline | tree solver |

- The **diagonal Hessian is an approximation of `K` itself**: it drops the
  off-diagonal part of the battery curvature `f_uᵀ V_xx f_u` (and the `v`-`ℓ`
  cross terms of the SOC rows). It changes the iterates: 0-37% more
  iterations, and it needs `C_B > 0`.
- The **reduction is an exact way of solving whichever `K` was built.** It
  does not change the iterates.

They also interact favourably: `f_uᵀ V_xx f_u` lives entirely in the
`P_B` x `P_B` block, i.e. inside `K_EE`. An exact-Hessian version would
therefore change only the battery block `S`, not the network elimination.
Not built; the feeder-block split of Section 7 assumes the diagonal form.

## 7. Doing it fast on a radial network

Applying `K_FF⁻¹` is the expensive part in general (it is why MUMPS's Schur
complement was not faster). On a radial network it is cheap:

- Group the unknowns by bus: the incoming line's `P, Q, ℓ`, SOC slack, the
  bus voltage, its DER control, and its balance, voltage-drop and SOC rows
  (at most 10 unknowns). Each group couples only to its parent, through 3
  unknowns (the parent's two balance rows and voltage). Checked against the
  sparsity pattern on all three systems.
- Eliminate groups from the leaves to the substation: a 10 x 10 factorization
  per bus and a 3 x 3 update to the parent, no fill.
- Feeders meet only at the substation, so
  `S = D - U M⁻¹ Uᵀ` with `D` block diagonal by feeder and `M` the 6 x 6
  substation block: Woodbury, one small solve per feeder.

## 8. What it is not

- Not the earlier *reduced-space* method (FilterDDP on battery variables with
  an inner OPF). That changed the optimization problem and its iterates. This
  changes neither.
- Not a network equivalent or model reduction: the full network is solved.

## 9. Limits

- The fast sweep needs a radial network. With a cycle the structure check
  fails and the code stops; it does not return a wrong answer. A few loops
  could be handled by keeping the loop-closing unknowns alongside `E`
  (untested).
- ieee123 is slower with it (fixed overhead exceeds UMFPACK's 4 ms).
- Retained storage per stage is not yet lean (PARALLEL_IN_TIME.md).
