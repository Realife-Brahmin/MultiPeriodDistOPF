# KKT matrix updates: magnitudes, stale Jacobians, block-diagonal Jacobians, HSL in Ipopt

Agenda for the meeting with R. Gupta on 2026-10-02. Branch
`ddp-kkt-matrix-updates-sep27`. Every experiment is on the paper's current
formulation unless stated:
- `T = 6`, periodic profile
- per-system `C_B`
- soft terminal SOC
- diagonal stage Hessian with the exact assembly rewrites

The companion study of how the matrix changes between iterations is
`KKT_EVOLUTION.md`.

## 2a. Are many KKT entries tiny?

**Setup.**
- **Captures.** Stage-1 KKT systems of all three systems at iteration 5, at
  iteration 40, and near the Table II near-optimality iteration (ieee123 67,
  med2522 68, large10k 100). Script: `capture_kkt_magnitude.sh`; the captures
  are gitignored.
- **Analysis** (`kkt_magnitude_analysis.jl`, driven by
  `run_kkt_magnitude.sh`):
  - The magnitude of every stored entry, absolute and after symmetric scaling
    `r_ij = |K_ij| / sqrt(d_i d_j)`, `d_i = max_k |K_ik|`.
  - For thresholds `1e-12 .. 1e-2`, every off-diagonal entry below the
    threshold is dropped (the diagonal is always kept). The thresholded matrix
    is then factored and solved, and the analysis reports:
    - the entries dropped
    - the factor fill
    - how far the `(n_x+1)`-column solution moves
- **Outputs.**
  - `ddp/results/kkt_magnitude/magnitude_threshold.csv`
  - `magnitude_blocks.csv`
  - `figures/kkt_magnitude_<system>_iter<k>.png`, which shows the full matrix,
    the matrix thresholded at `1e-4` and `1e-2`, the entries dropped at `1e-2`,
    and the magnitude histograms.

**What changes is the "primal" solution.** That is the control rows, from which
FilterDDP takes its feedforward and feedback gains. The multiplier rows are
far more sensitive and are reported in the CSV.

| threshold (absolute) | entries dropped | factor size | change in the control directions |
|---:|---|---|---|
| `1e-8` | 0.0-0.8% | -0.1% to -0.6% | <= 9e-5 (mostly <= 1e-6) |
| `1e-6` | 1.7-5.5% | -1.5% to -5% | ieee123: up to **13%** (iteration 40); med2522 <= 0.07%; large10k <= 0.04% |
| `1e-4` | 3.9-16% | -2% to -20% | ieee123 13%, med2522 up to **41%** (iteration 40), large10k <= 0.2% |
| `1e-2` | 16-37% | -19% to -55% | 2-43% everywhere |

**Where the small entries are.** They are the line-impedance coefficients:
- in the voltage-drop rows, `2r·P`, `2x·Q` and `(r²+x²)·ℓ`. At large10k 40% of
  those entries are below `1e-4` and 60% below `1e-2`.
- the `r·ℓ` and `x·ℓ` loss terms of the power-balance rows.
- a few Hessian diagonal entries sitting at the `1e-8` floor or with a tiny
  barrier term.

The SOCP and battery-energy rows are `O(1)`. Line impedances are small in
per-unit; that does not make them unimportant.

**Scaling is no guide.** Dropping entries that are small *relative to their
row and column scale* makes the matrix **singular** in late iterations, even
at a threshold of `1e-8` (ieee123 and med2522). Some constraint rows are tied
to their variables only through entries that look negligible next to a
barrier-dominated diagonal of up to `1e10`. Those entries are structurally
essential.

**Answer.**
- **Harmless:** only entries below about `1e-8` in absolute value can be
  dropped harmlessly, and that is under 1% of the matrix, so the sparsity plot
  does not change.
- **Useful but instance-dependent:** thresholds that visibly thin the matrix
  (`1e-4` .. `1e-2`) change the Newton directions by 13-43% on ieee123 and
  med2522.
- **large10k is the exception.** There `|K_ij| < 1e-3` removes 26% of the
  entries and 45% of the factor while moving the directions by only 2%. That
  would make an approximation worth a full run, but it cannot be assumed on
  other feeders.
