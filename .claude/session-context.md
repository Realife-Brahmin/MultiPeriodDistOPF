# Session Context — IAS Transactions Paper (2026-05-04)

## Paper location
`c:\Users\aryan\Documents\documents_general\documentsCreated\IAS-Trans-2025-Scaling-MPOPF-Computation-via-Temporal-Decomposition\`

## Results location
`c:\Users\aryan\Documents\documents_general\MultiPeriodDistOPF\envs\tadmm\processedData\`

---

## Current paper state

### What was done this session
- **`sections/simulation.tex`** — fully rewritten with real data (3 commits on master)
- **`main.tex` abstract** — updated to reflect actual speedup results (removed red-flag placeholder)

### Sections NOT yet touched (user not confident in old content)
- `sections/introduction.tex`, `sections/theory.tex`, `sections/mpopf_formulation.tex`
- **Do not edit these without being asked.**

---

## Three test systems (reader-facing names — no A/B/C suffixes ever)

| System | Buses | Branches | Batteries | PV |
|--------|-------|----------|-----------|-----|
| ieee123 | 128 | 127 | 51 | 51 |
| medium2552 | 2,522 | 2,521 | TBD | TBD |
| large10k | 10,321 | 10,320 | 1,020 | 1,022 |

Internally all "C" variants — readers never see this.

### Result availability
- **large10k**: results confirmed for T=6,12,24,48 ✓
- **ieee123**: no confirmed results yet (ieee123C T=96 hit max iters, didn't converge); planned T∈{24,48,96,144}
- **medium2552**: no results yet; planned T≥12

---

## large10k confirmed results (Ipopt/MUMPS, 16 threads, SOCP/BFM-NL)

| T | Δt(h) | BF Time(s) | BF Status | tADMM Iters | tADMM Eff(s) | Speedup | Near-Opt(s) | N-O Speedup | Obj Diff |
|---|-------|-----------|-----------|-------------|-------------|---------|-------------|-------------|---------|
| 6 | 4.0 | 664.3 | OK | 37 | 594.3 | 1.12× | — | — | <0.001% |
| 12 | 2.0 | 1,534.9 | OK | 35 | 2,161.4 | 0.71× | 1,113.9 | 1.38× | 0.002% |
| 24 | 1.0 | 4,332.0 | OK | 58 | 5,338.9 | 0.81× | 4,337.0 | 1.00× | 0.005% |
| 48 | 0.5 | 17,479 FAILED | OOM | 79 | 8,398.5 | 2.08× | 5,709.3 | 3.06× | ~0.40% |

T=48 BF: Ipopt OOM at iter 1,444; last feasible obj $2,981,032.92  
T=48 tADMM (rho₀=30,000): converged, obj $2,992,815.54, feasible

tADMM params: ρ₀=500 (T≤24), 30,000 (T=48); ε_pri=1e-4 (T=6), 2e-3 (T≥12); ε_dual=1e-2

---

## What still needs doing
1. Run ieee123 experiments (T=24,48,96,144) and insert results
2. Run medium2552 experiments (T≥12) and insert results
3. Review/rewrite intro, theory, formulation sections
4. Refine abstract once all results are in

---

## Practical notes
- Paper compiles cleanly: `latexmk -pdf` from paper directory
- Pre-existing overfull hbox warnings in theory.tex lines 37,53 — not our problem yet
- Global Claude permissions: `bypassPermissions` in `~/.claude/settings.json` (needs session restart)
- Result files: `envs\tadmm\processedData\<system>_1ph_T<N>\results_socp_{bf,tadmm}.txt`
