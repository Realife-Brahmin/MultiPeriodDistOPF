# Battery with charging and discharging efficiencies (2026-11-15)

Reviewer items AE-4 and R2-3. `battery_efficiency_probe.jl`, serial, 309 PC,
July battery coefficient. Model: P_B = P_dis - P_ch, both nonnegative,
B_t = B_{t-1} + dt*(eta*P_ch - P_dis/eta). Simultaneous charging and
discharging is not forbidden; it is measured on the result.

In the storage-coordinated scheme this model appears in Update II only (one
small convex problem per battery, with a 1e-6 charge on throughput to choose
among equally close schedules). The period OPFs, the penalty rule and the
repair are unchanged. The reference is the centralized problem with the same
battery model, as an explicit cone program solved by Clarabel ("almost
optimal" in every case here).

| case | eta | first iteration with a repaired dispatch within 0.1% | lossless model | largest simultaneous charge and discharge (p.u.) |
|---|---|---|---|---|
| ieee123 T24 | 0.95 | 11 | 10 | 2e-8 |
| ieee123 T24 | 0.90 | 11 | 10 | 1e-8 |
| ieee123 T96 | 0.95 | 11 | 10 | 2e-8 |
| med2522 T24 | 0.95 | 9 | 9 | 3e-8 |

The repaired cost is within 3e-7 of the centralized one by iteration 40 on
ieee123. No battery charges and discharges at once, in the centralized
solution or in any logged tADMM iterate, although nothing forbids it.

Limits of this probe: iteration counts only, no timing; the probe has no pinned
repair, so on med2522 the repair is again refused at later iterations (for
example iteration 30); med2522 at other horizons and large10k were not run;
no OpenDSS replay, since the replay runs storage at 100% efficiency in its
"model" mode and would need its "dss" mode extended to these efficiencies.
