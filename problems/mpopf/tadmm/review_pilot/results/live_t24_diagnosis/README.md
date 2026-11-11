# T=24 warm-up outcome plots

Snapshot: November 7, 2026, 9:35 a.m. Pacific. Numerical source and options are
unchanged. These are the initial convergence passes, not completed timing or
independent OpenDSS-validation results. The archive freezes whole CSV lines,
solver logs, reference objectives, metadata and source/input identities.

`t24-warmup-outcomes-31643647/T24-convergence.png` shows both systems. The upper
row divides each residual by its actual stopping threshold; all three must be
at or below one. The bottom row shows the controller's normalized residual
ratio and its no-update band, checked every five iterations.

med2522 MA57 starts rho=3800, doubles it after iteration 10 to 7600, and never
changes it again. It settles into a repeatable odd/even objective/feasibility
oscillation and hits the 5000-iteration cap. All-copy residual is about 2.39
times its threshold; dual residual is about 17.24 times its threshold. The
relative-controller ratio is about 2.68, inside [0.1,10], so no penalty change
is triggered despite stagnation. `med2522-tail.png` shows the last 100 steps.
The dual residual is essentially constant; its insignificant 1e-11 fluctuations
are not magnified by the axis. Objective and original-feasibility alternation
are materially larger and plotted directly.

This is not runaway penalty growth. It is an observed plateau/two-cycle with a
fixed penalty after the early adjustment. The plots alone do not establish
whether the cause is local solve accuracy, the selected penalty, or another
numerical/implementation effect. Raising the iteration cap would not diagnose
it. A controlled saved-state accuracy check and a deliberately chosen fixed-rho
or relaxation comparison are more informative than another unchanged long run.
No new tuning sweep is launched by this report.

large10k MA57 completes warm-up at 793 iterations. Its rho updates are 48->96
after iteration 10, then 96->192 after 345, 192->384 after 355, and 384->768
after 560. Final original violation is 4.009e-7, consensus 2.495e-7 and dual
1.087e-5, all within the strict thresholds. The measured restart has begun;
its completion and independent OpenDSS checks remain pending.

The 0.1% objective / 1e-5 original-feasibility criterion is a separate posthoc
optimization measure: med2522 meets it from iteration 65 through its cap;
large10k first meets it at 168 and remains qualified from 403 through 793.
Those observations do not turn the failed strict med2522 run into a success,
and do not establish OpenDSS feasibility at intermediate iterates.

Reproduction: run `plot_live_t24.py` on the extracted archive directory. Python
syntax checks and visual inspection of both PNGs passed. Raw data is unchanged.

Archive SHA256:
`9bf10612aac3771f91729baa6ff1e90ac692cf2711d782285bae72be224337d5`
