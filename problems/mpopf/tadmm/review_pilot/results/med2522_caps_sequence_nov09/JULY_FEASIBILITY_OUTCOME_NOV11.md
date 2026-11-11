# Med2522 T96 July-input outcome

The flat-start retry 31681491 reached centralized MA57 and exited with
LOCALLY_INFEASIBLE after 369 Ipopt iterations and 406.49 s inside
Ipopt. The unscaled constraint violation was 0.0185262 in the model
convention. This was the centralized warmup. There is no measured
centralized benchmark time or tADMM trajectory. Slurm cancelled the
dependent dispatch audit 31681492.

An independent voltage bound confirms infeasibility of the corrected
model on the frozen 4.1569 kV July release inputs. For each period,
the bound grants every battery maximum discharge with no energy
constraint and every inverter maximum reactive injection at its fixed
PV output. Summing downstream net demands gives lossless lower bounds
on active and reactive branch flow. In the lossy equations, write each
branch flow as that subtree's net demand plus its own and downstream
series losses. Substituting into the voltage drop leaves a negative
own-line term -(r^2+x^2)*ell and nonpositive contributions from
downstream losses. This follows from nonnegative r, x and ell, with
ell also implied nonnegative by the cone and positive parent voltage.
Discarding these losses therefore gives an upper bound on squared
voltage, even if a subtree exports power.

At bus 553 in period 24, the upper bound was 0.897704404788 p.u.
squared (0.947473 p.u.), below the required 0.9025 p.u. squared
(0.95 p.u.). The same calculation on the frozen current-input export
gave a positive minimum margin of 0.132047364305 p.u. squared. The
bound alone does not prove feasibility of the current case, whose
centralized solution independently passed its chronological replay.

The archive preserves the original solver and export logs, input and
source manifests, the initialization change, sampled memory, Slurm
accounting and the independent bound. SHA256:
6e876eeeaeba6073aa9a2dc8a38fc61ba00c7da056f89ef11a70d08ac6ab3c1f

The original release inputs remain unchanged. The result applies to
the corrected OPF equations and the archived voltage base. Altering
voltage, loads, or operating limits would define a separate experiment.
The current-input residual-gated diagnostic and both cap runs remain
archived separately. None provides a validated speedup.
