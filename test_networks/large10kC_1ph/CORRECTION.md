# large10kC electrical correction

The active feeder is a synthetic one-phase aggregate equivalent. On
2026-08-21 its source data were corrected before using it as a centralized or
FilterDDP benchmark.

The original conversion declared a 12.47-kV line-to-line system but wrote
7.1996 kV (12.47/sqrt(3)) into the one-phase source and all aggregate devices,
while retaining aggregate powers. The OPF parser additionally used a 2.4018-kV
base. The corrected aggregate equivalent uses 12.47 kV and 1 MVA throughout.
Aggregate load, PV, and storage powers are unchanged.

The original synthetic network also assigned every ordinary line section
`r1=0.07` ohm and `x1=0.01` ohm. Even after the voltage-base correction, the
exact SOCP and LinDistFlow models were infeasible under the stated 0.95--1.05
p.u. band. LinDistFlow with voltage bounds omitted reached a minimum voltage of
0.762 p.u. A uniform half-impedance diagnostic remained infeasible, while a
quarter-impedance diagnostic was feasible. The active synthetic line sections
therefore represent quarter-length sections:

- 10,219 ordinary sections: `(r1,x1)=(0.0175,0.0025)` ohm;
- 101 area-boundary sections: `(r1,x1)=(0.000025,0.000025)` ohm.

No voltage limits were relaxed, no powers or DER capacities were changed, and
the loss terms and exact second-order-cone constraint remain active. Cold-start
native-conic Gurobi solves are optimal and pass the independent validator at
both `T=1` and `T=3`; detailed results are recorded under
`ddp/results/network_filterddp/`.
