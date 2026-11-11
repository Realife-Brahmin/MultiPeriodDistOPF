# Med2522 T96 with July release inputs

The independent comparison uses the feeder and shape files from tag
260731-ias-trans-tadmm, commit 50b8901c324ba9eefd0d672c6e484e0ededc7a57.
The release feeder directory is rawData/ieee2552C_1ph. Its DSS files are
copied byte-for-byte into the isolated campaign under the current harness
system identifier. The export uses the release feeder's 4.1569 kV source,
endpoint-inclusive July load and price sampling, July daylight PV
sampling, 96 quarter-hour periods, and July C_B. The current inactive
substation-device policy is applied at export. The archived files are
unmodified. There are 2522 buses, 2521 branches, and 249 modeled batteries
and PV units.

The corrected OPF equations and exact MA57 local solves are shared with
current-input job 31681398. Both use 96 single-thread workers on one
128-core Xeon 6760P allocation, rho0=20000 and the same residual-gated,
primal-first legacy_phase implementation. The centralized solve is
fully converged. The comparison changes the input package as a group,
including feeder, voltage base and profile sampling. It is neither a
literal reproduction of the July solver nor a one-input-factor ablation.

Campaign: /home/aryan.jha/kamiak_jobs/tadmm_review_nov11_med2522_july_inputs.
Export 31681413 passed. Solve 31681415 waits for the current-input solve
to finish. Independent chronological and SOCP audit 31681416 follows the
July-input solve and uses the July feeder. It records model-matched
lossless storage replay separately from solve timing.

Release archive SHA256:
eca7a8b4b8a11ae62b5fcb4f84d0ec76f7ce3889a0806c3a840aaee8ac4bef76

Serialized export SHA256:
f8cac0d191ff9dbe0684b7dcdb93bcfb45966892c7a9816ab0bbb1d82a89ac11
