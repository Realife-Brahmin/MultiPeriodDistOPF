# Untimed value-representation inspection

See `../../notes/VALUE_REPRESENTATION_ASSESSMENT_OCT21.md` for interpretation,
assumptions and proposed next test. No compressed backend is implemented here.

Run with Julia 1.12.7 and standard libraries only:

```
julia --startup-file=no inspect.jl INPUT_DIRECTORY OUTPUT_DIRECTORY
```

INPUT_DIRECTORY contains stage_1.jls through stage_4.jls from the unchanged
original-1020 first-pass capture at
`/home/aryan.jha/kamiak_jobs/stage_b_replay_oct15/capture`, and original-battery
`network_data_large10kC_1ph_T96_periodic.jls` from that capture's source checkout.
The T3 metadata is the original (not X/Y) file from
`/home/aryan.jha/kamiak_jobs/battery_count_oct21_campaign/inputs`.
The script asserts matching Bset order and parent maps before using groups.
Four stage hashes were independently checked against the historical capture
manifest. All six copied input hashes are retained in input.sha256.

The script uses the existing `pint_combine`, does not run an optimizer, and
does not time anything. Rank thresholds are 1e-12, 1e-8, 1e-6 and 1e-4 times
the selected submatrix's largest singular value. In the CSV rank_1e12 means
threshold 1e-12, not 1e12. These submatrices are exceptionally small relative
to the full matrices; numerical rank alone is not a physical coupling measure.

Output covers two individual periods, a two-period composition and a four-period
composition, first optimization pass only. No later-pass evidence, speedup,
solver-feasibility claim or proof of fixed-rank closure is supplied. Binary
captures remain in their existing cluster locations; no HSL files are included.
