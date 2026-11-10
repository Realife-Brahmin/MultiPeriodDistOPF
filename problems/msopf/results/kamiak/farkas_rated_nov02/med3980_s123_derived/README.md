# med3980 S1-S3, no-battery certificates with derived current bounds (2026-11-03, 309 PC)

Phase C of the rated campaign (phase A from `../med3980_s123`) rerun with `FB_DERIVED=1 FB_CUTS=kcl0,jbound,lsec`
(commits 865e26d1 / e82585f2): every line's and source's current limit is min(rating, the bound the model's own
response implies at the point's angles). Certified at 1.0056-1.0094x the edge in all 16 directions (1.06-1.08x with
the ratings and kcl0, jbound alone), every accepted point exact-checked. About 22 min per direction, 6 in parallel.
`run.sh` / `med_par.sh`: the launchers.
