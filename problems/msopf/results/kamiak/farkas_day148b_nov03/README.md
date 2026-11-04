# sf6937 rated campaign, day 148 (2026-11-03)

Kamiak ~/msopf/jobs/farkas_day148b_nov03 (repo_rated at 335a30e3, later 9e31b135 / 36dbdf58 / 5712c808 / 4f29d098 /
865e26d1 for the audits and reruns), with the 309 PC's directions 69-78 (local_worker.sh). SMART-DS customer profiles
of day 148 (FB_LOADDAY=148), PV from sunrise to sunset, every line at its rating (BranchRatings.csv, 335a30e3: a switch
rated below the conductor in series takes the conductor's rating).

| files | phase |
|---|---|
| fb_A, fb_Acur | no-battery edge (78) |
| fb_C | no-battery certificate, ratings alone (78; 17 directions rerun with the exact check, earlier results in Kamiak C_unverified/) |
| fb_V | exact audit of fb_C (re-solve; 61 re-verified, the other 17 rerun) |
| fb_Cderived | no-battery certificate with derived current bounds + lsec (FB_DERIVED=1, FB_CTAG=d; on Kamiak fb_Cd, renamed here: Windows filenames are case-insensitive and fb_Cd_<i> is fb_CD_<i>) |
| fb_AD, fb_ADcur | battery reach over the day (78) |
| fb_CD | battery certificate (FB_CD_BISECT=0, grid 1.02/1.05/1.10; 39 of 78 as of 19:50, the rest running) |
| fb_VD | exact audit of the CD certificates taken before the exact check (13 of 13 passed) |
