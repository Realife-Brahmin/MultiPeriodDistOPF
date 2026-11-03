# Rated campaign, med3980 S1-S3 (2026-11-02/03)

Per-line ampacity limits (`FB_RATINGS=1`, `test_networks/ieee9500xf_1ph/BranchRatings.csv`), 16 directions,
T = 24, PV reactive power dispatched. Kamiak `~/msopf/jobs/farkas_rated_nov02`; `farkas_job.sh`, `submit.sh`.

| phase | what | result |
|---|---|---|
| A | idle edge (OPF) | equal to the unrated campaign in all 16 directions: ratings never bind without batteries |
| C | idle certificate | certified at 1.06-1.08x the edge in all 16 |
| AD | batteries over the day (OPF) | shrinks in 4 directions (90 deg 9.52->8.59, 112.5 7.43->7.23, 225 13.01->12.35, 337.5 7.96->7.64); the rest within 0.01 deg |
| CD | dispatched certificate | all 16 certified: 1.02x in the four rating-bound directions (relaxation INFEASIBLE, a Farkas ray), 1.032x in direction 1, 1.05x elsewhere |

CD ran with `FB_CD_BISECT=0` (the 1.02/1.05 grid, as in the unrated CD) except direction 1, which bisected.

Direction 10: on MA57 its hour-3 solves near 1 deg returned LOCALLY_SOLVED but failed the balance check, so A
came out 1.0 deg (`dir10_ma57_artefact/`). A and C were rerun with MA97 on the 309 PC (`dir10_ma97_309pc_run.sh`,
commit 335a30e3; ieee9500xf ratings are unchanged by the switch rule): A 5.5603 deg, C certified at 1.0625x.
AD direction 10 was started from the artefact idle edge (column `idle` = 1.0037); its edge 7.0878 deg is a separate
search over the day and matches the unrated campaign.
