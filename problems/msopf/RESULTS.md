# Paper-matched MSOPF results

These are the rated-current, T=24 cases used by the PESGM 2027 study. Begin with the [Farkas map viewer](https://realife-brahmin.github.io/MultiPeriodDistOPF/farkas_maps.html) or its [released HTML snapshot](../../docs/farkas_maps.html), then inspect the numerical evidence below. [Release provenance and validation coverage](RELEASE.md) apply to these links.

## IEEE123, S1/S3, 4.16 kV

[Campaign explanation](results/ieee123_416_nov03/README.md) and [derived-bound results](results/ieee123_416_nov03/derived_bounds/ieee123_s13/): two directions, with fb_A/fb_AD edge and reach files and fb_C/fb_CD certificates. Older 20 kV results and earlier fleet/rating trials remain historical evidence; they are not the paper's final IEEE123 case.

## med3980, S1–S3

[Rated campaign](results/kamiak/farkas_rated_nov02/README.md), [edges, battery reach and whole-day certificates](results/kamiak/farkas_rated_nov02/med3980_s123/), and [tighter idle certificates with derived bounds](results/kamiak/farkas_rated_nov02/med3980_s123_derived/): 16 directions. The campaign README identifies the superseded direction-10 MA57 artefact and the MA97 rerun; retain that distinction.

## sf6937, S1–S4, SMART-DS day 148

[Campaign explanation](results/kamiak/farkas_day148b_nov03/README.md) and [all per-direction results](results/kamiak/farkas_day148b_nov03/sfp1u_four/): 78 directions. fb_A contains idle edges, fb_Cderived the paper's tighter idle certificates, fb_AD the daily battery reach, fb_CD the numerical whole-day certificates, and fb_VD the post-run whole-day audits. [Exact-check coverage](results/kamiak/farkas_day148b_nov03/reference_config/cd_exact_status.csv) records 78/78 for fb_CD (59 in-run and 19 by VD), including the [six follow-up audits](results/kamiak/farkas_day148b_nov03/reference_config/vd_audits_nov09/README.md).

The [older polytope/region viewer](../../docs/four_regions.html) is supporting historical material. Its sensitivity predictions, sparse screens and older test configurations should not be read as the current rated certificate campaign or as globally certified polytopes.
