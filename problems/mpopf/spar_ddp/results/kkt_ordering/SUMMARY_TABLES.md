
### ieee123C_1ph, diag Hessian: n = 1,353, nnz(K) = 5,269, RHS = 52 columns

| config | ordering used | ordering s | factor s | 1-RHS solve s | wide solve s | factor+wide s | vs baseline | fill (stored) | fill (ordering only) | wide residual | pivots |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| UMFPACK_default | AMD/COLAMD (UNSYMMETRIC) | 0.000602 | 0.001376 | 3.14e-05 | 0.0009068 | 0.002283 | 1.00x | 2.09 | 1.91 | 1.5e-13 | offdiag_pivots=-1 |
| UMFPACK_unsym_colamd | AMD/COLAMD (UNSYMMETRIC) | 0.000576 | 0.001358 | 3.19e-05 | 0.0009037 | 0.002261 | 0.99x | 2.09 | 1.91 | 1.5e-13 | offdiag_pivots=-1 |
| UMFPACK_sym_amd | AMD/COLAMD (SYMMETRIC) | 0.0004526 | 0.001465 | 3.95e-05 | 0.00105 | 0.002515 | 1.10x | 2.30 | 1.90 | 1.0e-14 | offdiag_pivots=213 |
| UMFPACK_unsym_metis | METIS (UNSYMMETRIC) | 0.004318 | 0.001683 | 4.7e-05 | 0.001267 | 0.00295 | 1.29x | 2.50 | 2.48 | 5.7e-14 | offdiag_pivots=-1 |
| UMFPACK_sym_metis | METIS (SYMMETRIC) | 0.002767 | 0.001804 | 5.39e-05 | 0.001321 | 0.003125 | 1.37x | 2.71 | 2.32 | 1.6e-14 | offdiag_pivots=198 |
| UMFPACK_cholmod | AMD/COLAMD (UNSYMMETRIC) | 0.001018 | 0.001397 | 3.1e-05 | 0.0008919 | 0.002289 | 1.00x | 2.06 | 1.90 | 4.6e-14 | offdiag_pivots=-1 |
| UMFPACK_best | AMD/COLAMD (UNSYMMETRIC) | 0.008045 | 0.001319 | 3.11e-05 | 0.0008872 | 0.002206 | 0.97x | 2.08 | 1.91 | 1.0e-13 | offdiag_pivots=-1 |
| MUMPS_sym_auto | AMF (sym=2) | 0.001244 | 0.01428 | 0.00229 | 0.01265 | 0.02694 | 11.80x | 5.35 | 3.25 | 3.0e-16 | negative=562 delayed=10 null=0 est_entries=13965 |
| MUMPS_sym_amd | AMD (sym=2) | 0.0009747 | 0.01791 | 0.00227 | 0.01242 | 0.03033 | 13.28x | 5.96 | 3.58 | 2.4e-15 | negative=562 delayed=66 null=0 est_entries=14770 |
| MUMPS_sym_amf | AMF (sym=2) | 0.001233 | 0.01421 | 0.00219 | 0.01226 | 0.02647 | 11.59x | 5.35 | 3.25 | 3.0e-16 | negative=562 delayed=10 null=0 est_entries=13965 |
| MUMPS_sym_qamd | QAMD (sym=2) | 0.000978 | 0.01793 | 0.00221 | 0.0124 | 0.03033 | 13.28x | 5.96 | 3.58 | 2.4e-15 | negative=562 delayed=66 null=0 est_entries=14770 |
| MUMPS_sym_pord | PORD (sym=2) | 0.002514 | 0.01511 | 0.00228 | 0.01237 | 0.02748 | 12.03x | 6.62 | 3.91 | 2.6e-15 | negative=562 delayed=23 null=0 est_entries=17044 |
| MUMPS_sym_scotch | SCOTCH (sym=2) | 0.002587 | 0.01613 | 0.00216 | 0.01198 | 0.02811 | 12.31x | 7.59 | 4.16 | 3.0e-16 | negative=562 delayed=34 null=0 est_entries=19366 |
| MUMPS_sym_metis | METIS (sym=2) | 0.002905 | 0.0152 | 0.00222 | 0.01237 | 0.02757 | 12.08x | 6.74 | 3.93 | 2.2e-16 | negative=562 delayed=26 null=0 est_entries=17378 |
| MUMPS_sym_metis_compressed | METIS (sym=2) | 0.002887 | 0.01552 | 0.00221 | 0.0123 | 0.02783 | 12.19x | 6.74 | 3.93 | 2.2e-16 | negative=562 delayed=26 null=0 est_entries=17378 |
| MUMPS_sym_ipopt | AMF (sym=2) | 0.00126 | 0.0133 | 0.00227 | 0.01251 | 0.02581 | 11.30x | 5.30 | 3.25 | 4.6e-12 | negative=562 delayed=0 null=0 est_entries=13965 |
| MUMPS_sym_amd_plain | AMD (sym=2) | 0.0004275 | 0.01746 | 0.00231 | 0.01263 | 0.0301 | 13.18x | 5.86 | 3.27 | 4.3e-16 | negative=562 delayed=109 null=0 est_entries=14019 |
| MUMPS_sym_amf_plain | AMF (sym=2) | 0.0005182 | 0.01581 | 0.00224 | 0.01201 | 0.02782 | 12.18x | 6.42 | 3.39 | 3.8e-16 | negative=562 delayed=147 null=0 est_entries=14879 |
| MUMPS_sym_metis_plain | METIS (sym=2) | 0.004234 | 0.01666 | 0.00217 | 0.01218 | 0.02883 | 12.63x | 6.78 | 3.65 | 2.8e-16 | negative=562 delayed=88 null=0 est_entries=16576 |

Baseline wide-solve share of ordering+factor+wide: 31%. Fastest factor+wide: UMFPACK_best (0.97x baseline). Least ordering-only fill: UMFPACK_sym_amd (1.90).

### ieee123C_1ph, exact Hessian: n = 1,353, nnz(K) = 8,072, RHS = 52 columns

| config | ordering used | ordering s | factor s | 1-RHS solve s | wide solve s | factor+wide s | vs baseline | fill (stored) | fill (ordering only) | wide residual | pivots |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| UMFPACK_default | AMD/COLAMD (UNSYMMETRIC) | 0.001054 | 0.003747 | 5.88e-05 | 0.001445 | 0.005193 | 1.00x | 2.77 | 2.38 | 5.5e-11 | offdiag_pivots=-1 |
| UMFPACK_unsym_colamd | AMD/COLAMD (UNSYMMETRIC) | 0.001024 | 0.003818 | 5.88e-05 | 0.00146 | 0.005278 | 1.02x | 2.77 | 2.38 | 5.5e-11 | offdiag_pivots=-1 |
| UMFPACK_sym_amd | AMD/COLAMD (SYMMETRIC) | 0.0006937 | 0.002648 | 5.71e-05 | 0.001463 | 0.004111 | 0.79x | 2.95 | 1.97 | 1.6e-13 | offdiag_pivots=237 |
| UMFPACK_unsym_metis | METIS (UNSYMMETRIC) | 0.005418 | 0.003122 | 6.15e-05 | 0.001457 | 0.004579 | 0.88x | 2.67 | 2.41 | 4.2e-12 | offdiag_pivots=-1 |
| UMFPACK_sym_metis | METIS (SYMMETRIC) | 0.003348 | 0.003535 | 6.82e-05 | 0.00164 | 0.005175 | 1.00x | 3.31 | 2.39 | 3.8e-13 | offdiag_pivots=249 |
| UMFPACK_cholmod | AMD/COLAMD (UNSYMMETRIC) | 0.00185 | 0.003848 | 5.92e-05 | 0.001439 | 0.005286 | 1.02x | 2.98 | 2.50 | 8.3e-12 | offdiag_pivots=-1 |
| UMFPACK_best | AMD/COLAMD (UNSYMMETRIC) | 0.0127 | 0.004025 | 5.81e-05 | 0.001458 | 0.005483 | 1.06x | 2.76 | 2.37 | 1.9e-12 | offdiag_pivots=-1 |
| MUMPS_sym_auto | AMF (sym=2) | 0.001564 | 0.0139 | 0.0019 | 0.01088 | 0.02478 | 4.77x | 6.61 | 3.83 | 6.5e-16 | negative=562 delayed=0 null=0 est_entries=26685 |
| MUMPS_sym_amd | AMD (sym=2) | 0.001261 | 0.01502 | 0.0019 | 0.01113 | 0.02615 | 5.04x | 7.48 | 4.36 | 7.0e-16 | negative=562 delayed=33 null=0 est_entries=29408 |
| MUMPS_sym_amf | AMF (sym=2) | 0.001547 | 0.01386 | 0.00191 | 0.01068 | 0.02453 | 4.72x | 6.61 | 3.83 | 6.5e-16 | negative=562 delayed=0 null=0 est_entries=26685 |
| MUMPS_sym_qamd | QAMD (sym=2) | 0.00129 | 0.01423 | 0.00188 | 0.01055 | 0.02478 | 4.77x | 7.41 | 4.14 | 7.1e-16 | negative=562 delayed=21 null=0 est_entries=29448 |
| MUMPS_sym_pord | PORD (sym=2) | 0.003255 | 0.01363 | 0.00199 | 0.01111 | 0.02474 | 4.76x | 7.70 | 4.72 | 6.4e-16 | negative=562 delayed=17 null=0 est_entries=30623 |
| MUMPS_sym_scotch | SCOTCH (sym=2) | 0.003383 | 0.01467 | 0.0018 | 0.01016 | 0.02484 | 4.78x | 9.05 | 4.73 | 6.5e-16 | negative=562 delayed=14 null=0 est_entries=36102 |
| MUMPS_sym_metis | METIS (sym=2) | 0.003751 | 0.01427 | 0.00191 | 0.01097 | 0.02523 | 4.86x | 7.45 | 4.42 | 6.5e-16 | negative=562 delayed=23 null=0 est_entries=29531 |
| MUMPS_sym_metis_compressed | METIS (sym=2) | 0.00381 | 0.01431 | 0.00194 | 0.01088 | 0.02519 | 4.85x | 7.45 | 4.42 | 6.5e-16 | negative=562 delayed=23 null=0 est_entries=29531 |
| MUMPS_sym_ipopt | AMF (sym=2) | 0.001558 | 0.01312 | 0.00193 | 0.01089 | 0.02401 | 4.62x | 6.61 | 3.83 | 4.6e-13 | negative=562 delayed=0 null=0 est_entries=26685 |
| MUMPS_sym_amd_plain | AMD (sym=2) | 0.0005549 | 0.02201 | 0.00229 | 0.01256 | 0.03457 | 6.66x | 6.44 | 2.89 | 7.0e-16 | negative=562 delayed=295 null=0 est_entries=18350 |
| MUMPS_sym_amf_plain | AMF (sym=2) | 0.000618 | 0.02341 | 0.00236 | 0.01279 | 0.0362 | 6.97x | 6.84 | 3.12 | 8.1e-16 | negative=562 delayed=341 null=0 est_entries=19116 |
| MUMPS_sym_metis_plain | METIS (sym=2) | 0.005136 | 0.01496 | 0.002 | 0.01098 | 0.02594 | 5.00x | 6.69 | 4.06 | 7.6e-16 | negative=562 delayed=83 null=0 est_entries=25296 |

Baseline wide-solve share of ordering+factor+wide: 23%. Fastest factor+wide: UMFPACK_sym_amd (0.79x baseline). Least ordering-only fill: UMFPACK_sym_amd (1.97).

### ieee2522C_1ph, diag Hessian: n = 23,691, nnz(K) = 96,025, RHS = 250 columns

| config | ordering used | ordering s | factor s | 1-RHS solve s | wide solve s | factor+wide s | vs baseline | fill (stored) | fill (ordering only) | wide residual | pivots |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| UMFPACK_default | AMD/COLAMD (UNSYMMETRIC) | 0.01094 | 0.02333 | 0.000491 | 0.1125 | 0.1358 | 1.00x | 2.24 | 2.06 | 5.1e-15 | offdiag_pivots=-1 |
| UMFPACK_unsym_colamd | AMD/COLAMD (UNSYMMETRIC) | 0.01044 | 0.02289 | 0.000469 | 0.1221 | 0.145 | 1.07x | 2.24 | 2.06 | 5.1e-15 | offdiag_pivots=-1 |
| UMFPACK_sym_amd | AMD/COLAMD (SYMMETRIC) | 0.00845 | 0.0327 | 0.000687 | 0.1405 | 0.1732 | 1.28x | 2.89 | 1.97 | 2.7e-15 | offdiag_pivots=5162 |
| UMFPACK_unsym_metis | METIS (UNSYMMETRIC) | 0.1153 | 0.03464 | 0.000691 | 0.1663 | 0.2009 | 1.48x | 2.77 | 2.74 | 1.3e-14 | offdiag_pivots=-1 |
| UMFPACK_sym_metis | METIS (SYMMETRIC) | 0.09349 | 0.03776 | 0.000686 | 0.1703 | 0.2081 | 1.53x | 3.07 | 2.53 | 8.0e-15 | offdiag_pivots=4640 |
| UMFPACK_cholmod | AMD/COLAMD (UNSYMMETRIC) | 0.01866 | 0.03258 | 0.000455 | 0.1242 | 0.1568 | 1.15x | 2.19 | 2.03 | 8.8e-15 | offdiag_pivots=-1 |
| UMFPACK_best | AMD/COLAMD (UNSYMMETRIC) | 0.2379 | 0.02335 | 0.000469 | 0.1047 | 0.128 | 0.94x | 2.23 | 2.06 | 4.2e-15 | offdiag_pivots=-1 |
| MUMPS_sym_auto | AMF (sym=2) | 0.04795 | 0.2964 | 0.0343 | 1.115 | 1.411 | 10.39x | 5.73 | 3.19 | 1.5e-15 | negative=10336 delayed=304 null=0 est_entries=270685 |
| MUMPS_sym_amd | AMD (sym=2) | 0.01814 | 0.3658 | 0.0336 | 1.127 | 1.493 | 10.99x | 6.30 | 3.66 | 2.0e-15 | negative=10336 delayed=1739 null=0 est_entries=277207 |
| MUMPS_sym_amf | AMF (sym=2) | 0.04657 | 0.3039 | 0.036 | 1.166 | 1.469 | 10.82x | 5.73 | 3.19 | 1.5e-15 | negative=10336 delayed=304 null=0 est_entries=270685 |
| MUMPS_sym_qamd | QAMD (sym=2) | 0.01827 | 0.3431 | 0.0332 | 1.103 | 1.446 | 10.65x | 6.30 | 3.66 | 2.0e-15 | negative=10336 delayed=1739 null=0 est_entries=277207 |
| MUMPS_sym_pord | PORD (sym=2) | 0.07825 | 0.2811 | 0.035 | 1.159 | 1.441 | 10.61x | 6.99 | 4.06 | 2.0e-15 | negative=10336 delayed=629 null=0 est_entries=325618 |
| MUMPS_sym_scotch | SCOTCH (sym=2) | 0.05645 | 0.3372 | 0.0305 | 1.034 | 1.371 | 10.10x | 8.01 | 4.09 | 8.8e-16 | negative=10336 delayed=543 null=0 est_entries=374163 |
| MUMPS_sym_metis | METIS (sym=2) | 0.0998 | 0.3098 | 0.0346 | 1.102 | 1.412 | 10.40x | 7.43 | 4.31 | 2.2e-15 | negative=10336 delayed=545 null=0 est_entries=347568 |
| MUMPS_sym_metis_compressed | METIS (sym=2) | 0.09919 | 0.3041 | 0.0343 | 1.141 | 1.445 | 10.64x | 7.43 | 4.31 | 2.2e-15 | negative=10336 delayed=545 null=0 est_entries=347568 |
| MUMPS_sym_ipopt | AMF (sym=2) | 0.0472 | 0.2756 | 0.0348 | 1.111 | 1.387 | 10.21x | 5.64 | 3.19 | 3.3e-11 | negative=10336 delayed=4 null=0 est_entries=270685 |
| MUMPS_sym_amd_plain | AMD (sym=2) | 0.005928 | 0.3161 | 0.0353 | 1.063 | 1.379 | 10.15x | 5.95 | 3.22 | 1.5e-15 | negative=10336 delayed=2745 null=0 est_entries=249975 |
| MUMPS_sym_amf_plain | AMF (sym=2) | 0.007259 | 0.3195 | 0.0375 | 1.104 | 1.423 | 10.48x | 6.28 | 3.36 | 3.3e-15 | negative=10336 delayed=3651 null=0 est_entries=254839 |
| MUMPS_sym_metis_plain | METIS (sym=2) | 0.1222 | 0.3332 | 0.0356 | 1.119 | 1.452 | 10.69x | 7.09 | 4.06 | 1.4e-15 | negative=10336 delayed=1769 null=0 est_entries=312739 |

Baseline wide-solve share of ordering+factor+wide: 77%. Fastest factor+wide: UMFPACK_best (0.94x baseline). Least ordering-only fill: UMFPACK_sym_amd (1.97).

### ieee2522C_1ph, exact Hessian: n = 23,691, nnz(K) = 162,818, RHS = 250 columns

| config | ordering used | ordering s | factor s | 1-RHS solve s | wide solve s | factor+wide s | vs baseline | fill (stored) | fill (ordering only) | wide residual | pivots |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| UMFPACK_default | AMD/COLAMD (UNSYMMETRIC) | 0.03588 | 0.06561 | 0.000926 | 0.1594 | 0.225 | 1.00x | 2.65 | 2.41 | 5.8e-12 | offdiag_pivots=-1 |
| UMFPACK_unsym_colamd | AMD/COLAMD (UNSYMMETRIC) | 0.03457 | 0.067 | 0.000869 | 0.1553 | 0.2223 | 0.99x | 2.65 | 2.41 | 5.8e-12 | offdiag_pivots=-1 |
| UMFPACK_sym_amd | AMD/COLAMD (SYMMETRIC) | 0.01236 | 0.09174 | 0.001 | 0.1787 | 0.2704 | 1.20x | 3.16 | 1.81 | 3.6e-14 | offdiag_pivots=5090 |
| UMFPACK_unsym_metis | METIS (UNSYMMETRIC) | 0.221 | 0.05775 | 0.000903 | 0.1794 | 0.2372 | 1.05x | 2.35 | 2.25 | 2.4e-14 | offdiag_pivots=-1 |
| UMFPACK_sym_metis | METIS (SYMMETRIC) | 0.09357 | 0.09911 | 0.00111 | 0.1975 | 0.2966 | 1.32x | 3.24 | 2.22 | 3.0e-14 | offdiag_pivots=5235 |
| UMFPACK_cholmod | AMD/COLAMD (UNSYMMETRIC) | 0.1034 | 0.06639 | 0.000903 | 0.1927 | 0.259 | 1.15x | 2.73 | 2.37 | 7.4e-12 | offdiag_pivots=-1 |
| UMFPACK_best | AMD/COLAMD (UNSYMMETRIC) | 0.4782 | 0.07917 | 0.000881 | 0.1826 | 0.2617 | 1.16x | 2.48 | 2.22 | 2.6e-13 | offdiag_pivots=-1 |
| MUMPS_sym_auto | AMF (sym=2) | 0.04741 | 0.301 | 0.0358 | 1.188 | 1.489 | 6.62x | 4.63 | 2.89 | 1.2e-14 | negative=10336 delayed=160 null=0 est_entries=374262 |
| MUMPS_sym_amd | AMD (sym=2) | 0.02331 | 0.3232 | 0.035 | 1.174 | 1.497 | 6.65x | 4.99 | 3.19 | 1.2e-14 | negative=10336 delayed=1049 null=0 est_entries=387768 |
| MUMPS_sym_amf | AMF (sym=2) | 0.0481 | 0.2954 | 0.0355 | 1.184 | 1.479 | 6.57x | 4.63 | 2.89 | 1.2e-14 | negative=10336 delayed=160 null=0 est_entries=374262 |
| MUMPS_sym_qamd | QAMD (sym=2) | 0.02556 | 0.2877 | 0.0354 | 1.096 | 1.384 | 6.15x | 5.36 | 3.36 | 1.2e-14 | negative=10336 delayed=544 null=0 est_entries=424546 |
| MUMPS_sym_pord | PORD (sym=2) | 0.06938 | 0.3143 | 0.0404 | 1.154 | 1.468 | 6.52x | 5.20 | 3.36 | 1.2e-14 | negative=10336 delayed=706 null=0 est_entries=410691 |
| MUMPS_sym_scotch | SCOTCH (sym=2) | 0.08529 | 0.3197 | 0.0344 | 1.142 | 1.462 | 6.50x | 6.02 | 3.24 | 1.4e-14 | negative=10336 delayed=391 null=0 est_entries=481099 |
| MUMPS_sym_metis | METIS (sym=2) | 0.09684 | 0.3018 | 0.0353 | 1.124 | 1.426 | 6.33x | 5.52 | 3.51 | 1.2e-14 | negative=10336 delayed=715 null=0 est_entries=436790 |
| MUMPS_sym_metis_compressed | METIS (sym=2) | 0.1182 | 0.3254 | 0.0391 | 1.222 | 1.548 | 6.88x | 5.52 | 3.51 | 1.2e-14 | negative=10336 delayed=715 null=0 est_entries=436790 |
| MUMPS_sym_ipopt | AMF (sym=2) | 0.04742 | 0.2659 | 0.0358 | 1.2 | 1.466 | 6.51x | 4.60 | 2.89 | 7.1e-12 | negative=10336 delayed=11 null=0 est_entries=374262 |
| MUMPS_sym_amd_plain | AMD (sym=2) | 0.008215 | 0.3156 | 0.0399 | 1.018 | 1.333 | 5.92x | 4.60 | 2.64 | 1.6e-14 | negative=10336 delayed=2782 null=0 est_entries=328570 |
| MUMPS_sym_amf_plain | AMF (sym=2) | 0.009232 | 0.3159 | 0.0356 | 1.038 | 1.354 | 6.02x | 4.82 | 2.79 | 1.4e-14 | negative=10336 delayed=3577 null=0 est_entries=334072 |
| MUMPS_sym_metis_plain | METIS (sym=2) | 0.1314 | 0.3252 | 0.0458 | 1.148 | 1.473 | 6.54x | 5.35 | 3.17 | 1.3e-14 | negative=10336 delayed=1743 null=0 est_entries=397801 |

Baseline wide-solve share of ordering+factor+wide: 61%. Fastest factor+wide: UMFPACK_unsym_colamd (0.99x baseline). Least ordering-only fill: UMFPACK_sym_amd (1.81).

### large10kC_1ph, diag Hessian: n = 96,968, nnz(K) = 393,075, RHS = 1021 columns

| config | ordering used | ordering s | factor s | 1-RHS solve s | wide solve s | factor+wide s | vs baseline | fill (stored) | fill (ordering only) | wide residual | pivots |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| UMFPACK_default | AMD/COLAMD (UNSYMMETRIC) | 0.05339 | 0.1536 | 0.00253 | 2.159 | 2.312 | 1.00x | 2.22 | 2.45 | 3.7e-16 | offdiag_pivots=-1 |
| UMFPACK_unsym_colamd | AMD/COLAMD (UNSYMMETRIC) | 0.05113 | 0.148 | 0.00236 | 2.148 | 2.296 | 0.99x | 2.22 | 2.45 | 3.7e-16 | offdiag_pivots=-1 |
| UMFPACK_sym_amd | AMD/COLAMD (SYMMETRIC) | 0.05852 | 0.2509 | 0.0112 | 5.603 | 5.854 | 2.53x | 6.10 | 2.29 | 2.0e-12 | offdiag_pivots=19609 |
| UMFPACK_unsym_metis | METIS (UNSYMMETRIC) | 0.5705 | 0.1735 | 0.00327 | 2.965 | 3.139 | 1.36x | 2.89 | 8.08 | 5.9e-15 | offdiag_pivots=-1 |
| UMFPACK_sym_metis | METIS (SYMMETRIC) | 0.3864 | 0.2302 | 0.00379 | 3.61 | 3.84 | 1.66x | 3.33 | 4.02 | 7.9e-12 | offdiag_pivots=18774 |
| UMFPACK_cholmod | AMD/COLAMD (UNSYMMETRIC) | 0.08319 | 0.1408 | 0.00223 | 2.221 | 2.361 | 1.02x | 2.19 | 2.33 | 5.9e-16 | offdiag_pivots=-1 |
| UMFPACK_best | AMD/COLAMD (UNSYMMETRIC) | 1.125 | 0.1136 | 0.00231 | 2.185 | 2.299 | 0.99x | 2.17 | 2.45 | 3.5e-16 | offdiag_pivots=-1 |
| MUMPS_sym_auto | METIS (sym=2) | 0.4931 | 1.161 | 0.15 | 18.34 | 19.5 | 8.43x | 6.91 | 4.12 | 7.7e-15 | negative=42303 delayed=6577 null=0 est_entries=1266604 |
| MUMPS_sym_amd | AMD (sym=2) | 0.05897 | 1.151 | 0.167 | 18.93 | 20.08 | 8.69x | 5.74 | 3.33 | 8.5e-15 | negative=42303 delayed=10023 null=0 est_entries=1003721 |
| MUMPS_sym_amf | AMF (sym=2) | 0.06707 | 1.155 | 0.18 | 20.4 | 21.55 | 9.32x | 5.82 | 3.36 | 8.4e-15 | negative=42303 delayed=10517 null=0 est_entries=1012496 |
| MUMPS_sym_qamd | QAMD (sym=2) | 0.07779 | 1.2 | 0.157 | 20.95 | 22.15 | 9.58x | 5.94 | 3.38 | 1.2e-14 | negative=42303 delayed=9453 null=0 est_entries=1051610 |
| MUMPS_sym_pord | PORD (sym=2) | 0.3257 | 1.206 | 0.173 | 20.15 | 21.36 | 9.24x | 6.23 | 3.46 | 1.1e-14 | negative=42303 delayed=9882 null=0 est_entries=1097733 |
| MUMPS_sym_scotch | SCOTCH (sym=2) | 0.6036 | 1.483 | 0.141 | 17.2 | 18.69 | 8.08x | 7.74 | 3.38 | 2.3e-15 | negative=42303 delayed=7309 null=0 est_entries=1397525 |
| MUMPS_sym_metis | METIS (sym=2) | 0.5677 | 1.248 | 0.208 | 20.04 | 21.29 | 9.21x | 6.91 | 4.12 | 7.7e-15 | negative=42303 delayed=6577 null=0 est_entries=1266604 |
| MUMPS_sym_metis_compressed | METIS (sym=2) | 0.3109 | 1.183 | 0.131 | 18.06 | 19.24 | 8.32x | 7.54 | 4.40 | 6.8e-15 | negative=42303 delayed=2028 null=0 est_entries=1446605 |
| MUMPS_sym_amd_plain | AMD (sym=2) | 0.02483 | 1.183 | 0.161 | 19.88 | 21.07 | 9.11x | 5.74 | 3.33 | 9.2e-15 | negative=42303 delayed=10021 null=0 est_entries=1003721 |
| MUMPS_sym_amf_plain | AMF (sym=2) | 0.02842 | 1.261 | 0.194 | 20.09 | 21.35 | 9.23x | 5.82 | 3.36 | 9.2e-15 | negative=42303 delayed=10513 null=0 est_entries=1012496 |
| MUMPS_sym_metis_plain | METIS (sym=2) | 0.5493 | 1.282 | 0.168 | 19.34 | 20.63 | 8.92x | 6.91 | 4.12 | 1.4e-14 | negative=42303 delayed=6574 null=0 est_entries=1266604 |
| MUMPS_sym_ipopt | METIS (sym=2) | 0.5605 | 1.203 | 0.197 | 19.12 | 20.32 | 8.79x | 6.80 | 4.12 | 1.6e-10 | negative=42303 delayed=5118 null=0 est_entries=1266604 |

Baseline wide-solve share of ordering+factor+wide: 91%. Fastest factor+wide: UMFPACK_unsym_colamd (0.99x baseline). Least ordering-only fill: UMFPACK_sym_amd (2.29).

### large10kC_1ph, exact Hessian: n = 96,968, nnz(K) = 1,453,094, RHS = 1021 columns

| config | ordering used | ordering s | factor s | 1-RHS solve s | wide solve s | factor+wide s | vs baseline | fill (stored) | fill (ordering only) | wide residual | pivots |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| UMFPACK_default | AMD/COLAMD (UNSYMMETRIC) | 0.148 | 0.3096 | 0.00393 | 3.594 | 3.903 | 1.00x | 1.45 | 2.39 | 8.0e-16 | offdiag_pivots=-1 |
| UMFPACK_unsym_colamd | AMD/COLAMD (UNSYMMETRIC) | 0.1235 | 0.3314 | 0.00396 | 3.771 | 4.102 | 1.05x | 1.45 | 2.39 | 8.0e-16 | offdiag_pivots=-1 |
| UMFPACK_sym_amd | AMD/COLAMD (SYMMETRIC) | 0.1434 | 0.4801 | 0.00611 | 5.253 | 5.733 | 1.47x | 2.00 | 4.87 | 7.9e-12 | offdiag_pivots=20549 |
| UMFPACK_unsym_metis | METIS (UNSYMMETRIC) | 4.098 | 0.6471 | 0.00956 | 9.883 | 10.53 | 2.70x | 4.52 | 10.15 | 9.8e-12 | offdiag_pivots=-1 |
| UMFPACK_sym_metis | METIS (SYMMETRIC) | 0.6575 | 0.5065 | 0.00602 | 6.032 | 6.539 | 1.68x | 2.23 | 4.44 | 2.7e-12 | offdiag_pivots=20833 |
| UMFPACK_cholmod | AMD/COLAMD (UNSYMMETRIC) | 4.607 | 0.3442 | 0.00429 | 4.33 | 4.674 | 1.20x | 1.67 | 8.50 | 3.3e-11 | offdiag_pivots=-1 |
| UMFPACK_best | AMD/COLAMD (UNSYMMETRIC) | 11.6 | 0.3843 | 0.00464 | 4.48 | 4.864 | 1.25x | 1.78 | 8.72 | 1.4e-11 | offdiag_pivots=-1 |
| MUMPS_sym_auto | METIS (sym=2) | 0.7678 | 1.182 | 0.186 | 19.77 | 20.95 | 5.37x | 2.73 | 1.91 | 7.6e-15 | negative=42303 delayed=6218 null=0 est_entries=1877962 |
| MUMPS_sym_amd | AMD (sym=2) | 0.1287 | 1.185 | 0.167 | 17.27 | 18.46 | 4.73x | 2.55 | 1.74 | 2.2e-14 | negative=42303 delayed=9635 null=0 est_entries=1719745 |
| MUMPS_sym_amf | AMF (sym=2) | 0.1232 | 1.14 | 0.173 | 18.18 | 19.32 | 4.95x | 2.54 | 1.81 | 2.2e-14 | negative=42303 delayed=10759 null=0 est_entries=1702765 |
| MUMPS_sym_qamd | QAMD (sym=2) | 0.114 | 1.281 | 0.124 | 16.4 | 17.68 | 4.53x | 3.14 | 1.84 | 1.1e-14 | negative=42303 delayed=8472 null=0 est_entries=2125217 |
| MUMPS_sym_pord | PORD (sym=2) | 0.5781 | 1.125 | 0.164 | 19.45 | 20.57 | 5.27x | 2.72 | 1.88 | 8.7e-15 | negative=42303 delayed=9313 null=0 est_entries=1848907 |
| MUMPS_sym_scotch | SCOTCH (sym=2) | 0.7236 | 1.745 | 0.135 | 17.47 | 19.21 | 4.92x | 4.73 | 2.58 | 3.3e-15 | negative=42303 delayed=6982 null=0 est_entries=2838892 |
| MUMPS_sym_metis | METIS (sym=2) | 0.7963 | 1.137 | 0.173 | 19.6 | 20.74 | 5.31x | 2.73 | 1.91 | 7.6e-15 | negative=42303 delayed=6218 null=0 est_entries=1877962 |
| MUMPS_sym_metis_compressed | METIS (sym=2) | 0.6801 | 1.132 | 0.17 | 18.02 | 19.15 | 4.91x | 3.07 | 2.17 | 7.2e-15 | negative=42303 delayed=2194 null=0 est_entries=2191563 |
| MUMPS_sym_amd_plain | AMD (sym=2) | 0.05022 | 1.206 | 0.165 | 18.12 | 19.32 | 4.95x | 2.55 | 1.74 | 9.3e-15 | negative=42303 delayed=9635 null=0 est_entries=1719745 |
| MUMPS_sym_amf_plain | AMF (sym=2) | 0.04934 | 1.176 | 0.151 | 18.58 | 19.76 | 5.06x | 2.54 | 1.81 | 1.4e-14 | negative=42303 delayed=10759 null=0 est_entries=1702765 |
| MUMPS_sym_metis_plain | METIS (sym=2) | 0.7414 | 1.183 | 0.162 | 19.47 | 20.65 | 5.29x | 2.73 | 1.91 | 1.4e-14 | negative=42303 delayed=6219 null=0 est_entries=1877962 |
| MUMPS_sym_ipopt | METIS (sym=2) | 0.7855 | 1.068 | 0.159 | 19.7 | 20.77 | 5.32x | 2.71 | 1.91 | 4.7e-11 | negative=42303 delayed=5704 null=0 est_entries=1877962 |

Baseline wide-solve share of ordering+factor+wide: 89%. Fastest factor+wide: UMFPACK_default (1.00x baseline). Least ordering-only fill: MUMPS_sym_amd (1.74).
