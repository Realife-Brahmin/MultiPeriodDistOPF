/* Thin C shim over HSL_MA97 (double) for calling from Julia: one factorization
 * held in static state, plain-argument entry points, 1-based (Fortran) arrays.
 * Not part of HSL; links against the locally built hsl_ma97 objects. */
#include "hsl_ma97d.h"

static void *akeep = 0, *fkeep = 0;
static struct ma97_control_d ctl;
static struct ma97_info_d inf;

void s97_init(int ordering, int scaling, int solve_blas3, double u) {
    if (akeep || fkeep) ma97_finalise_d(&akeep, &fkeep);
    akeep = 0; fkeep = 0;
    ma97_default_control_d(&ctl);
    ctl.f_arrays = 1;
    ctl.ordering = ordering;
    ctl.scaling = scaling;
    ctl.solve_blas3 = solve_blas3;
    ctl.print_level = -1;
    if (u > 0) ctl.u = u;
}

int s97_analyse(int n, const int *ptr, const int *row) {
    ma97_analyse_d(0, n, ptr, row, 0, &akeep, &ctl, &inf, 0);
    return inf.flag;
}

int s97_factor(const int *ptr, const int *row, const double *val) {
    /* 4 = HSL_MATRIX_REAL_SYM_INDEF */
    ma97_factor_d(4, ptr, row, val, &akeep, &fkeep, &ctl, &inf, 0);
    return inf.flag;
}

int s97_solve(int nrhs, double *x, int ldx) {
    ma97_solve_d(0, nrhs, x, ldx, &akeep, &fkeep, &ctl, &inf);
    return inf.flag;
}

void s97_info(double *out) {
    out[0] = inf.flag;       out[1] = (double)inf.num_factor; out[2] = (double)inf.num_flops;
    out[3] = inf.num_neg;    out[4] = inf.num_delay;          out[5] = inf.num_two;
    out[6] = inf.maxfront;   out[7] = inf.ordering;           out[8] = inf.matrix_rank;
    out[9] = inf.maxdepth;   out[10] = inf.num_sup;           out[11] = inf.stat;
}

void s97_free(void) {
    if (akeep || fkeep) ma97_finalise_d(&akeep, &fkeep);
    akeep = 0; fkeep = 0;
}
