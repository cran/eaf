#include <R.h>
#include <Rinternals.h>
#include <R_ext/Error.h>
#include <R_ext/Memory.h>
#include "eaf.h"

#ifndef DEBUG
#define DEBUG 0
#endif

#if DEBUG >= 1
#define DEBUG1(X) X;(void)0
#else  
#define DEBUG1(X) (void)0
#endif

#define CHECK_ARG_IS_INT_VECTOR(A)					\
    if (!isInteger(A) || !isVector(A))					\
	error("Argument '" #A "' is not an integer vector");

/*
 * Unpack an integer vector stored in SEXP S.
 */
#define SEXP_2_INT_VECTOR(S, I, N)               \
    CHECK_ARG_IS_INT_VECTOR(S);                  \
    int *I = INTEGER(S);                         \
    const R_len_t N = length(S);

#define SEXP_2_INT(S,var)                                               \
    int var = asInteger(S);                                             \
    if (var == NA_INTEGER)                                              \
        error ("Argument '" #S "' is not an integer");

#define SEXP_2_STRING(S,var)                                            \
    if (!isString(S) || length(S) != 1)                                 \
        error ("Argument '" #S "' is not a string");                    \
    const char * var = CHAR(STRING_ELT(S,0));

SEXP compute_eaf_C(SEXP DATA, SEXP NOBJ, SEXP CUMSIZES, SEXP NRUNS, SEXP PERCENTILES);

SEXP compute_eaf_C(SEXP DATA, SEXP NOBJ, SEXP CUMSIZES, SEXP NRUNS, SEXP PERCENTILE)
{
    eaf_t **eaf;
    int k;
    SEXP_2_INT(NOBJ, nobj);
    SEXP_2_INT_VECTOR(CUMSIZES, cumsizes, cumsizes_len);
    SEXP_2_INT(NRUNS, nruns);
    SEXP_2_INT_VECTOR(PERCENTILE, percentile, nlevels);

    if (cumsizes_len < nruns)
        error("length of cumsizes (%d) is less than nruns (%d)",
              cumsizes_len, nruns);

    int *level = malloc(sizeof(int) * nlevels);
    for (k = 0; k < nlevels; k++)
        level[k] = percentile2level(percentile[k], nruns);

    double *data = REAL(DATA);

    DEBUG1(
        Rprintf ("attsurf ({(%f,%f)...}, %d, { %d",
                 data[0], data[1], nobj, cumsizes[0]);
        for (k = 1; k < nruns; k++) {
            Rprintf (", %d", cumsizes[k]);
        }
        Rprintf ("}, %d, { %d", nruns, level[0]);
        for (k = 1; k < nlevels; k++) {
            Rprintf (", %d", level[k]);
        }
        Rprintf ("}, %d)\n", nlevels);
        );

    eaf = attsurf (data, nobj, cumsizes, nruns, level, nlevels);

    DEBUG1(
        Rprintf ("eaf computed\n");
        for (k = 0; k < nlevels; k++) {
            Rprintf ("eaf[%d] = %d\n", k, eaf[k]->size);
        });

    int totalpoints = 0;
    for (k = 0; k < nlevels; k++) {
        totalpoints += eaf[k]->size;
    }

    SEXP mat;
    double *rmat;
    PROTECT(mat = allocMatrix(REALSXP, totalpoints, nobj + 1));
    rmat = REAL(mat);

    int pos = 0;
    for (k = 0; k < nlevels; k++) {
        int npoints = eaf[k]->size;
        int i;

        DEBUG1(
            int totalsize = npoints * nobj;
            Rprintf ("totalpoints eaf[%d] = %d\n", k, totalsize)
            );

        for (i = 0; i < npoints; i++) {
            rmat[pos] = eaf[k]->data[i * nobj];
            rmat[pos + totalpoints] = eaf[k]->data[1 + i * nobj];
            rmat[pos + 2 * totalpoints] = percentile[k];
            pos++;
        }
        eaf_delete (eaf[k]);
    }
    free(eaf);
    UNPROTECT (1);
    return mat;
}

SEXP compute_eafdiff_C(SEXP DATA, SEXP NOBJ, SEXP CUMSIZES, SEXP NRUNS, SEXP DIVISION, SEXP INTERVALS);

SEXP compute_eafdiff_C(SEXP DATA, SEXP NOBJ, SEXP CUMSIZES, SEXP NRUNS, SEXP DIVISION, SEXP INTERVALS)
{
    eaf_t **eaf;
    int k;
    SEXP_2_INT(NOBJ, nobj);
    SEXP_2_INT_VECTOR(CUMSIZES, cumsizes, cumsizes_len);
    SEXP_2_INT(NRUNS, nruns);
    SEXP_2_INT(INTERVALS, intervals);

    if (cumsizes_len < nruns)
        error("length of cumsizes (%d) is less than nruns (%d)",
              cumsizes_len, nruns);

    int nsets1 = nruns / 2;
    int nsets2 = nruns - nsets1;

    int *level = malloc(sizeof(int) * nruns);
    for (k = 0; k < nruns; k++)
        level[k] = k + 1;

    double *data = REAL(DATA);

    DEBUG1 (
        Rprintf ("attsurf ({(%f,%f)...}, %d, { %d",
                 data[0], data[1], nobj, cumsizes[0]);
        for (k = 1; k < nruns; k++) {
            Rprintf (", %d", cumsizes[k]);
        }
        Rprintf ("}, %d, ALL, %d)\n", nruns, nruns);
        );

    eaf = attsurf (data, nobj, cumsizes, nruns, level, nruns);

    DEBUG1(
        Rprintf ("eaf computed\n");
        for (k = 0; k < nruns; k++) {
            Rprintf ("eaf[%d] = %d\n", k, eaf[k]->size);
        });

    int totalpoints = 0;
    for (k = 0; k < nruns; k++) {
        totalpoints += eaf[k]->size;
    }
    SEXP mat;
    double *rmat;
    PROTECT(mat = allocMatrix(REALSXP, totalpoints, nobj + 1));
    rmat = REAL(mat);

    int pos = 0;
    for (k = 0; k < nruns; k++) {
        int npoints = eaf[k]->size;
        int i;

        DEBUG1(
            int totalsize = npoints * nobj;
            Rprintf ("totalpoints eaf[%d] = %d\n", k, totalsize)
            );

        for (i = 0; i < npoints; i++) {
            rmat[pos] = eaf[k]->data[i * nobj];
            rmat[pos + totalpoints] = eaf[k]->data[1 + i * nobj];
            pos++;
        }
    }

    pos += totalpoints;
    for (k = 0; k < nruns; k++) {
        int i;
        int npoints = eaf[k]->size;
        for (i = 0; i < npoints; i++) {
            int count_left;
            int count_right;
            attained_left_right (eaf[k]->attained + i * eaf[k]->nruns,
                                 nruns/2, nruns, &count_left, &count_right);
            rmat[pos] = intervals * (double) ((count_left / (double) nsets1) - 
                                              (count_right / (double) nsets2));
            pos++;
        }
        eaf_delete (eaf[k]);
    }
    free(eaf);
    UNPROTECT (1);
    return mat;
}

SEXP read_data_sets (SEXP FILENAME);

SEXP
read_data_sets (SEXP FILENAME)
{
    SEXP_2_STRING(FILENAME, filename);
    objective_t *data = NULL;
    int* cumsizes = NULL;
    int nobj = 0, nruns = 0;

    /* Rprintf ("filename: %s\n", filename); */

    read_objective_t_data (filename, &data, &nobj, &cumsizes, &nruns);

    const int ntotal = cumsizes[nruns - 1];
    int * runtab = malloc (ntotal * sizeof(int));
    int k, j, i;
    for (k = 0, j = 0; k < ntotal; k++) {
        if (k == cumsizes[j])
            j++;
        runtab[k] = j + 1;
    }

    SEXP DATA;
    PROTECT(DATA = allocMatrix(REALSXP, cumsizes[nruns-1], nobj + 1));
    double *rdata = REAL(DATA);
    int pos = 0;
    for (j = 0; j < nobj; j++) {
        for (i = 0; i < ntotal; i++) {
            rdata[pos] = data[j + i * nobj];
            pos++;
        }
    }
    for (j = 0; j < ntotal; j++, pos++) {
        rdata[pos] = runtab[j];
    }
    free(data);
    free(cumsizes);
    UNPROTECT (1);
    return DATA;
}
