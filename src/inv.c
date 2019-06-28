/* Copyright 2019, Alistair Boyle, 3-clause BSD License */
#include "config.h"
#include <stdlib.h> /* malloc, free */
#include <assert.h> /* assert */
#include <stdio.h> /* printf */
#include <string.h> /* memset */
#include <float.h> /* DBL_MIN */
#include <lapacke.h>

#include "matrix.h"
#include "fwd.h"
#include "inv.h"

#ifdef UNIT_TESTING
extern void * _mock_test_malloc(const size_t size, const char * file, const int line);
extern void _test_free(void * const ptr, const char * file, const int line);
#define malloc(size) _mock_test_malloc(size, __FILE__, __LINE__)
#define free(ptr) _test_free(ptr, __FILE__, __LINE__)
#endif

#define SUCCESS 1
#define FAILURE 0

// #define DEBUG_SVD_MATRICES 1

#ifdef DEBUG_SVD_MATRICES
static void set_matrix(matrix * X, char symbol[], char name[], char units[], int m, int n, double * x)
{
    memset(X, 0, sizeof(matrix));
    X->m = m;
    X->n = n;
    X->type = DENSE;
    X->symmetric = false;
    X->scale = 1.0;
    X->symbol = &(symbol[0]);
    X->name = &(name[0]);
    X->units = &(units[0]);
    X->x.dense = x;
}
#endif

/*
static int svd(matrix * J, matrix * U, matrix * sv, matrix * Vt)
{
    assert(J != NULL);
    const int m = J->m;
    const int n = J->n;
    return m*n;
}
*/

/* max/min: watch out for double evaluation of arguments e.g. a++ == a+=2 */
#define min(a,b) (((a) < (b)) ? (a) : (b) )
#define max(a,b) (((a) > (b)) ? (a) : (b) )
/* dense matrix indexing */
#define colmaj(i,j,ld) ((i) + (j)*(ld)) /* column major indexing */

/* solve Δx=(JᵀJ+λ²I)⁻¹Jb for b=b₂-b₁
 * where J=UΣVᵀ (by the SVD)
 * and
 * Δx=VΣ⁺Uᵀb=VϕΣ⁻¹Uᵀb with Σ⁺=(Σ²+λ²I)⁻¹Σ or ϕₖ=σₖ²/(σₖ²+λ²)
 * for filter factor ϕ
 */
#define phi(s,hp) ((s*s)/((s*s)+(hp*hp)))
static int inv_solve_svd(model * const mdl, double const * b, double * x)
{
    int ret = FAILURE;
    assert(mdl != NULL);
    assert(b != NULL);
    assert(x != NULL);
    const int m = mdl->n_stimmeas;
    const int n = mdl->n_params[0];
    /* init J */
    // printf("J is %d x %d = %d --> %d B\n", m, n, m * n, m * n * (int) sizeof(double));
    double * Ux = NULL, * sx = NULL, * Vtx = NULL, * superb = NULL;
    double * Jx = malloc(sizeof(double) * m * n);
    if(Jx == NULL) {
        printf("error: %s: out of memory\n", "jacobian");
        return FAILURE;
    }
    if(!calc_jacobian(mdl, &Jx[0])) {
        printf("error: inv_solve/calc_jacobian failed\n");
        free(Jx);
        return FAILURE;
    }
#ifdef DEBUG_SVD_MATRICES
    matrix J;
    set_matrix(&J, "J", "Jacobian", "Vm/S", m, n, &Jx[0]);
    printf_matrix(&J);
#endif
    /* calculate SVD of J; J is destroyed in the process */
    /* SVD config */
    const int layout = LAPACK_COL_MAJOR;
    /* A= m cols of U --> u; S= min(m,n) cols U --> u; O=min(m,n) cols U --> a (overwrite); 'N'=no U */
    const char jobu = 'A';
    const char jobvt = 'A';
    const int ldJ = (layout == LAPACK_ROW_MAJOR) ? n : m;
    const int ldU = ((jobu == 'S') && (layout == LAPACK_ROW_MAJOR)) ? min(m, n) : m;
    const int ldVt = ((jobu == 'S') && (layout == LAPACK_COL_MAJOR)) ? min(m, n) : n;
    /* SVD outputs */
    const int Usz = ((jobu == 'S') && (layout == LAPACK_COL_MAJOR)) ? ldU * min(m, n) : ldU * m;
    const int ssz = min(m, n);
    const int Vtsz = ((jobu == 'S') && (layout == LAPACK_COL_MAJOR)) ? ldVt * min(m, n) : ldVt * n;
    const int sbsz = min(m, n) - 1;
    Ux = malloc(sizeof(double) * Usz);
    sx = malloc(sizeof(double) * ssz);
    Vtx = malloc(sizeof(double) * Vtsz);
    superb = malloc(sizeof(double) * sbsz);
    assert(Ux != NULL);
    assert(sx != NULL);
    assert(Vtx != NULL);
    assert(superb != NULL);
    int info = LAPACKE_dgesvd(layout, jobu, jobvt,
                              m, n, Jx, ldJ,
                              sx,
                              Ux, ldU,
                              Vtx, ldVt,
                              superb);
    //printf("LAPACKE_dgesvd info = %d\n", info);
    assert(info == 0);
    free(Jx);
    free(superb);
#if DEBUG_SVD_MATRICES
    matrix U, s, Vt;
    set_matrix(&U, "U", "left singular vectors", NULL, m, Usz / ldU, &Ux[0]);
    set_matrix(&s, "Σ", "singular values", NULL, ssz, 1, &sx[0]);
    set_matrix(&Vt, "Vᵀ", "right singular vectors", NULL, n, Vtsz / ldVt, &Vtx[0]);
    printf_matrix(&U);
    printf_matrix(&s);
    printf_matrix(&Vt);
#endif
    /* calculate Utb = Uᵀb */
    double Utb[ssz];
    memset(Utb, 0, sizeof(double)*ssz);
    for(int i = 0; i < m; i++) { /* i: b row = Uᵀ col = U row */
        for(int j = 0; j < ssz; j++) { /* j: Uᵀb row = Uᵀ row = U col */
            Utb[j] += Ux[colmaj(i, j, ldU)] * b[i];
        }
    }
#if DEBUG_SVD_MATRICES
    matrix Utbm;
    set_matrix(&Utbm, "Uᵀb", NULL, NULL, ssz, 1, &Utb[0]);
    printf_matrix(&Utbm);
#endif
    /* calculate Utb = ϕΣ⁻¹(Uᵀb) */
    const double hp = mdl->hp;
    for(int j = 0; j < ssz; j++) {
        Utb[j] = (phi(sx[j], hp) * Utb[j]) / (sx[j] + DBL_EPSILON); /* filtered Uᵀb */
    }
#ifdef DEBUG_SVD_MATRICES
    for(int j = 0; j < ssz; j++) {
        printf("ϕ[%d]=%8g Utb[%d]=%8g\n", j, phi(sx[j], hp), j, Utb[j]);
    }
#endif
    /* complete the solution as Δx=V(ϕΣ⁻¹Uᵀb) */
    memset(x, 0, sizeof(double) * n);
    for(int i = 0; i < ssz; i++) { /* i: Uᵀb row = V col = Vᵀ row */
        const double Utbi = Utb[i];
        for(int j = 0; j < n; j++) { /* j: x row = V row = Vᵀ col */
            x[j] += Vtx[colmaj(i, j, ldVt)] * Utbi;
        }
    }
#ifdef DEBUG_SVD_MATRICES
    matrix xm;
    set_matrix(&xm, "x", NULL, NULL, n, 1, &x[0]);
    printf_matrix(&xm);
#endif
    free(Ux);
    free(sx);
    free(Vtx);
    ret = SUCCESS;
    return ret;
}

static int inv_solve_diff(model * const mdl, double const * b1, double const * b2, double * x)
{
    assert(mdl != NULL);
    assert(b1 != NULL);
    assert(b2 != NULL);
    const int m = mdl->n_stimmeas;
    double b[m];
    for(int i = 0; i < m; i++) {
        b[i] = b2[i] - b1[i];
    }
#ifdef DEBUG_SVD_MATRICES
    printf("b = ");
    for(int i = 0; i < m; i++) {
        printf(" %g", b[i]);
    }
    printf("\n");
#endif
    return inv_solve_svd(mdl, &b[0], x);
}

int inv_solve(model * const mdl, double * x)
{
    assert(mdl != NULL);
    assert(x != NULL);
    const int m = mdl->n_stimmeas;
    const int n = mdl->n_params[0];
    for(int frame_idx = 1; frame_idx < mdl->n_data[1]; frame_idx++) {
        int ret = inv_solve_diff(mdl, &(mdl->data[0]), &(mdl->data[frame_idx * m]), &x[(frame_idx - 1) * n]);
        if(!ret) {
            return FAILURE;
        }
        if(mdl->n_params[1] > 1) {
            break;
        }
    }
    // printf_model(mdl);
    return SUCCESS;
}
