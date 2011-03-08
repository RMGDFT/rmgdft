/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*

     C interface for BLAS3 calls

*/
#include "md.h"







void dgemm_c(char *transa, char *transb, int *m, int *n, int *k,
             double *alpha, double *a, int *lda, double *b,
             int *ldb, double *beta, double *c, int *ldc)
{
#if CRAY_T3E
    _fcd transa_fcd, transb_fcd;
    double *cptr;
    int idx, nrows, leftover, maxpts, aoffset, coffset, boffset;

    transa_fcd = _cptofcd(transa, 1);
    transb_fcd = _cptofcd(transb, 1);

    SGEMM(transa_fcd, transb_fcd, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#endif
}






void dsymm_c(char *side, char *uplo,
             int *m, int *n, double *alpha, double *a, int *lda,
             double *b, int *ldb, double *beta, double *c, int *ldc)
{
#if CRAY_T3E
    _fcd side_fcd, uplo_fcd;
    double *cptr;
    int idx, nrows, leftover, maxpts, coffset, boffset;

    side_fcd = _cptofcd(side, 1);
    uplo_fcd = _cptofcd(uplo, 1);

    SSYMM(side_fcd, uplo_fcd, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
#endif
}



void dgemv_c(char *trans, int *m, int *n,
             double *alpha, double *a, int *lda, double *x,
             int *incx, double *beta, double *y, int *incy)
{
#if CRAY_T3E
    _fcd trans_fcd;

    trans_fcd = _cptofcd(trans, 1);

    SGEMV(trans_fcd, m, n, alpha, a, lda, x, incx, beta, y, incy);

#else
    dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
#endif
}

void dsymv_c(char *uplo, int *m,
             double *alpha, double *a, int *lda, double *x,
             int *incx, double *beta, double *y, int *incy)
{
#if CRAY_T3E
    _fcd uplo_fcd;

    uplo_fcd = _cptofcd(uplo, 1);

    SSYMV(uplo_fcd, m, alpha, a, lda, x, incx, beta, y, incy);

#else
    dsymv(uplo, m, alpha, a, lda, x, incx, beta, y, incy);
#endif
}


void ssyrk_c(char *uplo, char *trans, int *n, int *k, double *alpha,
             double *a, int *lda, double *beta, double *c, int *ldc)
{
#if CRAY_T3E
    _fcd trans_fcd, uplo_fcd;

    trans_fcd = _cptofcd(trans, 1);
    uplo_fcd = _cptofcd(uplo, 1);

    SSYRK(uplo_fcd, trans_fcd, n, k, alpha, a, lda, beta, c, ldc);

#else

    dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);

#endif
}

void dtrmm_c(char *side, char *uplo, char *trans, char *diag,
             int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb)
{
#if CRAY_T3E
    _fcd side_fcd, diag_fcd, trans_fcd, uplo_fcd;

    side_fcd = _cptofcd(side, 1);
    diag_fcd = _cptofcd(diag, 1);
    trans_fcd = _cptofcd(trans, 1);
    uplo_fcd = _cptofcd(uplo, 1);

    STRMM(side_fcd, uplo_fcd, trans_fcd, diag_fcd, n, k, alpha, a, lda, b, ldb);

#else

    dtrmm(side, uplo, trans, diag, n, k, alpha, a, lda, b, ldb);

#endif
}
