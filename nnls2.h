/******************************************************************************

  nnls2.h
  Copyright (c) 2004, Pete R. Jemian

  
  Version:
  2004-12-13 Pete Jemian - direct implementation of NNLS algorithm.
       This routine is based on the text in Chapter 23, Section 3 (NNLS) of
       C.L. Lawson and R.J. Hanson, Solving Least Squares Problems,
       SIAM, Philadelphia, PA, 19104, 1995, ISBN 0-89871-356-0.
       NOTE:  This routine is not derived from the NNLS.FOR
       provided by Lawson and Hanson.

******************************************************************************/
#ifndef _NNLS_H
#define _NNLS_H
/*****************************************************************************/
  /* definitions of provided functions */
int nnls2(double **E, int M, int N, double *f, double *x, double *w, int *PZ);
/*****************************************************************************/
  /* nnls2 return results */
#define NNLS_SUCCESS             0
#define NNLS_ITMAX_EXCEEDED      1
#define NNLS_MEMORY_PROBLEM      2
#define NNLS_M_LT_N              3

/* ------------------------------------------------------------------- */
#define DEBUG_MARKER {printf("<%s,%d>\n", __FILE__, __LINE__); fflush(stdout);}

#define SHOW_VECTOR(_label_, _vec_, _lo_, _hi_, _fmt_)        \
{                                                             \
  int i; char fmt[50];                                        \
  for (i = _lo_; i <= _hi_; i++) {                            \
    sprintf(fmt, "vector %s[%d] = %s\n", _label_, i, _fmt_);  \
    printf(fmt, _vec_[i]);                                    \
  }                                                           \
}
/* ------------------------------------------------------------------- */

#endif
