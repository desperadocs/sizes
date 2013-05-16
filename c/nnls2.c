/*****************************************************************************

  nnls2.c
  Copyright (c) 2004, Pete R. Jemian

  Algorithm NNLS (Non-negative least-squares)
 
  Given an m by n matrix A, and an m-vector B, computes an n-vector x,
  that solves the least squares problem
      A x = B     subject to x >= 0
  
  Version:
  2004-12-13 Pete Jemian - direct implementation of NNLS algorithm.
       This routine is based on the text in Chapter 23, Section 3 (NNLS) of
       C.L. Lawson and R.J. Hanson, Solving Least Squares Problems,
       SIAM, Philadelphia, PA, 19104, 1995, ISBN 0-89871-356-0.
       NOTE:  This routine is not derived from the NNLS.FOR
       provided by Lawson and Hanson.

     Algorithm NNLS:

     set P (the set of solution vector components unconstrained, x[] >= 0)
     set Z (the set of solution vector components constrained to x[] = 0)

     step  description
      1    assign all elements to set Z (PZ[] = 0) and x[] = 0
      2    Compute the n-Vector w[] = E_tr(f-Ex)
           where E_tr denotes the transpose of matrix E
      3    If set Z is empty (PZ[] != 0 for all elements)
             or w[j] <= 0 for all j in set Z,
             then computation is complete (go to step 12)
      4    Find an index t (from set Z) such that w[t] = MAX(w[j])
      5    Move index t from set Z to set P (PZ[t]=1)
      6    Define an m by n matrix Ep:
               column j of Ep[j][] = {
                     j in set P: E[j][]
                     j in set Z: 0
               }
           Compute the n-vector z as a solution 
           of the least squares problem
              Ep z ~= f
           Note that only the components z[j], 
           with j in set P,
           are determined by this problem.
           Define z[j] = 0 for j in set Z
      7    If z[j] > 0 for all j in set P
           then assign x[] = z[] and go to step 2
      8    Find an index q in set P, such that
           x[q] / (x[q] - z[q]) = MIN( x[j] / (x[j] - z[j]) ) 
           for z[j] <= 0 and j in set P
      9    Assign alpha = x[q] / (x[q] - z[q])
     10    Increment x += alpha * ( z[] - x[] )
     11    move all indices j in set P to set Z
           where x[j] = 0, go to step 6
     12    Comment: computation is complete

  All vectors and matrices used within these routines rely on the 
  Numerical Recipes library functions for dynamic allocation
  and disposal.  They are 1-based, rather than the typical 0-based
  arrays created by direct calls to the C malloc function.
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

   /* ---------- */

#include "recipes.h"
#include "nnls2.h"

/****************************************************************************/

  /* Local function definitions */
static int _end_outer_loop(
    double **E, double *x, double *f, int M, int N,
    double *resid, double *w, int *PZ);

static int _end_inner_loop(
    double **E, int M, int N, double *f, double *x,
    int *PZ, double *z, double *w, int t
);

/* static */
void _least_squares_solution(
    double **a, int M, int N, double *b, double *x
);

/****************************************************************************/

#define SV_CUTOFF   (1.0e-6)


/*****************************************************************************
 
   nnls2(double **E, int M, int N, double *f, double *x, double *w, int *PZ)
     returns
       0   (NNLS_SUCCESS) successful
       1   (NNLS_ITMAX_EXCEEDED) iteration count exceeded
       2   (NNLS_MEMORY_PROBLEM) memory allocation problems
       3   (NNLS_M_LT_N) invalid dimensions (M < N is not allowed)
     definition of parameter list
     E     MxN matrix, input design matrix, where E x = f
     M,N   rows,columns of input design matrix
     f     response M-vector, where E x = f
     x     solution N-vector, where E x = f (this is the result)
     w     dual N-vector
     PZ    index (integer) N-vector, 
           PZ=1 for unconstrained elements (set P), 
           PZ=0 for elements with value constrained to zero (set Z)

 ****************************************************************************/


int nnls2(double **E, int M, int N, double *f, double *x, double *w, int *PZ)
{
  int i, j, k, jp, q, t;
  double sum, alpha, alpha_j, alpha_q, *z, *resid, rnorm;
  int outerLoop, innerLoop;

  if (N > M) {
    return (NNLS_M_LT_N);
  }

  z     = vector(1, N);
  resid = vector(1, M);
  if (!z || !resid) {
    /* dump any memory already allocated */
    free_vector (resid, 1, M);
    free_vector (z,     1, N);
    return (NNLS_MEMORY_PROBLEM);
  }

  /* step 1
     assign all elements of x[] to set Z, assign all x[] = 0
   */
  for (j=1; j<=N; j++) x[j] = PZ[j] = 0;
  /* DEBUG_MARKER;  printf ("assigning all elements to set Z\n"); */

  /* step 2 (& 3)
     Compute the n-Vector w[] = E_tr(f-Ex)
     Typical number of iterations to complete: n/2
     Certainly, 3n is an outrageous number of iterations.
   */
  for (outerLoop = 3*N; outerLoop; outerLoop--) {
    if ( _end_outer_loop(E, x, f, M, N, resid, w, PZ) ) break;
    /* DEBUG_MARKER; /

    /* step 4
       Find an index t (from set Z) such that MAX(w[j])
     */
    t = 0;
    for (j=1; j<=N; j++) {
      if (PZ[j] == 0) {      /* j in set Z */
        if (t == 0)      t = j;  /* 1st time */
        if (w[j] > w[t]) t = j;  /* search ... */
      }
    }

    /* step 5
       Move index t from set Z to set P (PZ[t]=1)
     */
    PZ[t] = 1;
    /* DEBUG_MARKER;  printf ("moving element %d to set P\n", t); */
    
    /* 
       count how many points in set P
         Exit from inner loop must occur at or before jp-1 iterations 
         (where jp is number in set P when inner loop was entered)
     */
    jp = 0;
    for (j=1; j<=N; j++) jp += PZ[j];
    jp = MAX(jp-1, 1);  /* always allow at least once through step 6 */
    for (innerLoop = jp; innerLoop; innerLoop--) {
      /* DEBUG_MARKER; printf ("outer loop = %d  inner loop %d\n", outerLoop, innerLoop); */
      if ( _end_inner_loop(E, M, N, f, x, PZ, z, w, t) == 1) break;
      /* DEBUG_MARKER; */
      t = -1;  /* indicate that step 6 entered from inner loop */

      /* step 8
         Find an index q in set P, such that
             x[q] / (x[q] - z[q]) = MIN( x[j] / (x[j] - z[j]) )
             for z[j] <= 0 and j in set P
       */
      /* 
for (j=1; j<=N; j++)
  printf(
    "<%s,%d> x[%d] = %lg, z[%d] = %lg, PZ[%d] = %d\n", 
    __FILE__, __LINE__, j, x[j], j, z[j], j, PZ[j]
  );
       */
      q = 0;
      for (j=1; j<=N; j++) {
        if (PZ[j] == 1) {                       /*  j in set P */
          if (z[j] <= 0) {                      /* z[j] not positive */
            alpha_j = x[j] / (x[j] - z[j]);
            if (q == 0) {q = j; alpha_q = alpha_j;}             /* 1st time */
            if (alpha_j < alpha_q) {q = j; alpha_q = alpha_j;}  /* search */
/* printf("<%s,%d> alpha_%d = %lg\n", __FILE__, __LINE__, q, alpha_q); */
          }
        }
      }

      /* step 9
         Assign alpha = x[q] / (x[q] - z[q])
       */
      alpha = alpha_q;

      /* step 10
         Increment x += alpha * ( z[] - x[] )
       */
      for (j=1; j<=N; j++) x[j] += alpha * ( z[j] - x[j] );

      /* step 11
         move all indices j in set P to set Z
             where x[j] = 0, go to step 6
       */
      for (j=1; j<=N; j++) {
        if (PZ[j] == 1) {               /*  j in set P */
          /* algorithm says x[j] == 0 */
          if (x[j] <= 0.0) {            /* chap 23, pp. 164-5 */
            PZ[j] = 0;                  /* move j to set Z */
            /* DEBUG_MARKER;  printf ("moving element %d to set Z\n", j); */
          }
        }
      }

    } /* inner loop */
    /* if (innerLoop == 0) return (NNLS_ITMAX_EXCEEDED); */

  } /* outer loop */
  /* if (outerLoop == 0) return (NNLS_ITMAX_EXCEEDED); */

  /* step 12
     Comment: computation is complete
   */

  /*
    return all the memory allocated by this routine
   */
  free_vector (resid, 1, M);
  free_vector (z,     1, N);

  return(0);
}


/****************************************************************************/


static int _end_outer_loop(
  double **E, double *x, double *f, int M, int N, 
  double *resid, double *w, int *PZ)
{
  int i, j, Z_count;
  double sum, rnorm;

  /* step 2
     Compute the n-Vector w[] = E_tr(f-Ex)
   */
  rnorm = 0.0;  /* useful for diagnostics, not reported to calling routine */
  for (i=1; i<=M; i++) {
    sum = 0;
    for (j=1; j<=N; j++) sum += E[i][j] * x[j];
    resid[i] = f[i] - sum;                           /* resid = f - E x */
    rnorm += pow(resid[i], 2);
  }
  rnorm = sqrt(rnorm);  /* magnitude of resid vector */
  /* printf("||E x - f|| = %lg\t", rnorm); DEBUG_MARKER; */
  Z_count = 0;  /* zero this counter */
  for (j=1; j<=N; j++) {
    sum = 0;
    for (i=1; i<=M; i++) sum += E[i][j] * resid[i];
    w[j] = sum;                                      /* w = E_tr resid */
    /* step 3
       Check if computation (outer loop) is complete
       If set Z is empty (PZ[] != 0 for all elements)
             or w[j] <= 0 for all j in set Z
       Here, we check for loop completion 
       and then return the complement as the result.
     */
    if (PZ[j] == 0) {             /* found a member of set Z */
      Z_count++;                  /* count # items in set Z (but w[j]>0) */
      if (w[j] <= 0) Z_count--;   /* count number of non-pos w[j] in set Z */
    }
  }

  return (!Z_count);
}


static int _end_inner_loop(
        double **E, int M, int N, double *f, double *x, 
        int *PZ, double *z, double *w, int t)
{
  int i, j, jp, endLoopTest;
  double **Ep;

  /* step 6
     Define an m by n matrix Ep:
         column j of Ep[j][] = {
               j in set P: E[j][]
               j in set Z: 0
         }

     NOTE: 
       The NNLS code of Lawson and Hanson did this more
       efficiently, without need to copy values to an Ep matrix.
       This was done using features of FORTRAN matrix indexing.
       Here, we suffer the inefficiency but gain code portability.

     Compute the n-vector z as a solution 
     of the least squares problem
        Ep z ~= f
     Note that only the components z[j], 
     with j in set P,
     are determined by this problem.
     Define z[j] = 0 for j in set Z
   */
  Ep = matrix(1,M, 1,N);
  for (i=1; i<=M; i++) {
    for (j=1; j<=N; j++) {
      if (PZ[j] == 1) {
        Ep[i][j] = E[i][j];  /* in set P */
      } else {
        Ep[i][j] = 0;        /* in set Z */
      }
    }
  }
  /* 
     Solve the matrix equation:       Ep z = f
  */
  _least_squares_solution(Ep, M, N, f, z);
  free_matrix (Ep, 1,M, 1,N);  /* don't need this any more */

  /*
     Here, we check for a round-off error situation
     which can occur if we have entered step 6 directly from step 5.
     Our new candidate point (subscript t determined in step 4)
     will theoretically be positive.  
     If it is not, it is surely due to round-off.  
     We then zero that element of the dual vector w[], 
     (without moving the element to set Z),
     and break to the outer loop.
     See page 164 for details.
   */
  if (t>0) {          /* step 6 entered directly from step 5 */
    if (z[t] <= 0) {  /* can only happen due to round-off error */
      w[t] = 0;
      return (1);     /* break from inner loop */
    }
  }

  /* ... continuing ...
      Define z[j] = 0 for j in set Z
   */
  for (j=1; j<=N; j++) {
    if (PZ[j] == 0) {
      z[j] = 0;        /* in set Z */
    }
  }

  /* step 7
     If z[j] > 0 for all j in set P
         then assign x[] = z[] and go to step 2 (outer loop).
   */
  endLoopTest = 1;
  for (j=1; j<=N; j++) {
    if (PZ[j] == 1) {      /* in set P */
      if (z[j] <= 0) {
        endLoopTest = 0;   /* don't end loop */
        break;             /* no need to keep searching */
      }
    }
  }
  if (endLoopTest == 1) {  /* good, did not lose any z[] from set P */
    for (j=1; j<=N; j++) x[j] = z[j];   /* x[] = z[] */
  }

  return (endLoopTest);
}


/* static */
void _least_squares_solution(double **a, int M, int N, double *b, double *x)
{
    /* 
       Solve the matrix equation
         Ep z = f
       Here we use SVD from Numerical Recipes since
       we cannot ensure that M==N.  
       In other languages, such as IgorPro, 
       could use QR or LU decomposition instead 
       if the chosen routine allows the M!=N case.
    */
   double *w, **v, wMax, wMin;
   int j;
   w = vector(1,N);
   v = matrix(1,N, 1,N);

   /*
       Singular Value decomposition a[][] = u[][] w[] v[][] 
       a[][] is changed into u[][] by SVD
    */
   svdcmp(a, M, N, w, v);
   /* find and zero the insignificant singular values */
   wMax = w[1];
   for (j=1; j<=N; j++) if (w[j]>wMax) wMax = w[j];
   wMin = wMax * SV_CUTOFF;
   /* DEBUG_MARKER; SHOW_VECTOR("SVD w", w, 1, N, "%lg"); */
   /* printf ("singular value cutoff = %lg\n", wMin); */
   for (j=1; j<=N; j++) if (w[j]<wMin) w[j] = 0.0;

   /*
       Singular Value backsubstitution solution of a z = b
    */
   svbksb(a, w, v, M, N, b, x);
   /* DEBUG_MARKER; SHOW_VECTOR("SVD bksb x", x, 1, N, "%lg"); */

   free_matrix(v, 1,N, 1,N);
   free_vector(w, 1,N);
}

/************************************************************************/
/************************************************************************/
/************************************************************************/
